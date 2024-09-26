#!/usr/bin/env python3

"""
aligator annotate

This script takes assembled genomic contigs from the Ig loci and looks for and
    annotates V/D/J/C genes based on comparison to a previous set of genomic
    annotations from the same or a related species. Requires BLAST and pybedtools.

Usage: annotator.py CONTIGS RSS12 RSS23 LOCUS TARGETGENOME TARGETBED [ --alleledb coding.fa --outgff genes.gff --outfasta genes.fa --blast blastn --debug log ]

Options:
   CONTIGS                    - Assembled genomic contigs to be annotated in fasta format.
   RSS12                      - Predicted RSS12 sequences in bed format.
   RSS23                      - Predicted RSS23 sequences in bed format.
   LOCUS                      - Which locus is being annotated: IGH, IGK, IGL, TRA, or TRB.
   TARGETGENOME               - The genomic contigs for comparison in fasta format.
   TARGETBED                  - A bed file containing the annotations for the comparison.
                                    Expects gene names in the format of "IGKV", "IGLJ",
                                    "IGHCA", etc.
   --alleledb coding.fa       - An optional database of known coding alleles other than 
                                    those in the reference genome. Will be used to match
                                    and name discovered alleles.
   --outgff genes.gff         - Where to save the final annotations - will use GFF3 format.
                                    [default: genes.gff]
   --outfasta genes.fa        - Where to save the extracted sequences of the annotated genes.
                                    [default: genes.fa]
   --blast blastn             - Path to the `blastn` executable. [Default: blastn]
   --debug log                - Mangage program messages. Options are `quiet` (no messages),
                                    `log` (save everything to aligator.log file), or `verbose`
                                    (print everything to screen). [default: log]

Created by Chaim A Schramm on 2019-07-16.
Updated and documented by CA Schramm 2019-09-21.
Many updates and tweaks by Simone Olubo, 2022-2024.
Clean up and tweaks by CA Schramm 2024-04-15.
Pulled exon types through by S Olubo 2024-04-15.
Moved functionality checks to separate script by CA Schramm 2024-04-16.
Added debug option by CA Schramm 2024-05-28.
Added check for missing terminal bases by S Olubo & CA Schramm 2024-09-26

Copyright (c) 2019-2024 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

import sys, os, re, csv, shutil
from docopt import docopt
from pybedtools import BedTool
from collections import defaultdict, Counter
from Bio import BiopythonWarning
import warnings
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

SOURCE_DIR = os.environ.get('ALIGATOR_PATH', "")
if not os.path.exists( SOURCE_DIR ):
	find_path = re.match("(.+)/annotator.py", sys.argv[0])
	if find_path:
		SOURCE_DIR = find_path.group(1)
	else:
		sys.exit( "Can't find the code directory, please set ALIGATOR_PATH environmental variable." )
sys.path.append(f"{SOURCE_DIR}/..")
from aligator import *



def main():

	# 0a. Make a temporary working directory
	try:
		os.makedirs( "annoTemp" )
	except FileExistsError:
		sys.exit( "Error: Are you trying to run multiple annotations at the same time?\n\tIf not, please remove existing 'annoTemp' directory and try again.")

	# 0b. Prep output
	sequences = []
	mode = 'w'
	if arguments['LOCUS'] == "TRD":
		mode = "a"
	gffhandle = open( arguments['--outgff'], mode )
	gffwriter = csv.writer( gffhandle, delimiter="\t", lineterminator="\n" )


	## THE PIPELINE GETS RUN SEPARATELY FOR EACH GENE IN THE LOCUS
	for gene in evalues:

		#slight kludge
		if gene == "D" and not arguments['LOCUS'] in ['IGH','TRB','TRD']:
			continue

		# 1a. Generate a search database from target genome
		targets     = BedTool( arguments['TARGETBED'] )
		geneTargets = targets.filter(
									lambda x: f"{arguments['LOCUS']}{gene}" in x.name and "gene" in x.name
									).sequence( fi=arguments['TARGETGENOME'],
												fo=f"annoTemp/targets_{arguments['LOCUS']}{gene}.fa",
					 							name=True, s=True)

		#read in those fasta seqeuences to create a length dictionary for checking below
		targetLengths = dict()
		with open( f"annoTemp/targets_{arguments['LOCUS']}{gene}.fa", 'r' ) as fh:
			for seq in SeqIO.parse( fh, "fasta" ):
				targetLengths[ seq.id ] = len( seq.seq )
		
		# 1b. Run blast
		blast2bed( arguments['--blast'], arguments['CONTIGS'],
							f"annoTemp/targets_{arguments['LOCUS']}{gene}.fa",
							f"annoTemp/rawHits_{arguments['LOCUS']}{gene}.bed",
							evalue=evalues[gene] )

		# 1c. Uniquify the blast hits
		#for constant regions the main concern is assembly gaps that might cause each domain to come up as separate blast hit.
		#    Since genes are much further apart than exons, we can use bedtools merge to account for this
		if gene == "C":
			try:
				blastHits = BedTool( f"annoTemp/rawHits_{arguments['LOCUS']}{gene}.bed" ).merge( s=True, d=5000, c="4,5,6", o="distinct,max,first" ).saveas( f"annoTemp/uniqueHits_{arguments['LOCUS']}{gene}.bed" )
			except:
				print(f"Warning: bedtools merge failed for {arguments['LOCUS']}. Maybe no blast hits were found on this contig?\nSkipping...\n\n", file=sys.stderr)
				continue
		else:
			# For V genes, on the other hand, the most important thing is to get a template of the right family
			#    otherwise, the L-part1 exon might get lost and/or the splice sites might not align properly.
			# To account for this, we use bedtools cluster to group overlapping hits, sort by length and then take
			#    the longest, highest scoring, match for downstream processing.
			# This is really mostly relevant for V genes, but including D and J here for now, as well.
			try:
				BedTool( f"annoTemp/rawHits_{arguments['LOCUS']}{gene}.bed" ).cluster(s=True).saveas( f"annoTemp/clusteredHits_{arguments['LOCUS']}{gene}.bed" )
				#Parse though unique hits output for max score
				with open(f"annoTemp/clusteredHits_{arguments['LOCUS']}{gene}.bed",'r') as input:
					maxBlastHits = list()
					dic = defaultdict(list)
					reader = csv.reader(input, delimiter="\t")
					for row in reader:
						row[8]=int(row[8])					
						x = row[8]
						dic[x].append(row)
				for clust in dic.values():
					sortedHits=sorted( clust,reverse=True, key=lambda r: (int(r[2])-int(r[1]), float(r[4])) )#sort blasthits in descending order first by length then by score
					#check the best hit to make sure it covers the whole query
					#if 1-3 bases are missing from either end, claw them back manually
					#more than that it's probably a pseudogene?
					if int(sortedHits[0][6]) > 1 and int(sortedHits[0][6]) < 4:
						sortedHits[0][1] = str( int(sortedHits[0][1]) - int(sortedHits[0][6]) + 1 )
					if int(sortedHits[0][7]) < targetLengths[ sortedHits[0][3] ] and targetLengths[ sortedHits[0][3] ]-int(sortedHits[0][7]) < 4:
						sortedHits[0][2] = str( int(sortedHits[0][2]) + targetLengths[ sortedHits[0][3] ] - int(sortedHits[0][7]) )

					maxBlastHits.append(sortedHits[0]) #append first hit from each dictionary list to maxBlastHits list 
				with open(f"annoTemp/maxBlastHits_{arguments['LOCUS']}{gene}.bed", 'w') as output:
					writer = csv.writer(output, delimiter="\t")
					for m in sorted( maxBlastHits, key=lambda r: (r[0],r[1]) ):
						writer.writerow(m[0:6])
				blastHits =BedTool(f"annoTemp/maxBlastHits_{arguments['LOCUS']}{gene}.bed")
			except:
				print(f"Warning: `bedtools cluster` failed for {arguments['LOCUS']}{gene}. Maybe no blast hits were found on this contig?\nSkipping...\n\n", file=sys.stderr)
				continue							

		# 2. Match hits with RSS predictions and do some sanity checking
		if (arguments['LOCUS']=="IGH" and (gene=="V" or gene=="J")) or (arguments['LOCUS']=="IGK" and gene=="J") or (arguments['LOCUS']=="IGL" and gene=="V") or (arguments['LOCUS']=="TRA" and gene=="V") or (arguments['LOCUS']=="TRD" and gene=="V") or (arguments['LOCUS']=="TRB" and gene=="V"):
			rss23 = BedTool( arguments['RSS23'] )
			selectedRSS = parseRSS( blastHits, rss23, 'rss23', gene, arguments['LOCUS'] )
		elif (arguments['LOCUS'] in ['TRB','TRD'] and gene =='D'):
			#5' RSS12
			rss12 = BedTool( arguments['RSS12'] )
			selectedRSS = parseRSS( blastHits, rss12, 'rss12', gene, arguments['LOCUS'] )

			#3' RS23
			rss23 = BedTool( arguments['RSS23'] )
			tempRSS = parseRSS( blastHits, rss23, 'rss23', gene, arguments['LOCUS'] )
			for hit in tempRSS:
				selectedRSS[ hit ].append(tempRSS[hit][1])
		else:
			rss12 = BedTool( arguments['RSS12'] )
			selectedRSS = parseRSS( blastHits, rss12, 'rss12', gene, arguments['LOCUS'] )
			
	
		# 3. Check splice sites and recover exons
		mappedExons, geneStatus, spliceNotes = checkSplice( blastHits, arguments['TARGETBED'], arguments['TARGETGENOME'], 
																				arguments['CONTIGS'], gene, arguments['--blast'], 
																				arguments['--alleledb'] )

		# 4. Check functionality
		geneStatus2, splicedSequences, stopCodon, mutatedInvar = checkFunctionality( mappedExons, arguments['CONTIGS'], 
																												SOURCE_DIR, arguments['LOCUS'], gene )

		# 5. Figure out naming
		finalNames, novelG, novelA = assignNames( mappedExons, arguments['CONTIGS'], targets, arguments['TARGETGENOME'], 
																	arguments['LOCUS'], gene, blast=arguments['--blast'], 
																	codingDB=arguments['--alleledb'] )

		funcNg = 0
		funcNa = 0

		# 6. Check if this needs to be labeled as a pseudogene and write GFF output
		#        I don't see an obvious way to mark genes with missing RSS as "ORF" using the Sequence Ontology as required by
		#        GFF3 format, but --on the other hand-- if splice sites are conserved and there are no stop codons, then that
		#        actually seems like pretty good evidence it's a functional gene and most likely a false negative of the RSS
		#        prediction. So I'm going to override the Ramesh et al category definitions on this one.
		for b in blastHits:

			stringhit = "\t".join(b[0:6])

			isPseudo	= False
		
			if stringhit in spliceNotes:
				print(f"{finalNames.get(stringhit, stringhit)}: {spliceNotes[stringhit]}", file=sys.stderr)

			if stringhit in geneStatus:
				print(f"{finalNames.get(stringhit, stringhit)} marked as a pseudogene due to {geneStatus[stringhit]}", file=sys.stderr)
				isPseudo = True
			elif stringhit in geneStatus2:
				print(f"{finalNames.get(stringhit, stringhit)} marked as a {geneStatus2[stringhit]}", file=sys.stderr)
				isPseudo = True

			# 6a. GFF output
			gType = f"{arguments['LOCUS']}_{gene}_gene" 
			eType = "exon"
			if isPseudo:
				gType = f"{arguments['LOCUS']}_{gene}_pseudogene" 
				eType = "pseudogenic_exon"

			gffwriter.writerow( [ b[0], "ALIGaToR", gType, int(b[1])+1, b[2], ".", b[5], ".", f"ID={finalNames.get(stringhit, 'NA')}" ] )
			for rss in selectedRSS.get( stringhit, [] ):
				if len(rss)==0: continue # D gene with 3' RSS only
				gffwriter.writerow( [ rss[0], "ALIGaToR", rss[6], int(rss[1])+1, rss[2], rss[4], rss[5], ".", f"parent={finalNames[stringhit]}" ] )
			for exon in mappedExons.get( stringhit, [] ):
				exon_name=exon[3].split()
				gffwriter.writerow( [ exon[0], "ALIGaToR", f"{exon_name[1]}-{eType}", int(exon[1])+1, exon[2], ".", exon[5], ".", f"parent={finalNames[stringhit]}" ] )

			# 6b. Fasta output - functional coding sequences only
			if not isPseudo:
				if stringhit in novelG: funcNg += 1
				if stringhit in novelA: funcNa += 1

				#create and save a SeqRecord
				sequences.append( SeqRecord( Seq(splicedSequences[stringhit]), id=finalNames[stringhit], description="") )

		num_pseudo = len(set(geneStatus).union(set(geneStatus2)))
		# 6c. Print some statistics
		print( f"{arguments['LOCUS']}{gene}: {len(blastHits)} genes found; {len(selectedRSS)} had predicted RSSs.")
		print( f"      {num_pseudo} are labeled as pseudogene")
		print( f"      Of {len(mappedExons)-num_pseudo} functional genes, {funcNg} functional novel ASCs were detected and {funcNa} functional new alleles were reported" )
	
	# 7. Finish outputs and clean up
	gffhandle.close()
	with open(arguments['--outfasta'], 'w') as fasta_handle:
		SeqIO.write( sequences, fasta_handle, 'fasta' )
		
	# if we were annotating TRD, add it back to previous TRA results
	if arguments['LOCUS'] == 'TRD':
		newAnnotations = parseTRA(arguments['--outgff'])
		newAnnotations.saveas(arguments['--outgff'])

	# clean up
	shutil.rmtree("annoTemp")

	# If this was TRA, we now need to go back and hit TRD, since it is within the TRA locus
	if arguments['LOCUS'] == 'TRA':
		arguments['LOCUS'] = "TRD"
		arguments['--outfasta'] = f"TRGenes_TRD.fa"
		main()



if __name__ == '__main__':

	arguments = docopt(__doc__)

	if arguments['LOCUS'] not in ['IGH','IGK','IGL','TRA','TRB']:
		sys.exit("Valid choices for LOCUS are IGH, IGK, IGL, TRA, or TRB only")

	if arguments['--debug'] not in ['quiet','log','verbose']:
		sys.exit("Valid choices for --debug are quiet, log, or verbose only")
	if arguments['--debug'] == 'quiet':
		sys.stderr = open('/dev/null', 'w')
	elif arguments['--debug'] == 'log':
		sys.stderr = open('aligator.log', 'a')

	evalues = { "V":"1e-20", "D":"1e-10", "J":"1e-10", "C":"1e-100" }

	main()
