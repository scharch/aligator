#!/usr/bin/env python3

"""
aligator annotate

This script takes assembled genomic contigs from the Ig loci and looks for and
    annotates V/D/J/C genes based on comparison to a previous set of genomic
    annotations from the same or a related species. Requires BLAST and pybedtools.

Usage: annotator.py CONTIGS RSS12 RSS23 LOCUS TARGETGENOME TARGETBED [ options ]

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
   --psfasta pseudo.fa        - Optional fasta file for saving sequences of genes annotated
                                    as pseudogenes.
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
Added check for missing terminal bases by S Olubo & CA Schramm 2024-09-26.
Minor tweaks for D handling by CA Schramm 2024-10-09.
Prevent "V-Region" from being labeled as an exon  by CA Schramm 2024-10-09.
Refactored and rationalized functionality calls by CA Schramm 2024-11-05.
Sorted GFF and fasta output by CASchramm 2024-11-05.
Fixed GFF3 format compatibility by CASchramm 2024-11-06.
Added optional ORF/P fasta output by CASchramm 2024-11-06.
Added checking of gene/exon/RSS boundaries by CASchramm 2024-11-06.
Changed --nonfunctional to --psfasta.

Copyright (c) 2019-2024 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

import sys, os, re, csv, shutil
from docopt import docopt
from pybedtools import BedTool, featurefuncs
from collections import defaultdict, Counter
from operator import itemgetter
import warnings
from Bio import BiopythonWarning
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
		sys.stderr = sys.__stderr__
		sys.exit( "Error: Are you trying to run multiple annotations at the same time?\n\tIf not, please remove existing 'annoTemp' directory and try again.")

	# 0b. Prep output
	sequences = []
	nfseqs = []
	gffRows = []
	mode = 'w'
	if arguments['LOCUS'] == "TRD":
		mode = "a"
	gffhandle = open( arguments['--outgff'], mode )
	gffwriter = csv.writer( gffhandle, delimiter="\t", lineterminator="\n" )
	gffwriter.writerow(["##gff-version 3"])
	with open( arguments['CONTIGS'], 'r' ) as seqIn:
		for seq in SeqIO.parse( seqIn, 'fasta' ):
			gffwriter.writerow([f"##sequence-region {seq.id} 1 {len(seq.seq)}"])


	## THE PIPELINE GETS RUN SEPARATELY FOR EACH GENE IN THE LOCUS
	for gene in evalues:

		#slight kludge
		if gene == "D" and not arguments['LOCUS'] in ['IGH','TRB','TRD']:
			continue

		# 1a. Generate a search database from target genome
		targets     = BedTool( arguments['TARGETBED'] )
		target2  = targets.filter( lambda x: f"{arguments['LOCUS']}{gene}" in x.name and "gene" in x.name ).saveas()
		if targets.file_type == "gff":
			target2  = targets.filter( lambda x: x['functionality']=="F" and f"{arguments['LOCUS']}{gene}" in x.name and x[2].endswith("gene") ).each( featurefuncs.gff2bed ).saveas()

		geneTargets = target2.sequence( fi=arguments['TARGETGENOME'],
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
		if ( (gene=="V" and arguments['LOCUS'] in ["IGH","IGL","TRA","TRB","TRD"]) or
			  (gene=="J" and arguments['LOCUS'] in ["IGH","IGK"]) ):
			rss23 = BedTool( arguments['RSS23'] )
			selectedRSS, statusDict = parseRSS( blastHits, rss23, 'rss23', gene, arguments['LOCUS'] )
		elif (arguments['LOCUS'] in ['TRB','TRD'] and gene =='D'):
			#5' RSS12
			rss12 = BedTool( arguments['RSS12'] )
			selectedRSS, statusDict = parseRSS( blastHits, rss12, 'rss12', gene, arguments['LOCUS'] )

			#3' RS23
			rss23 = BedTool( arguments['RSS23'] )
			tempRSS, tempStatus = parseRSS( blastHits, rss23, 'rss23', gene, arguments['LOCUS'] )
			for hit in tempRSS:
				selectedRSS[ hit ].append(tempRSS[hit][1])
				if tempStatus[ hit ]['type'] == "ORF":
					statusDict[ hit ]['type'] = "ORF"
					statusDict[ hit ]['notes'].append( tempStatus[hit]['notes'] )
		elif gene != "C":
			rss12 = BedTool( arguments['RSS12'] )
			selectedRSS, statusDict = parseRSS( blastHits, rss12, 'rss12', gene, arguments['LOCUS'] )
		else:
			#C gene, no RSS but have to initiate statusDict
			statusDict = { "\t".join(hit[0:6]):{'type':'F','notes':[]} for hit in blastHits }
	
		# 3. Check splice sites and recover exons
		mappedExons, geneBoundaries, statusDict = checkSplice( blastHits, arguments['TARGETBED'], arguments['TARGETGENOME'], 
																						arguments['CONTIGS'], gene, arguments['--blast'], 
																						arguments['--alleledb'], selectedRSS, statusDict )

		# 4. Check functionality
		splicedSequences, statusDict = checkFunctionality( mappedExons, arguments['CONTIGS'], SOURCE_DIR,
																				arguments['LOCUS'], gene, statusDict )

		# 5. Figure out naming
		finalNames, novelG, novelA = assignNames( splicedSequences, arguments['CONTIGS'], targets, arguments['TARGETGENOME'], 
																	arguments['LOCUS'], gene, blast=arguments['--blast'], 
																	codingDB=arguments['--alleledb'] )

		funcNg = 0
		funcNa = 0

		# 6. Collect output
		for b in blastHits:

			stringhit = "\t".join(b[0:6])
			if stringhit not in statusDict:
				print(f"Missing all info for {finalNames.get(stringhit,stringhit)}...?", file=sys.stderr)
				continue
			elif statusDict[ stringhit ][ 'type' ] == "drop":
				print(f"Discarding {finalNames.get(stringhit,stringhit)} due to {'; '.join(statusDict[stringhit]['notes'])}", file=sys.stderr)
				continue
			elif statusDict[ stringhit ][ 'type' ] == "ORF":
				print(f"{finalNames.get(stringhit,stringhit)} marked as ORF due to {'; '.join(statusDict[stringhit]['notes'])}", file=sys.stderr)
			elif statusDict[ stringhit ][ 'type' ] == "P":
				print(f"{finalNames.get(stringhit,stringhit)} marked as pseudogene due to {'; '.join(statusDict[stringhit]['notes'])}", file=sys.stderr)


			# 6a. GFF output
			gType = f"{arguments['LOCUS'][0:2]}_{gene}_gene" 
			eType = "exon"
			if statusDict[stringhit]['type']=="P":
				gType = f"{arguments['LOCUS'][0:2]}_{gene}_pseudogene" 
				eType = "pseudogenic_exon"

			gffRows.append( [ b[0], "ALIGaToR", gType, geneBoundaries[stringhit]['start'], geneBoundaries[stringhit]['stop'], ".", b[5], ".", f"ID={finalNames.get(stringhit, 'NA')};locus={arguments['LOCUS']};functionality={statusDict[stringhit]['type']}" ] )

			for rss in selectedRSS.get( stringhit, [] ):
				if len(rss)==0: continue # D gene with 3' RSS only
				gffRows.append( [ rss[0], "ALIGaToR", rss[6], int(rss[1])+1, int(rss[2]), rss[4], rss[5], ".", f"Parent={finalNames[stringhit]};locus={arguments['LOCUS']};functionality={statusDict[stringhit]['type']}" ] )
			for exon in mappedExons.get( stringhit, [] ):
				exon_name=exon[3].split()
				if exon_name[1] == "V-Region":
					gffRows.append( [ exon[0], "ALIGaToR", "V_region", int(exon[1])+1, int(exon[2]), ".", exon[5], "0", f"Parent={finalNames[stringhit]};locus={arguments['LOCUS']};functionality={statusDict[stringhit]['type']}" ] )
				else:
					gffRows.append( [ exon[0], "ALIGaToR", eType, int(exon[1])+1, int(exon[2]), ".", exon[5], ".", f"Parent={finalNames[stringhit]};exontype={exon_name[1]};locus={arguments['LOCUS']};functionality={statusDict[stringhit]['type']}" ] )

			# 6b. Fasta output - functional coding sequences only
			if statusDict[stringhit]['type'] in ["F", "ORF"]:
				if statusDict[stringhit]['type']=="F" and stringhit in novelG: funcNg += 1
				if statusDict[stringhit]['type']=="F" and stringhit in novelA: funcNa += 1

				#create and save a SeqRecord
				sequences.append( { 'pos':geneBoundaries[stringhit]['start'], 'seq':SeqRecord( Seq(splicedSequences[stringhit]), id=finalNames[stringhit], description="") } )

			elif arguments['--psfasta'] is not None and stringhit in splicedSequences:
				nfseqs.append( { 'pos':geneBoundaries[stringhit]['start'], 'seq':SeqRecord( Seq(splicedSequences[stringhit]), id=finalNames[stringhit], description="") } )

		# 6c. Print some statistics
		num_pseudo = len( [ g for g in statusDict if statusDict[g]['type']=="P" ] )
		num_orf    = len( [ g for g in statusDict if statusDict[g]['type']=="ORF" ] )
		num_func   = len( [ g for g in statusDict if statusDict[g]['type']=="F" ] )
		print( f"{arguments['LOCUS']}{gene}: {num_pseudo+num_orf+num_func} genes found; {num_func} are predicted to be functional.")
		print( f"      {num_pseudo} are labeled as pseudogenes; {num_orf} are labeled as ORFs")
		print( f"      Of {num_func} functional genes, {funcNg} functional novel ASCs were detected and {funcNa} functional new alleles were reported" )
	
	# 6d. Print outputs to file
	gffRows.sort(key=itemgetter(4), reverse=True)
	gffRows.sort(key=itemgetter(3))
	for row in gffRows:
		gffwriter.writerow(row)
	gffhandle.close()
	with open(arguments['--outfasta'], 'w') as fasta_handle:
		sequences.sort( key=itemgetter('pos') )
		SeqIO.write( [ s['seq'] for s in sequences ], fasta_handle, 'fasta' )
	if arguments['--psfasta'] is not None:
		with open(arguments['--psfasta'], 'w') as fasta_handle:
			nfseqs.sort( key=itemgetter('pos') )
			SeqIO.write( [ s['seq'] for s in nfseqs ], fasta_handle, 'fasta' )
		
	# 7. if we were annotating TRD, add it back to previous TRA results
	if arguments['LOCUS'] == 'TRD':
		newAnnotations = parseTRA(arguments['--outgff'])
		newAnnotations.saveas(arguments['--outgff'])

	# 8. clean up
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
