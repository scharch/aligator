#!/usr/bin/env python3

"""
annotator.py

This script takes assembled genomic contigs from the Ig loci and looks for and
    annotates V/D/J/C genes based on comparison to a previous set of genomic
    annotations from the same or a related species. Requires BLAST and pybedtools.

Usage: annotator.py CONTIGS RSS12 RSS23 LOCUS TARGETGENOME TARGETBED [ --alleledb coding.fa --outgff annotations.gff --outfasta IgGenes.fa --blast blastn ]

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
   --outgff annotations.gff   - Where to save the final annotations - will use GFF3 format.
                                    [default: annotations.gff]
   --outfasta IgGenes.fa      - Where to save the extracted sequences of the annotated genes.
                                    [default: IgGenes.fa]
   --blast blastn             - Path to the `blastn` executable. [Default: blastn]

Created by Chaim A Schramm on 2019-07-16.
Updated and documented by CA Schramm 2019-09-21.
Many updates and tweaks by Simone Olubo, 2022-2024.
Clean up and tweaks by CA Schramm 2024-04-15.

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
		sys.exit( "Can't find the code directory, please try calling the script using the full absolute path." )
sys.path.append(f"{SOURCE_DIR}/..")
from aligator import *


def checkInvariants( align, gene ):
	#remove gaps in reference (one at a time so we don't get tripped up by shifting indices)
	gap = re.search( "-+", align['ref'] )
	while (gap):
		align['ref']  = align['ref'][0:gap.start()]  + align['ref'][gap.end():]
		align['test'] = align['test'][0:gap.start()] + align['test'][gap.end():]
		gap = re.search( "-+", align['ref'] )

	#look up invariant positions in reference and check
	invars = { "IGH":{ "V":{21:'C', 95:'C'}, "J":{6:'W'} },
			   "IGK":{ "V":{22:'C', 87:'C'}, "J":{3:'F'} },
			   "IGL":{ "V":{21:'C', 88:'C'}, "J":{3:'F'} },
			   "TRA":{ "V":{21:'C', 87:'C'}, "J":{12:'F'} },
			   "TRB":{"V":{22:'C', 90:'C'}, "J":{6:'F'} },
			   "TRD":{ "V":{22:'C', 91:'C'}, "J":{9:'F'}}}


	for pos in invars[ arguments['LOCUS'] ][ gene ]:
		if align[ 'test' ][ pos ] != invars[arguments['LOCUS']][gene][pos]:
			return False

	return True



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
						row[6]=int(row[6])					
						x = row[6]
						dic[x].append(row)
				for clust in dic.values():
					sortedHits=sorted( clust,reverse=True, key=lambda r: (int(r[2])-int(r[1]), float(r[4])) )#sort blasthits in descending order first by length then by score
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
		#       This will print an error message for incomplete hits
		#       but I'm not trying to extend them, because blasting with multiple related genes should hopefully do the trick
		mappedExons, geneStatus, spliceNotes = checkSplice( blastHits, arguments['TARGETBED'], arguments['TARGETGENOME'], arguments['CONTIGS'], gene, arguments['--blast'], arguments['--alleledb'] )

		# 4a. Figure out existing naming
		#     Go through raw databases and track the max known allele number for naming
		alleleMax = defaultdict( int )
		with open( f"annoTemp/targets_{arguments['LOCUS']}{gene}.fa", 'r' ) as dbhandle:
			for seq in SeqIO.parse( dbhandle, 'fasta' ):
				info = re.match("(.+)\*(\d+)", seq.id)
				if info:
					if int(info.group(2))>alleleMax[info.group(1)]:
						alleleMax[ info.group(1) ] = int( info.group(2) )
				else:
					print( f"Unrecognized gene name {seq.id}", file=sys.stderr )
		if arguments['--alleledb'] is not None:
			with open( arguments['--alleledb'], 'r' ) as dbhandle:
				for seq in SeqIO.parse( dbhandle, 'fasta' ):
					info = re.match("(.+?)(?:-.)?\*(\d+)(?:_.+)?", seq.id)
					if info:
						if info.group(1) in alleleMax and int(info.group(2))>alleleMax[info.group(1)]:
							alleleMax[ info.group(1) ] = int( info.group(2) )
					else:
						print( f"Unrecognized gene name {seq.id}", file=sys.stderr )


		novelG = 0
		novelA = 0
		funcNg = 0
		funcNa = 0
		used = defaultdict( int )
		mutatedInvar = 0
		stopCodon = 0

		for b in blastHits:

			localG   = False
			localA   = False
			isPseudo = False

			# 4b. Find closest known sequence and name based on that
			#     Unfortunately, the `merge` operation doesn't preserve the relationship between hit names
			#         scores, and the possibility of the boundaries changing based on overlapping hits means
			#         we can't easily go back to the raw BLAST output. So, get the final sequence from the
			#         post-merge coordinates and BLAST again.

			s = BedTool([b]).sequence(fi=arguments['CONTIGS'],fo="annoTemp/findGene.fa",s=True)
			blastOnly(arguments['--blast'], "annoTemp/findGene.fa", f"annoTemp/targets_{arguments['LOCUS']}{gene}.fa", "annoTemp/blastNames.txt", outformat="6 qseqid pident", minPctID='95' )
			namelist = dict()
			with open("annoTemp/blastNames.txt", 'r') as handle:
				reader = csv.reader( handle, delimiter="\t")
				for row in reader:
					namelist[ row[0] ] = float(row[1])
			if arguments['--alleledb'] is not None:
				blastOnly(arguments['--blast'], "annoTemp/findGene.fa", arguments['--alleledb'], "annoTemp/blastCoding.txt", outformat="6 qseqid pident", minPctID='95' )
				with open("annoTemp/blastCoding.txt", 'r') as handle:
					reader = csv.reader( handle, delimiter="\t")
					for row in reader:
						if row[0] in namelist:
							continue #use full-length identity, rather than coding region only identity if available
						namelist[ row[0] ] = float(row[1])

			if len(namelist) == 0: #if dictionary is empty for namelist then its a novel gene
				novelG += 1
				localG = True
				#print( f"{str(b).strip()} looks like a novel gene, will be named IG{arguments['LOCUS']}{gene}-novel{novelG}" )
				finalName = f"{arguments['LOCUS']}{gene}-novel{novelG}"
			else: #if it isnt novel sort names in dictionary
				byID   = sorted( list(namelist.keys()), key=lambda x: namelist[x], reverse=True )
				geneID = re.sub( "(.+?)(?:-.)?\*\d+(?:_.+)?", "\\1", byID[0] )
				used[ geneID ] += 1
				if namelist[ byID[0] ] == 100: #100 could the combination of blasthit, percentage identify, blast paramerter, allele frequency 
					#exact match
					finalName = byID[0]
					print(f"|\texact match for {finalName}\t|")
				else:
					#treat this as a new allele of a known gene
					# TODO: check if there are multiple possible gene matches and assign to minimize conflicts with other genes???
					novelA += 1
					localA = True
					alleleMax[ geneID ] += 1
					finalName = geneID + f"*{alleleMax[geneID]:02}"

			# 5. Check if this needs to be labeled as a pseudogene and write GFF output
			#        I don't see an obvious way to mark genes with missing RSS as "ORF" using the Sequence Ontology as required by
			#        GFF3 format, but --on the other hand-- if splice sites are conserved and there are no stop codons, then that
			#        actually seems like pretty good evidence it's a functional gene and most likely a false negative of the RSS
			#        prediction. So I'm going to override the Ramesh et al category definitions on this one.
			stringhit = "\t".join(b[0:6])
		
			if stringhit in spliceNotes:
				print(f"{finalName}: {spliceNotes[stringhit]}")

			if stringhit in geneStatus:
				print(f"{finalName} marked as a pseudogene due to {geneStatus[stringhit]}")
				isPseudo = True
			else:

				# 5a. Get spliced sequence to check for stop codons and invariant residues
				splicedSeq = ""
				if b.strand == "+":
					for exon in mappedExons[ stringhit ]:
						splicedSeq += BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), arguments['CONTIGS'] )
				else:
					for exon in reversed(mappedExons[ stringhit ]):
						rc = BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), arguments['CONTIGS'] )
						splicedSeq += str( Seq(rc).reverse_complement() )

				# 5b. V gene: already in frame, translate and align
				if gene == "V":
					with warnings.catch_warnings():
						warnings.simplefilter('ignore', BiopythonWarning)
						splicedAA = Seq( splicedSeq ).translate(table=GAPPED_CODON_TABLE)
					if "*" in splicedAA:
						stopCodon += 1
						print(f"{finalName} marked as a pseudogene due to an internal stop codon")
						isPseudo = True
					else:
						with open( f"{SOURCE_DIR}/{arguments['LOCUS']}{gene}.fa", 'r' ) as refHandle:
							refSeq = SeqIO.read(refHandle, 'fasta')
						align = quickAlign( refSeq, SeqRecord(splicedAA) )
						invar = checkInvariants( align, gene )
						if not invar:
							mutatedInvar += 1
							print(f"{finalName} marked as a pseudogene due to a missing invariant residue")
							isPseudo =True

				# 5c. J gene, align in nt space, then translate
				elif gene == "J":
					with open( f"{SOURCE_DIR}/{arguments['LOCUS']}{gene}.fa", 'r' ) as refHandle:
						refSeq = SeqIO.read(refHandle, 'fasta')
					align = quickAlign( refSeq, SeqRecord(splicedSeq) )
					with warnings.catch_warnings():
						warnings.simplefilter('ignore', BiopythonWarning)
						align['ref']  = str( Seq(align['ref'] ).translate(table=GAPPED_CODON_TABLE) )
						align['test'] = str( Seq(align['test']).translate(table=GAPPED_CODON_TABLE) )
					if "*" in align['test']:
						stopCodon += 1
						print(f"{finalName} marked as a pseudogene due to an internal stop codon")
						isPseudo = True
					else:
						invar = checkInvariants( align, gene )
						if not invar:
							mutatedInvar += 1
							print(f"{finalName} marked as a pseudogene due to a missing invariant residue")
							isPseudo = True

				# 5d. C gene, need to check secreted and membrane-bound CDSs separately.
				#          Push each into frame and translate to check for stop codons.
				elif gene == "C":
					hasStop  = False
					boundary = len(mappedExons[stringhit]) - 2
					if "IGHCA" in finalName: #?
						boundary = len(mappedExons[stringhit]) - 1 #only one M exon for IGA
					toCheck = [ (0,boundary), (boundary, None) ]
					if len(mappedExons[stringhit]) == 1: #For light chains
							checkSeq = BedTool.seq( (mappedExons[stringhit][0][0],int(mappedExons[stringhit][0][1]),int(mappedExons[stringhit][0][2])), arguments['CONTIGS'] ) #this needs to be fixed because exon is not defined outside fo the for loop
							if b.strand=="-":
								checkSeq = str( Seq(checkSeq).reverse_complement() ) 
							checkSeq = "G" + checkSeq
							with warnings.catch_warnings():
								warnings.simplefilter('ignore', BiopythonWarning)
								splicedAA = Seq( checkSeq ).translate(table=GAPPED_CODON_TABLE)
							if "*" in splicedAA:
								stopCodon += 1
								hasStop = True
								print(f"{finalName} marked as a pseudogene due to an internal stop codon")

					else:
						#only loop if there are multiple exons found
						for cds in toCheck:
							checkSeq = "G"
							if b.strand == "+":
								for exon in mappedExons[ stringhit ][ cds[0]:cds[1] ]:
									checkSeq += BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), arguments['CONTIGS'] )
							else:
								for exon in list(reversed(mappedExons[ stringhit ]))[ cds[0]:cds[1] ]:
									rc = BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), arguments['CONTIGS'] )
									checkSeq += str( Seq(rc).reverse_complement() )

							with warnings.catch_warnings():
								warnings.simplefilter('ignore', BiopythonWarning)
								splicedAA = Seq( checkSeq ).translate(table=GAPPED_CODON_TABLE)
							if "*" in splicedAA:
								stopCodon += 1
								hasStop = True
								print(f"{finalName} marked as a pseudogene due to an internal stop codon")
								break #if secreted has stop codon, don't also check M, so it doesn't end up listed twice

							if hasStop:
								isPseudo = True 

			# 5e. GFF output
			gType = f"{arguments['LOCUS']}_{gene}_gene" 
			eType = "exon"
			if isPseudo:
				gType = f"{arguments['LOCUS']}_{gene}_pseudogene" 
				eType = "pseudogenic_exon"

			gffwriter.writerow( [ b[0], "ALIGaToR", gType, int(b[1])+1, b[2], ".", b[5], ".", f"ID={finalName}" ] )
			for rss in selectedRSS.get( stringhit, [] ):
				if len(rss)==0: continue # D gene with 3' RSS only
				gffwriter.writerow( [ rss[0], "ALIGaToR", rss[6], int(rss[1])+1, rss[2], rss[4], rss[5], ".", f"parent={finalName}" ] )
			for exon in mappedExons.get( stringhit, [] ):
				gffwriter.writerow( [ exon[0], "ALIGaToR", eType, int(exon[1])+1, exon[2], ".", exon[5], ".", f"parent={finalName}" ] )

			# 5f. Fasta output - functional coding sequences only
			if not isPseudo:
				if localG: funcNg += 1
				if localA: funcNa += 1

				#create and save a SeqRecord
				sequences.append( SeqRecord( Seq(splicedSeq), id=finalName) )

		# 5h. Print some statistics
		totals = Counter( geneStatus.values() )
		print( f"{arguments['LOCUS']}{gene}: {len(blastHits)} genes found; {len(selectedRSS)} had predicted RSSs.")
		print( f"      {len(geneStatus)+stopCodon+mutatedInvar} are labeled as pseudogenes: {totals['no target CDS found']} without target CDSs, {totals['a bad splice donor']+totals['a bad splice acceptor']} bad splice sites,")
		print( f"          {totals['an invalid start codon']} missing start codons, {stopCodon} internal stop codons, {mutatedInvar} mutated invariants.")
		print( f"      Of {len(mappedExons)-len(geneStatus)-stopCodon-mutatedInvar} functional genes, {funcNg} novel genes were detected and {funcNa} new alleles were reported" )
		print(  "      Genes with more than 2 alleles found: " + ",".join([ g for g in used if used[g]>2 ]) + "\n" )
	
	# 6. Finish outputs and clean up
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

	#log command line
	logCmdLine(sys.argv)


	evalues = { "V":"1e-20", "D":"1e-10", "J":"1e-10", "C":"1e-100" }

	main()
