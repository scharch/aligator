#!/usr/bin/env python3

"""
assignNames.py

This script assigns names to the genes found by blast.
Note that while naming of novel "alleles" (within 5% of a known allele) should
	be consistent across runs, the naming of novel "genes" (>5% divergence to
	everything known) will only be consistent if the coding sequence is exactly
	the same. Otherwise, post-processing will be required to discover the 
	similarity between hypothetical novel genes IGHV4-kj82leq and IGHV4-2lk493k.

Split out from annotator.py by Chaim A Schramm on 2024-04-17.

Copyright (c) 2024 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

from pybedtools import BedTool
from Bio import SeqIO
from Bio.Seq import Seq
import re, sys, csv
from collections import defaultdict
from hashlib import sha1
from aligator import blastOnly



def assignNames( toName, targets, genomeFile, locus, gene, blast="blastn", codingDB=None ):

	#intialize sequence database
	seqDB = dict()

	#initialize some tracking variables and results database
	novelG = list()
	novelA = list()
	newIDs = dict()


	#start by extracting genomic sequences
	#first get a list of unique gene names in the input bedfile
	target_names = [ t.name for t in targets ]
	gene_ids     = set( re.match( f"{locus}{gene}\S+", n).group() for n in target_names )


	#now extract the exons for each one in turn
	for gene in gene_ids:

		#some kludge to save them since the BedTool generator object produced by `filter` is weird
		geneExons = targets.filter( lambda x: f"{locus}{gene}" in x.name and "exon" in x.name )
		exonList = []
		for e in geneExons:
			exonList.append( e )

		#check strand and get sequence
		splicedSeq = ""
		if exonList[0].strand == "+": #assuming all exons are in the same direction
			for exon in exonList:
				splicedSeq += BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), genomeFile )
		else:
			for exon in reversed(exonList):
				rc = BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), genomeFile )
				splicedSeq += str( Seq(rc).reverse_complement() )

		seqDB[ gene ] = splicedSeq


	#parse DB of coding alleles, if present
	if codingDB is not None:

		cExonOrder = "CH1 H H1 H2 H3 H4 H-CH2 CH2 CH3 CH3-CHS CH4 CH4-CHS CHS M M1 M2 EX1 EX2 EX3 EX4".split()
		cSeqs      = defaultdict( dict )

		with open( codingDB, 'r' ) as dbhandle:
			for seq in SeqIO.parse( dbhandle, 'fasta' ):

				#someone might have downloaded an allele database directly from IMGT, so parse accordingly
				imgtID = re.search( f"{locus}({gene}|[DMAGE])[^\s|]+", seq.description)
				if not imgtID:
					print( f"Cannot parse gene ID {seq.description}, skipping...")
					continue

				if imgtID.group() in seqDB:
					#it's in the genome reference, move on
					continue
				else:
					#IGH and TR* C genes have multiple exons that might need to be put together manually
					if gene == "C" and locus not in ['IGK','IGL']:
						exonID = [ e for e in cExonOrder if re.search(f"\\b{e}\\b",seq.description) ]

						#stupid kludge to account for three types of exons with an internal '\b'
						if "H-CH2" in exonID:
							exonID.remove("H")
							exonID.remove("CH2")
						elif "CH3-CHS" in exonID:
							exonID.remove("CH3")
							exonID.remove("CHS")
						elif "CH4-CHS" in exonID:
							exonID.remove("CH4")
							exonID.remove("CHS")

						if len(exonID) == 1:
							print( f"{imgtID.group()} {exonID[0]}")
							cSeqs[ imgtID.group() ][ exonID[0] ] = str( seq.seq )
						else:
							#multiple or no matches - just assume it is full-length and move on
							seqDB[ imgtID.group() ] = str( seq.seq )
					else:
						seqDB[ imgtID.group() ] = str( seq.seq )

		#put any C exons back together
		for allele in cSeqs:
			seqDB[ allele ] = "".join([ cSeqs[allele].get(e,"") for e in cExonOrder ])



	#now iterate over the hits
	for stringhit, exonList in toName.items():
		
		#Get spliced sequence
		#  TODO - change dictionary structure in checkSplice/checkFunctionality so I don't have to do this again
		splicedSeq = ""
		if exonList[0].strand == "+": #assuming all exons are in the same direction
			for exon in exonList:
				splicedSeq += BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), contigs )
		else:
			for exon in reversed(exonList):
				rc = BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), contigs )
				splicedSeq += str( Seq(rc).reverse_complement() )

		#check for exact match in database
		exactMatches = { k:v if v in splicedSeq or splicedSeq in v for k,v in seqDB.items() }

		#multiple exact matches probably means a duplication of some sort in the reference
		#  print a warning and then choose arbitrarily
		if len(exactMatches) > 1:
			print( f"Warning: {','.join(exactMatches)} all appear to have the same sequence. Arbitrarily using {exactMatches.keys()[0]} for naming purposes...", file=sys.stderr)

		if len(exactMatches) > 0:
			newIDs[ stringhit ] = exactMatches.keys()[0]

			#check is this sequence is a superset of something from a coding allele
			#  if so replace with the longer version
			if len( splicedSeq ) > len( exactMatches[ exactMatches.keys()[0] ] ):
				exactMatches[ exactMatches.keys()[0] ] = splicedSeq

		else:

			#no exact match, will name using seq hash
			seqHash = sha1( splicedSeq.encode() ).hexdigest()[0:4]

			#write sequences to disk and blast for closest match
			fastaString = "\n".join( [ f">{k}\n{v}" for k,v in seqDB ] )
			with open( "annoTemp/nameDB.fa", 'w' ) as fh:
				fh.write( fastaString )
			with open( "annoTemp/nameQuery.fa", 'w' ) fh2:
				fh2.write( f">query\n{splicedSeq}" )
			blastOnly( blast, "annoTemp/nameDB.fa", "annoTemp/nameQuery.fa",
							"annoTemp/nameResult.txt", outformat="6 sseqid pident", 
							minPctID=75, maxTarget=1, maxHSP=1 )

			#read in blast hit
			with open("annoTemp/nameResult.txt", 'r') as handle:

				#should only be one line based on parameters - only need the best match
				try:
					row = next( csv.reader( handle, delimiter="\t") )
				except StopIteration:
					#this shouldn't happen, means blast couldn't find anything with at least 75% id
					print( f"Warning! No sequence in the database has at least 75% identity with {seqHash} ({stringhit})..")
					newIDs[ stringhit ] = f"{locus}{gene}{seqHash}"
					novelG.append(stringhit)

				#parse the gene name of the hit
				#this regex probably still isn't right...
				geneParse = re.match( f"(({locus}{gene}?[DMAGE]?\w*)-?[^*]*)", row[0] )

				if not geneParse:
					#it's probably not great that we have gotten this far without checking this assumption...
					sys.exit(f"ERROR: cannot parse database gene {row[0]} from blast hits for {stringhit}. Please use canonical IMGT naming eg IGHV1-2*02.")

				#check how similar query is to subject
				if row[1] > 0.95:
					newIDs[ stringhit ] = f"{geneParse.match(1)}*{seqHash}"
					novelA.append(stringhit)

				elif gene=="C":
					#special case, since there's no family number
					newIDs[ stringhit ] = f"{geneParse.match(1)}*{seqHash}"
					novelA.append(stringhit)

				else:
					newIDs[ stringhit ] = f"{geneParse.match(0)}-{seqHash}*01"
					novelG.append(stringhit)
					

	return (newIDs, novelG, novelA)
