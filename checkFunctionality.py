#!/usr/bin/env python3

"""
checkFunctionality.py

This script checks for start/stop codons and expected invariants.

Split out from annotator.py by Chaim A Schramm on 2024-04-16.
Correction for TRxC splicing by CA Scrhamm on 2014-04-17

Copyright (c) 2024 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

from pybedtools import BedTool
import warnings
from Bio import BiopythonWarning
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re, sys


def checkInvariants( align, locus, gene ):
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


	for pos in invars[ locus ][ gene ]:
		if align[ 'test' ][ pos ] != invars[locus][gene][pos]:
			return False

	return True



def checkFunctionality( exonDict, contigs, directory, locus, gene):

	from aligator import quickAlign, GAPPED_CODON_TABLE

	status  = dict()
	seqdict = dict()
	mutatedInvar = 0
	stopCodon = 0

	for stringhit, exonList in exonDict.items():
		
		# 4a. Get spliced sequence
		splicedSeq = ""
		if exonList[0].strand == "+": #assuming all exons are in the same direction
			for exon in exonList:
				splicedSeq += BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), contigs )
		else:
			for exon in reversed(exonList):
				rc = BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), contigs )
				splicedSeq += str( Seq(rc).reverse_complement() )

		seqdict[ stringhit ] = splicedSeq

		# 4b. V gene: already in frame
		# 		Check start codon
		#		Translate and check for stop codons
		#		Align and check for invariants
		if gene == "V":
			if splicedSeq[0:3] != "ATG":
				status[ stringhit ] = "pseudogene due to missing/incorrect start codon"
			else:
				with warnings.catch_warnings():
					warnings.simplefilter('ignore', BiopythonWarning)
					splicedAA = Seq( splicedSeq ).translate(table=GAPPED_CODON_TABLE)
				if "*" in splicedAA:
					stopCodon += 1
					status[ stringhit ] = "pseudogene due to an internal stop codon"
				else:
					with open( f"{directory}/{locus}{gene}.fa", 'r' ) as refHandle:
						refSeq = SeqIO.read(refHandle, 'fasta')
					align = quickAlign( refSeq, SeqRecord(splicedAA) )
					invar = checkInvariants( align, locus, gene )
					if not invar:
						mutatedInvar += 1
						status[ stringhit ] = "pseudogene due to a missing invariant residue"

		# 4c. J gene: no guarantee of frame so align in nt space
		#		Then translate and check for stop codons/invariant
		elif gene == "J":
			with open( f"{directory}/{locus}{gene}.fa", 'r' ) as refHandle:
				refSeq = SeqIO.read(refHandle, 'fasta')
			align = quickAlign( refSeq, SeqRecord(splicedSeq) )
			with warnings.catch_warnings():
				warnings.simplefilter('ignore', BiopythonWarning)
				align['ref']  = str( Seq(align['ref'] ).translate(table=GAPPED_CODON_TABLE) )
				align['test'] = str( Seq(align['test']).translate(table=GAPPED_CODON_TABLE) )
			if "*" in align['test']:
				stopCodon += 1
				status[ stringhit ] = "pseudogene due to an internal stop codon"
			else:
				invar = checkInvariants( align, locus, gene )
				if not invar:
					mutatedInvar += 1
					status[ stringhit ] = "pseudogene due to a missing invariant residue"

		# 4d. C gene, need to check secreted and membrane-bound CDSs separately.
		#          Push each into frame and translate to check for stop codons.
		elif gene == "C":
			boundary = len(exonList) - 2
			if locus in ["TRA","TRB","TRD"]:
				#no secreted isoform for TR
				boundary = len(exonList)
			elif locus=="IGH" and any( [ a in stringhit for a in ["IGHCA","IGHA","IGA","IgA"] ] ):
				boundary = len(exonList) - 1 #only one M exon for IGA

			toCheck = [ (0,boundary), (boundary, None) ]
			if len(exonList) == 1: #For light chains
					checkSeq = BedTool.seq( (exonList[0][0],int(exonList[0][1]),int(exonList[0][2])), contigs ) #this needs to be fixed because exon is not defined outside fo the for loop
					if exonList[0].strand=="-":
						checkSeq = str( Seq(checkSeq).reverse_complement() ) 
					checkSeq = "G" + checkSeq
					with warnings.catch_warnings():
						warnings.simplefilter('ignore', BiopythonWarning)
						splicedAA = Seq( checkSeq ).translate(table=GAPPED_CODON_TABLE)
					if "*" in splicedAA:
						stopCodon += 1
						status[ stringhit ] = "pseudogene due to an internal stop codon"

			else:
				#only loop if there are multiple exons found
				for cds in toCheck:
					checkSeq = "G"
					if exonList[0].strand == "+":
						for exon in exonList[ cds[0]:cds[1] ]:
							checkSeq += BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), contigs )
					else:
						for exon in list(reversed(exonList))[ cds[0]:cds[1] ]:
							rc = BedTool.seq( (exon[0],int(exon[1]),int(exon[2])), contigs )
							checkSeq += str( Seq(rc).reverse_complement() )

					with warnings.catch_warnings():
						warnings.simplefilter('ignore', BiopythonWarning)
						splicedAA = Seq( checkSeq ).translate(table=GAPPED_CODON_TABLE)
					if "*" in splicedAA:
						stopCodon += 1
						status[ stringhit ] = "pseudogene due to an internal stop codon"
						break #if secreted has stop codon, don't also check M, so it doesn't end up listed twice

	return status, seqdict, stopCodon, mutatedInvar
