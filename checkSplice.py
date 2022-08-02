#!/usr/bin/env python3

"""
checkSplice.py

This script checks that the expected splice sites are present in V, J, and C
    genes found by blast.

Created by Chaim A Schramm on 2019-08-06.

Copyright (c) 2019 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

from pybedtools import BedTool
import warnings
from Bio import BiopythonWarning
from Bio import SeqIO
from Bio.Seq import Seq
import re, sys


def gapFinder( exons, seqFile ):
	for e in exons:
		s = BedTool.seq( (e.chrom,e.start,e.stop), seqFile )
		#don't think it pays to worry about gaps of 1 or 2 bases, arbitrarily setting threshold at 10 bases
		if "NNNNNNNNNN" in s:
			return True

	return False


def reindexExons(e, minpos, maxpos):

	if e.strand == "+":
		e.start = e.start - minpos
		e.stop  = e.stop  - minpos
	else:
		temp     = e.start
		e.start  = maxpos - e.stop
		e.stop   = maxpos - temp
		e.strand = "+"

	return e


def mapExons(e, posDict, source ):

	#convert exon coordinates to contig being annotated
	e.chrom  = source.chrom
	e.name   = source.name
	e.strand = source.strand

	if source.strand == "+":
		e.start = posDict[ e.start ]['test'] + source.start
		e.stop  = posDict[ e.stop ]['test'] + source.start
	else:
		temp    = e.start
		e.start = source.stop - posDict[ e.stop ]['test']
		e.stop  = source.stop - posDict[ temp ]['test']

	return e


def checkSplice( hits, bedfile, targetSeq, contigs, gene, blast_exec, codingSeq ):

	from annotateIgLoci import blast2bed, quickAlign, GAPPED_CODON_TABLE

	results = dict()
	reasons = dict()
	notes   = dict()
	targetBed = BedTool(bedfile)

	#iterate through hits
	for h in hits:

		stringhit =	"\t".join(h[0:6])
		pseudo = False

		#check if at least one of the hits has CDS annotations
		names = re.sub(" .*$","",h.name).split(",")
		for n in names:
			exons = targetBed.filter(lambda x: x.name == n + " CDS").saveas()
			if len(exons) > 0:
				if not gapFinder(exons, targetSeq): #avoid things that might cause trouble
					break
		if len(exons) == 0 and codingSeq is not None:
			#no annotation: find one that does with BLAST
			s = BedTool([h]).sequence(fi=contigs,fo="annoTemp/temp.fa",s=True)
			blast2bed(blast_exec, codingSeq, "annoTemp/temp.fa", "annoTemp/temp.bed")
			blastHits = BedTool("annoTemp/temp.bed")
			for b in blastHits:
				correctedName = re.sub( "(-...)-.", "\\1", b.chrom ) #occasional extra letters in coding database names for some reason
				correctedName = re.sub( "\*\d\d", "*01", correctedName ) #change allele designation to 01 to match targetBed
				exons = targetBed.filter( lambda x: x.name == correctedName + " CDS" ).saveas()
				if len(exons) > 0:
					break

		#no matches whatsoever
		if len(exons) == 0:
			reasons[ stringhit ] = "no target CDS found"
			continue

		#now get sequences and align them
		s = BedTool([h]).sequence(fi=contigs,fo="annoTemp/temp.fa",s=True)
		with open("annoTemp/temp.fa",'r') as handle:
			testSeq = SeqIO.read(handle, "fasta")
		r = exons.merge(d=10000,c="4,5,6",o="first").sequence(fi=targetSeq,fo="annoTemp/ref.fa",s=True)
		with open("annoTemp/ref.fa",'r') as handle:
			refSeq = SeqIO.read(handle, "fasta")
		align = quickAlign(refSeq, testSeq)

		#make sure the feature numbering corresponds to the extracted sequence
		#    by rev-comp-ing if necessary and subtracting off the start position
		#    of the first exon. Makes sure we know which is the splice donor and
		#    which is the splice acceptor.
		finalExons = exons.each( reindexExons, exons[0].start, exons[len(exons)-1].stop ).saveas()
		finalExons = finalExons.sort() #if they were negative strand, this puts the flipped exons in order

		#extract and copy CDS/domain descriptions
		#kludgy, but start by making a dictionary of corresponding positions
                    #can use this to try and extend incomplete hits, if need be
		posDict = dict()
		refPos  = -1
		testPos = -1
		for p in range(len(align['ref'])):
			if not align['test'][p] == "-": testPos += 1
			if not align['ref'][p] == "-":
				refPos += 1
				posDict[ refPos ] = { 'align':p, 'test':testPos, 'gap':False }
				if testPos == -1:
					posDict[ refPos ]['test'] = 0 #special case where blast hit is missing 5' end; don't try to claw it back
			if align['test'][p] == "-": posDict[ refPos ][ 'gap' ] = True
		#add closed boundaries:
		if align['ref'][-1] == "-":
			#hit went past last exon, set boundary at testPos of last non gap
			posDict[ refPos+1 ] = { 'test':posDict[refPos]['test']+1, 'gap':False }
		else:
			posDict[ refPos+1 ] = { 'test':testPos+1, 'gap':False }
			if align['test'][-1] == "-": posDict[ refPos+1 ][ 'gap' ] = True

		#check for splice sites for gaps and proper motifs
		# also check if the ends of the alignment are missing
		incomplete5 = False
		incomplete3 = False
		for i in range(len(finalExons)):
			if i > 0:
				acceptor = re.sub("-","",align['test'][ posDict[finalExons[i].start]['align']-10 : posDict[finalExons[i].start]['align'] ]) #if there are gaps here, it's probably bad anyway, but trying for a safety margin
				if not acceptor.endswith("AG"):
					pseudo = True
					reasons[ stringhit ] = "a bad splice acceptor"
					break

			if posDict[ finalExons[i].start ]['gap']:
				if i==0:
					notes[ stringhit ] = f"incomplete 5' end"
					incomplete5 = True
				else:
					notes[ stringhit ] = f"splice acceptor in alignment gap in exon {i+1}"

			if i < len(finalExons)-1:
				donor = re.sub("-","",align['test'][ posDict[finalExons[i].stop]['align'] : posDict[finalExons[i].stop]['align']+10 ]) #if there are gaps here, it's probably bad anyway, but trying for a safety margin
				if not donor.startswith("GT"):
					with warnings.catch_warnings():
						warnings.simplefilter('ignore', BiopythonWarning)
						if not (gene=="C" and Seq(donor).translate().startswith("*")):
							#assume a stop codon at an exon boundary is always CHS...
							pseudo = True
							reasons[ stringhit ] = "a bad splice donor"
							break

			if posDict[ finalExons[i].stop ]['gap']:
				if i==len(finalExons)-1:
					notes[ stringhit ] = f"incomplete 3' end"
					incomplete3 = True
				else:
					notes[ stringhit ] = f"splice donor in alignment gap in exon {i+1}"

		# update coordinates/names to contig being annotated
		mapped = finalExons.each( mapExons, posDict, h ).saveas()
		mapped = mapped.sort() #in case h is negative strand

		#J and C: extract post/pre nt to verify splicing
		#  but don't bother if we've already marked the ends as missing
		if gene == "J" and not incomplete3:
			if h.strand == "+":
				s = BedTool.seq((mapped[0].chrom,mapped[0].stop,mapped[0].stop+2), contigs )
				if not s == "GT":
					pseudo = True
					reasons[ stringhit ] = "a bad splice donor"
			else:
				s = BedTool.seq((mapped[0].chrom,mapped[0].start-2,mapped[0].start), contigs )
				if not s == "AC":
					pseudo = True
					reasons[ stringhit ] = "a bad splice donor"
		elif gene == "C" and not incomplete5:
			if h.strand == "+":
				s = BedTool.seq((mapped[0].chrom,mapped[0].start-2,mapped[0].start), contigs )
				if not s == "AG":
					pseudo = True
					reasons[ stringhit ] = "a bad splice acceptor"
			else:
				s = BedTool.seq((mapped[len(mapped)-1].chrom,mapped[len(mapped)-1].stop,mapped[len(mapped)-1].stop+2), contigs )
				if not s == "CT":
					pseudo = True
					reasons[ stringhit ] = "a bad splice acceptor"

		#check first exon of V gene for start codon
		elif gene == "V":
			if h.strand == "+":
				s = BedTool.seq((mapped[0].chrom, mapped[0].start, mapped[0].start+3), contigs)
				if not s == "ATG":
					pseudo = True
					reasons[ stringhit ] = "an invalid start codon"
			else:
				s = BedTool.seq((mapped[len(mapped)-1].chrom, mapped[len(mapped)-1].stop-3, mapped[len(mapped)-1].stop), contigs)
				if not s == "CAT":
					pseudo = True
					reasons[ stringhit ] = "an invalid start codon"

		results[ stringhit ] = mapped

	return results, reasons, notes
