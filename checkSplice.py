#!/usr/bin/env python3

"""
checkSplice.py

This script checks that the expected splice sites are present in V, J, and C
    genes found by blast.

Created by Chaim A Schramm on 2019-08-06.
Switched from CDS to exon and pulled exon labels through by Simone Olubo 2024-04-15.
Some quick clean up by CA Schramm 2024-04-16.
Moved check for start codon to new script by CA Schramm 2024-04-16.
Removed lines 89-110 and pulled full gene through for the alignment by Simone Olubo & CA Schramm 2024-09-26
Added debugging for V-region by S Olubo & CA Schramm 2024-10-01
Fixed start position numbering for D/J genes and removed splice checking for
    D genes by CA Schramm 2024-10-09.
Fixed exon type regex by CA Schramm 2024-10-09.
Refactored and rationalized functionality calls by CA Schramm 2024-11-05.
Fixed GFF3 format compatibility by CASchramm 2024-11-06.
Added checking of gene/exon/RSS boundaries by CASchramm 2024-11-06.

Copyright (c) 2019-2024 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

from pybedtools import BedTool, featurefuncs
import warnings
from Bio import BiopythonWarning
from Bio import SeqIO
from Bio.Seq import Seq
import re, sys

#returns if there is a gap or not
def gapFinder( exons, seqFile ):
	for e in exons:
		s = BedTool.seq( (e.chrom,e.start,e.stop), seqFile ).upper()
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

	exon_type = re.sub("-exon", "", e.name)

	#convert exon coordinates to contig being annotated
	e.chrom  = source.chrom
	e.name   = exon_type
	e.strand = source.strand

	if source.strand == "+":
		e.start = posDict[ e.start ]['test'] + source.start
		e.stop  = posDict[ e.stop - 1 ]['test'] + source.start + 1
	else:
		temp    = e.start
		e.start = source.stop - posDict[ e.stop - 1 ]['test'] - 1
		e.stop  = source.stop - posDict[ temp ]['test']

	return e


def exons2bed( gff ):
	bed = featurefuncs.gff2bed( gff )

	if gff[2] == "V_Region":
		bed.name = bed.name + " V-Region"
	elif gff[2] == "exon":
		bed.name = bed.name + f" {gff['exontype']}-exon"

	return bed


def checkSplice( hits, bedfile, targetSeq, contigs, gene, locus, blast_exec, codingSeq, rss, status ):

	from aligator import blast2bed, quickAlign, GAPPED_CODON_TABLE

	results = dict()
	geneBoundaries = dict()
	targetBed = BedTool(bedfile)
	#iterate through hits
	for h in hits:

		stringhit =	"\t".join(h[0:6])

		#check if at least one of the hits has annotated exons
		names = re.sub(r"[:\s].*$","",h.name).split(",")
		exons = targetBed.filter(lambda x: names[0] in x.name and ("exon" in x.name or "V-Region" in x.name)).saveas()
		fullGene = targetBed.filter(lambda x: names[0] in x.name and "gene" in x.name).saveas()
		if targetBed.file_type == "gff":
			exons = targetBed.filter(lambda x: x.name==names[0] and x[2] in ["exon","V-Region"]).each( exons2bed ).saveas()
			fullGene = targetBed.filter(lambda x: x.name==names[0] and x[2].endswith("gene")).each( featurefuncs.gff2bed ).saveas()	

		if len(fullGene) == 0:
			print( f"Warning: could not extract gene {names[0]} for {stringhit}, skipping...", file=sys.stderr )
			continue

		#now get sequences and align them
		s = BedTool([h]).sequence(fi=contigs,fo="annoTemp/temp.fa",s=True)
		with open("annoTemp/temp.fa",'r') as handle:
			testSeq = SeqIO.read(handle, "fasta")
		r = fullGene.sequence(fi=targetSeq,fo="annoTemp/ref.fa",s=True)
		try:
			with open("annoTemp/ref.fa",'r') as handle:
				refSeq = SeqIO.read(handle, "fasta")
		except ValueError:
			sys.stderr = sys.__stderr__
			sys.exit( f"Error: multiple reference sequences with identical names: {names[0]}")
		align = quickAlign(refSeq, testSeq)

		#make sure the feature numbering corresponds to the extracted sequence
		#    by rev-comp-ing if necessary and subtracting off the start position
		#    of the first exon. Makes sure we know which is the splice donor and
		#    which is the splice acceptor.
		finalExons = exons.each( reindexExons, fullGene[0].start, fullGene[0].stop ).saveas()
		finalExons = finalExons.sort() #if they were negative strand, this puts the flipped exons in order

		#extract and copy exon/domain descriptions
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
		incomplete = False
		if gene != "D":
			for i in range(len(finalExons)):

				if "V-Region" in finalExons[i].name:
					continue

				#note: muscle returns alignments in caps, even if inputs were lower case
				#    so as long as we pull from `align`, don't have to worry about that

				if i > 0: #C acceptor handled below
					acceptor = re.sub("-","",align['test'][ posDict[finalExons[i].start]['align']-10 : posDict[finalExons[i].start]['align'] ]) #if there are gaps here, it's probably bad anyway, but trying for a safety margin
					if not acceptor.endswith("AG"):# or acceptor.endswith("AC"):
						status[ stringhit ][ 'type' ] = "ORF"
						status[ stringhit ][ 'notes' ].append( f"noncanonical splice acceptor {acceptor[-2:]}" )

				if i==0:
					if posDict[ finalExons[i].start ]['gap']:
						#is it just a few bases we can claw back?
						test = [ posDict[ base ]['gap'] for base in reversed(range(finalExons[i].start,finalExons[i].start+4)) ]
						clawback = 4 - test.index(True)
						#really should check that we don't go off the begining of the contig...
						if clawback > 3:
							#not fixable
							status[ stringhit ] = { 'type':'drop', 'notes':["incomplete 5' end"] }
							incomplete = True
							break
				else:
					if posDict[ finalExons[i].start - 1 ]['gap']:
						status[ stringhit ] = { 'type':'drop', 'notes':[f"splice acceptor in alignment gap in exon {i+1}"] }
						incomplete = True
						break

				if i < len(finalExons)-1 and not "V-exon" in finalExons[i].name: #J donor handled below, V-exon has V-Region after but is not a splice donor
					donor = re.sub("-","",align['test'][ posDict[finalExons[i].stop]['align'] : posDict[finalExons[i].stop]['align']+10 ]) #if there are gaps here, it's probably bad anyway, but trying for a safety margin
					if not donor.startswith("GT"):# or donor.startswith("GC"):
						with warnings.catch_warnings():
							warnings.simplefilter('ignore', BiopythonWarning)
							if not (gene=="C" and Seq(donor).translate().startswith("*")):
								#assume a stop codon at an exon boundary is always CHS...
								status[ stringhit ][ 'type' ] = "ORF"
								status[ stringhit ][ 'notes' ].append( f"noncanonical splice donor {donor[0:2]}" )

				if i==len(finalExons)-1 or "V-exon" in finalExons[i].name:
					if posDict[ finalExons[i].stop -1 ]['gap']:
						#is it just a few bases we can claw back?
						test = [ posDict[ base ]['gap'] for base in range(finalExons[i].stop-4,finalExons[i].stop) ]
						clawback = 4 - test.index(True)
						#really should check that we don't go off the end of the contig...
						if clawback > 3:
							#not fixable
							status[ stringhit ] = { 'type':'drop', 'notes':["incomplete 3' end"] }
							incomplete = True
							break
				else:
					if posDict[ finalExons[i].stop ]['gap']:
						status[ stringhit ] = { 'type':'drop', 'notes':[f"splice donor in alignment gap in exon {i+1}"] }
						incomplete = True
						break

		if incomplete:
			continue

		# update coordinates/names to contig being annotated
		mapped = finalExons.each( mapExons, posDict, h ).saveas()
		mapped = mapped.sort() #in case h is negative strand

		#remove exons that are empty
		mapped = [ mapped[i] for i in range(len(mapped)) if mapped[i].stop - mapped[i].start > 1 ]
		if len(mapped)==0:
			status[ stringhit ] = { 'type':'drop', 'notes':["no exons found"] }
			continue

		#check/fix gene and exon boundaries
		if stringhit in rss:
			if gene == "V":
				if h.strand == "+":
					for i in range(len(mapped)):
						if "V-Region" in mapped[i].name or "V-exon" in mapped[i].name:
							if mapped[i].stop != int(rss[stringhit][0][1]):
								print( f"Adjusting {mapped[i].name} stop position from {mapped[i].stop} to {rss[stringhit][0][1]} to account for RSS detection", file=sys.stderr )
								mapped[i].stop = int(rss[stringhit][0][1])
					geneBoundaries[ stringhit ] = { 'start':min([e.start for e in mapped])+1, 'stop':int(rss[stringhit][0][2]) }
				else:
					for i in range(len(mapped)):
						if "V-Region" in mapped[i].name or "V-exon" in mapped[i].name:
							if mapped[i].start != int(rss[stringhit][0][2]):
								print( f"Adjusting {mapped[i].name} start position from {mapped[i].start} to {rss[stringhit][0][2]} to account for RSS detection", file=sys.stderr )
								mapped[i].start = int(rss[stringhit][0][2])
					geneBoundaries[ stringhit ] = { 'start':int(rss[stringhit][0][1])+1, 'stop':max([e.stop for e in mapped]) }

			elif gene == "D":
				geneStart = mapped[0].start+1
				geneStop  = mapped[0].stop
				if h.strand == "+":
					if len(rss[stringhit][0]) > 0: #make sure 5' RSS was detected
						if mapped[0].start != int(rss[stringhit][0][2]):
							print( f"Adjusting {mapped[0].name} start position from {mapped[0].start} to {rss[stringhit][0][2]} to account for RSS detection", file=sys.stderr )
							mapped[0].start = int(rss[stringhit][0][2])
						geneStart = int(rss[stringhit][0][1])+1
					if len(rss[stringhit]) > 1: #3' RSS detected
						if mapped[0].stop != int(rss[stringhit][1][1]):
							print( f"Adjusting {mapped[0].name} stop position from {mapped[0].stop} to {rss[stringhit][1][1]} to account for RSS detection", file=sys.stderr )
							mapped[0].stop = int(rss[stringhit][1][1])
						geneStop = int(rss[stringhit][1][2])
				else:
					if len(rss[stringhit][0]) > 0: #make sure 5' RSS was detected
						if mapped[0].stop != int(rss[stringhit][0][1]):
							print( f"Adjusting {mapped[0].name} stop position from {mapped[0].stop} to {rss[stringhit][0][1]} to account for RSS detection", file=sys.stderr )
							mapped[0].stop = int(rss[stringhit][0][1])
						geneStop = int(rss[stringhit][0][2])
					if len(rss[stringhit]) > 1: #3' RSS detected
						if mapped[0].start != int(rss[stringhit][1][2]):
							print( f"Adjusting {mapped[0].name} start position from {mapped[0].start} to {rss[stringhit][1][2]} to account for RSS detection", file=sys.stderr )
							mapped[0].start = int(rss[stringhit][1][2])
						geneStart = int(rss[stringhit][1][1])+1
				geneBoundaries[ stringhit ] = { 'start':geneStart, 'stop':geneStop }

			elif gene == "J":
				if h.strand == "+":
					if mapped[0].start != int(rss[stringhit][0][2]):
						print( f"Adjusting {mapped[0].name} start position from {mapped[0].start} to {rss[stringhit][0][2]} to account for RSS detection", file=sys.stderr )
						mapped[0].start = int(rss[stringhit][0][2])
					geneBoundaries[ stringhit ] = { 'start':int(rss[stringhit][0][1])+1, 'stop':mapped[0].stop }
				else:
					if mapped[0].stop != int(rss[stringhit][0][1]):
						print( f"Adjusting {mapped[0].name} stop position from {mapped[0].stop} to {rss[stringhit][0][1]} to account for RSS detection", file=sys.stderr )
						mapped[0].stop = int(rss[stringhit][0][1])
					geneBoundaries[ stringhit ] = { 'start':mapped[0].start+1, 'stop':int(rss[stringhit][0][2]) }

		else: #C/no rss
			#just use min/max of exon(s)
			geneBoundaries[ stringhit ] = { 'start':min([e.start for e in mapped])+1, 'stop':max([e.stop for e in mapped]) }


		#J and C: extract post/pre nt to verify splicing
		#    **here, though, we are pulling from the genome, so need to correct for possible
		#      case mismatch (ask me how I know...)
		if gene == "J":
			if h.strand == "+":
				s = BedTool.seq((mapped[0].chrom,mapped[0].stop,mapped[0].stop+2), contigs )
				if not s.upper() == "GT":
					status[ stringhit ][ 'type' ] = "ORF"
					status[ stringhit ][ 'notes' ].append( "noncanonical J splice donor" )
			else:
				s = BedTool.seq((mapped[0].chrom,mapped[0].start-2,mapped[0].start), contigs )
				if not s.upper() == "AC":
					status[ stringhit ][ 'type' ] = "ORF"
					status[ stringhit ][ 'notes' ].append( "noncanonical J splice donor" )
		elif gene == "C":
			if h.strand == "+":
				s = BedTool.seq((mapped[0].chrom,mapped[0].start-2,mapped[0].start), contigs )
				if not s.upper() == "AG":
					status[ stringhit ][ 'type' ] = "ORF"
					status[ stringhit ][ 'notes' ].append( "noncanonical C splice acceptor" )
			else:
				s = BedTool.seq((mapped[len(mapped)-1].chrom,mapped[len(mapped)-1].stop,mapped[len(mapped)-1].stop+2), contigs )
				if not s.upper() == "CT":
					status[ stringhit ][ 'type' ] = "ORF"
					status[ stringhit ][ 'notes' ].append( "noncanonical C splice acceptor" )

		results[ stringhit ] = mapped

	return results, geneBoundaries, status
