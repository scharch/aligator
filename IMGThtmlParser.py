#!/usr/bin/env python3

"""
aligator extract

This script parses reference annotations from the IMGT website to produce a
    bedfile that can be used with `aligator annoate`

Usage: IMGThtmlParser.py IMGTREFNAME OUTPUTBED [--pseudo] [--fasta=FASTA]

Arguments:
	IMGTREFNAME     - IMGT Reference name
	OUTPUTBED       - output bed file

Options:	
	--pseudo        - optional flag to include pseudogenes
	--fasta=FASTA   - optional file to include fasta sequence

Created by S Olubo 2023.
Add docopt documentation by CA Schramm 2024-05-28.
Included optional pseudogene flag S Olubo 2025-02-24
Included optional fasta sequence S Olubo 2025-02-24
Changed command line arguments to include output bed file S Olubo 2025-02-24
Overrode IMGT fasta def line to prevent mismatches with bed file
     and a few other minor fixes by CA Schramm 2025-02-25.

Copyright (c) 2023-2025 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

from docopt import docopt
import requests
import sys
import re, csv
from bs4 import BeautifulSoup


def main():
	# Fetch and parse page content
	url = f"https://imgt.org/ligmdb/view.action?id={arguments['IMGTREFNAME']}"
	page = requests.get(url)
	soup = BeautifulSoup(page.content, "html.parser")

	# Initialize variables and lists for tracking gene features
	inGeneVDJ = False
	inGeneCunit = False
	inVgene = False
	inJgene = False
	inCgene = False
	geneName = ""
	start = ""
	end = ""
	geneType = ""
	geneNameList = []
	CexonList = []
	pseudoList = []
	ORFList= []
	RSList = []
	VSlist = []

	# Store rows to write to the output file
	rows = []

	for x in soup.select('div',attrs={'class':'features'}):
		for tr in x.select('tr'):
			row = [re.sub("\n",":",i.text.strip()) for i in tr]
			
			#skip the second chunk of the extraction
			if row[1]=="Database": 
				break

			#check for functional, pseudo, and ORF V-genes
			if inVgene or inJgene or inCgene and geneType in ["V-GENE", "J-GENE","C-GENE"]:
				if row[3] == "IMGT_allele":
					geneName = row[5]
				#if pseudogene flag is true then add line for pseudogene, if there is no pseudogene flag, skip pseudogenes
				elif row[3] == "pseudo":
					if not arguments['--pseudo']:
						inVgene = False
						inJgene = False
						inCgene = False
					else:
						if re.match("complement",splitStartEnd):
							strand = "-"
							splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
						splitStartEnd = re.sub( "[<>]", "", splitStartEnd )
						start, end = map(int, splitStartEnd.split(".."))
						geneName.replace('(I)','').replace('(II)','').replace('(III)','')
						if geneName not in pseudoList:
							rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " pseudogene", "0", strand])
					pseudoList.append(geneName)
				elif row[3] == "ORF":
					ORFList.append(geneName)
			#check for functional and pseudo V-genes, D-genes, and J-genes, keep track of gene names
			#If coordinates are "complement", the strand should be "-"
			#Subtract one from start and end coordinates
			#output IMGT reference name, start and end coordinates, and gene name
			if inGeneVDJ:
				if row[3] == "pseudo":
					if not arguments['--pseudo']:
						inGeneVDJ = False
					if geneName not in pseudoList:
						pseudoList.append(geneName)
				elif row[3] == "IMGT_allele":
					geneName = row[5]
					if geneName not in geneNameList:
						geneNameList.append(geneName)
						if geneName in pseudoList and not arguments['--pseudo']:
							inGeneVDJ = False
						elif geneType!="V-GENE": #only do this for L-V-GENE-UNIT, not V-GENE
							if re.match("complement",splitStartEnd):
								strand = "-"
								splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
							splitStartEnd = re.sub( "[<>]", "", splitStartEnd )
							start, end = map(int, splitStartEnd.split(".."))
							if geneName in ORFList:
								rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " ORF gene", "0", strand])
							elif geneName not in pseudoList:
								rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " gene", "0", strand])
				#fetch V,D,and J CDS
				elif row[1] in ["L-PART1","V-EXON","V-REGION","D-REGION","J-REGION"]:
					splitStartEnd = row[5]
					if geneName in pseudoList:
						inGeneVDJ = False
					else:
						if re.match("complement",splitStartEnd):
							strand = "-"
							splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
						splitStartEnd = re.sub( "[<>]", "", splitStartEnd )
						start, end = map(int, splitStartEnd.split(".."))
						if row[1] == "L-PART1":
							rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " L-part1-exon", "0", strand])
						if row[1] == "V-EXON":
							rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " V-exon", "0", strand])
						if row[1] == "V-REGION":
							rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " V-region", "0", strand])
						if row[1] == "D-REGION":
							rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " D-exon", "0", strand])
						if row[1] == "J-REGION":
							rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " J-exon", "0", strand])
							inGeneVDJ = False
			#check for RS regions of V-genes, D-genes, and J-genes, keep track of gene names
			#If coordinates are "complement", strand should be "-"
			#Subtract one from start and end coordinates
			#output IMGT reference name, start and end coordinates, and gene name			
			if re.match("(V|5'D|3'D|J)-RS", row[1]) and row[5] not in RSList and geneName not in VSlist:
				#fetch V-gene and J-gene RSS
				if row[1] == "V-RS" or row[1] == "J-RS":
					VSlist.append(geneName)
				if row[1] == "V-RS" or row[1] == "3'D-RS":
					inGeneVDJ = False
				RSList.append(row[5])
				splitStartEnd = row[5]
				if geneName in pseudoList:
					inGeneVDJ = False
				else:
					if re.match("complement",splitStartEnd):
						strand = "-"
						splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
					splitStartEnd = re.sub( "[<>]", "", splitStartEnd )
					start, end = map(int, splitStartEnd.split(".."))
					rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " RS", "0", strand])
			
			#Check for C genes, and their CDS, keep track of gene names
			#filter out pseudo C genes
			if inGeneCunit and geneType == "C":
				if geneName in pseudoList or row[1] == "3'UTR":
					inGeneCunit = False
				elif row[3] == "IMGT_allele":
					geneName = row[5].replace("IGH","IGHC")
					if geneName in geneNameList:
						inGeneCunit = False
					else:
						geneNameList.append(geneName)
						if re.match("complement",splitStartEnd):
							strand = "-"
							splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
						splitStartEnd = re.sub( "[<>]", "", splitStartEnd )
						start, end = map(int, splitStartEnd.split(".."))
						rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " gene", "0", strand])
					#fetch C-gene CDS
				if re.match("CL",row[1]) or re.match("CH|H",row[1]) or re.match("M[12]|M(?!ISC)",row[1]) or re.match(r"^EX([1234])$", row[1]):
					splitStartEnd = row[5]
					if re.match("CH([1234])-CH([2S])|H-CH([1-4])",row[1]):
						CexonList.append(geneName)
					if re.match("CH|H|CL|EX([1234])",row[1]) and geneName not in CexonList or re.match("M[12]|M(?!ISC)",row[1]) or re.match("CH([1234])-CH([2S])|H-CH([1-4])",row[1]):
						if re.match("complement",splitStartEnd):
							strand = "-"
							splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
						splitStartEnd = re.sub( "[<>]", "", splitStartEnd )
						start, end = map(int, splitStartEnd.split(".."))
						rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), f"{geneName} {row[1]}-exon", "0", strand])

			if re.match("(?:L-)?([VDJ])-GENE(-UNIT)?", row[1]):

				if row[1]=="D-GENE" or row[1]=="J-GENE":
					continue

				splitStartEnd = row[5]
				strand="+"

				if row[1]=="L-V-GENE-UNIT":
					geneType="V"
					if inGeneVDJ: #already have the gene name from V-GENE
						if re.match("complement",splitStartEnd):
							strand = "-"
							splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
						splitStartEnd = re.sub( "[<>]", "", splitStartEnd )
						start, end = map(int, splitStartEnd.split(".."))
						if geneName in ORFList:
							rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " ORF gene", "0", strand])
						elif geneName not in pseudoList:
							rows.append([arguments['IMGTREFNAME'], str(start - 1), str(end), geneName + " gene", "0", strand])

				elif row[1]=="D-GENE-UNIT":
					geneType="D"
				elif row[1]=="J-GENE-UNIT":
					geneType="J"

				inGeneVDJ = True
			if re.match("V-GENE",row[1]):
				inVgene = True
				geneType = "V-GENE"
				strand="+"
				splitStartEnd = row[5]
			if re.match("J-GENE",row[1]):
				inJgene = True
				geneType = "J-GENE"
				strand="+"
				splitStartEnd = row[5]
			if re.match("C-GENE",row[1]):
				inCgene = True
				geneType = "C-GENE"
			if row[1] == "C-GENE-UNIT":
				inGeneCunit = True
				splitStartEnd = row[5]
				geneType = "C"

		#don't process additional `div`s after break
		else:
			continue
		break

	# Sort the rows
	rows.sort(key=lambda x: (x[0], int(x[1])))

	# Write the sorted rows to the output file
	with open(arguments['OUTPUTBED'], 'w', newline='') as output:
		writer = csv.writer(output, delimiter="\t", lineterminator="\n")
		for row in rows:
			writer.writerow(row)

	if arguments['--fasta']:
		fastaUrl = f"https://imgt.org/ligmdb/view.action?format=FASTA&id={arguments['IMGTREFNAME']}"
		fastaPage = requests.get(fastaUrl)
		fastaSoup = BeautifulSoup(fastaPage.content, "html.parser")
		with open(arguments['--fasta'], "w") as fastaFile:
			for y in fastaSoup.select('pre'):
				for row in y:
					fixName = re.sub( ">.+?\n", f">{arguments['IMGTREFNAME']}\n", row.text.strip() )
					fastaFile.write(fixName)


if __name__ == '__main__':

	arguments = docopt(__doc__)

	main()

