import requests
import sys
import re, csv
from bs4 import BeautifulSoup

# Command line arguments for URL and reference name
URL = sys.argv[1]
IMGTrefName = sys.argv[2]

# Fetch and parse page content
page = requests.get(URL)
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
RSList = []
VSlist = []

# Store rows to write to the output file
rows = []

for x in soup.select('div',attrs={'class':'features'}):
	for tr in x.select('tr'):
		row = [re.sub("\n",":",i.text.strip()) for i in tr]
		
		#check for functional and pseudo V-genes
		if inVgene or inJgene or inCgene and geneType in ["V-GENE", "J-GENE","C-GENE"]:
			if row[3] == "IMGT_allele":
				geneName = row[5]
			elif row[3] == "pseudo":
				pseudoList.append(geneName)
				inVgene = False
				inJgene = False
				inCgene = False
		#check for functional and pseudo V-genes, D-genes, and J-genes, keep track of gene names
		#If coordinates are "complement", the strand should be "-"
		#Subtract one from start and end coordinates
		#output IMGT reference name, start and end coordinates, and gene name
		if inGeneVDJ:
			if row[3] == "pseudo":
				if geneName not in pseudoList:
					pseudoList.append(geneName)
				inGeneVDJ = False
			elif row[3] == "IMGT_allele":
				geneName = row[5]
				if geneName in geneNameList:
					inGeneVDJ = False
				else:
					geneNameList.append(geneName)
					if geneName in pseudoList:
						inGeneVDJ = False
					else:
						if re.match("complement",splitStartEnd):
							strand = "-"
							splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
						start, end = map(int, splitStartEnd.split(".."))
						rows.append([IMGTrefName, str(start - 1), str(end), geneName + " gene", "0", strand])
			#fetch V,D,and J CDS
			elif row[1] in ["L-PART1","V-EXON","D-REGION","J-REGION"]:
				splitStartEnd = row[5]
				if geneName in pseudoList:
					inGeneVDJ = False
				else:
					if re.match("complement",splitStartEnd):
						strand = "-"
						splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
					start, end = map(int, splitStartEnd.split(".."))
					if row[1] == "L-PART1":
						rows.append([IMGTrefName, str(start - 1), str(end), geneName + " L-part1-exon", "0", strand])
					if row[1] == "V-EXON":
						rows.append([IMGTrefName, str(start - 1), str(end), geneName + " V-exon", "0", strand])
					if row[1] == "D-REGION":
						rows.append([IMGTrefName, str(start - 1), str(end), geneName + " D-exon", "0", strand])
					if row[1] == "J-REGION":
						rows.append([IMGTrefName, str(start - 1), str(end), geneName + " J-exon", "0", strand])
		#check for RS regions of V-genes, D-genes, and J-genes, keep track of gene names
		#If coordinates are "complement", strand should be "-"
		#Subtract one from start and end coordinates
		#output IMGT reference name, start and end coordinates, and gene name			
		elif re.match("(V|5'D|3'D|J)-RS", row[1]) and row[5] not in RSList and geneName not in VSlist:
			#fetch V-gene and J-gene RSS
			if row[1] == "V-RS" or row[1] == "J-RS":
				VSlist.append(geneName)
			RSList.append(row[5])
			splitStartEnd = row[5]
			if geneName in pseudoList:
				inGeneVDJ = False
			else:
				if re.match("complement",splitStartEnd):
					strand = "-"
					splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
				start, end = map(int, splitStartEnd.split(".."))
				rows.append([IMGTrefName, str(start - 1), str(end), geneName + " RS", "0", strand])
		
		#Check for C genes, and their CDS, keep track of gene names
		#filter out pseudo C genes
		elif inGeneCunit and geneType == "C":
			if geneName in pseudoList or row[1] == "3'UTR":
				inGeneCunit = False
			elif row[3] == "IMGT_allele":
				geneName = row[5]
				if geneName in geneNameList:
					inGeneCunit = False
				else:
					geneNameList.append(geneName)
					Cgene = geneName.replace("IGH","IGHC")
					if re.match("complement",splitStartEnd):
						strand = "-"
						splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
					start, end = map(int, splitStartEnd.split(".."))
					rows.append([IMGTrefName, str(start - 1), str(end), Cgene + " gene", "0", strand])
				#fetch C-gene CDS
			if re.match("CL",row[1]) or re.match("CH|H",row[1]) or re.match("M[12]|M(?!ISC)",row[1]) or re.match("EX([1234])",row[1]):
				splitStartEnd = row[5]
				if re.match("CH([1234])-CH([2S])|H-CH([1-4])",row[1]):
					CexonList.append(geneName)
				if re.match("CH|H|CL|EX([1234])",row[1]) and geneName not in CexonList or re.match("M[12]|M(?!ISC)",row[1]) or re.match("CH([1234])-CH([2S])|H-CH([1-4])",row[1]):
					if re.match("complement",splitStartEnd):
						strand = "-"
						splitStartEnd = re.sub(r"complement|\(|\)", "", splitStartEnd)
					start, end = map(int, splitStartEnd.split(".."))
					rows.append([IMGTrefName, str(start - 1), str(end), f"{Cgene} {row[1]}-exon", "0", strand])

		else:
			#If V,D,J gene names have not been recorded, flag gene type, record start and end coordinates, set strand to default "+"
			if re.match("(?:L-)?([VDJ])-GENE-UNIT", row[1]):
				splitStartEnd = row[5]
				strand="+"
				inGeneVDJ = True
				if row[1]=="L-V-GENE-UNIT":
					geneType="V"
				elif row[1]=="D-GENE-UNIT":
					geneType="D"
				elif row[1]=="J-GENE-UNIT":
					geneType="J"
			if re.match("V-GENE",row[1]):
				inVgene = True
				geneType = "V-GENE"
			if re.match("J-GENE",row[1]):
				inJgene = True
				geneType = "J-GENE"
			if re.match("C-GENE",row[1]):
				inCgene = True
				geneType = "C-GENE"
			if row[1] == "C-GENE-UNIT":
				inGeneCunit = True
				splitStartEnd = row[5]
				geneType = "C"

# Sort the rows
rows.sort(key=lambda x: (x[0], int(x[1])))

# Write the sorted rows to the output file
with open(f"{IMGTrefName}.bed", 'w') as output:
    writer = csv.writer(output, delimiter="\t")
    for row in rows:
        writer.writerow(row)




		
