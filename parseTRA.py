#!/usr/bin/env python3

from pybedtools import BedTool
import csv, re
import pandas as pd

def parseTRA(gffFile):
	#open the TRD annotations output file
	with open(f"annoTemp/TRDalleles.gff", 'w') as output:
		writer = csv.writer(output, delimiter="\t")
		with open(gffFile, 'r') as input:
			reader = csv.reader(input, delimiter="\t")
			for row in reader:
				if re.match("#",row[0]):
					continue
				#write row with TRD annotations to TRDalleles.gff
				if re.match("TRD",row[8]):
					writer.writerow(row)
	TRDalleles = BedTool("annoTemp/TRDalleles.gff").sort()
	allAlleles = BedTool(gffFile).sort()

	#use BedTool to cluster the alleles
	BedTool(allAlleles).cluster(s=True).saveas(f"annoTemp/clusteredAlleles.gff")

	columnNames =["contig","source","feature","start","end","score","strand","frame","alleleName","cluster"]
	df = pd.read_csv("annoTemp/clusteredAlleles.gff",delimiter="\t",on_bad_lines="skip",header=None,comment="#")
	df.columns = columnNames

	grouped = df.groupby('cluster')
	df2 = pd.DataFrame()
	#create window from functional TRDV to TRDC, and remove fake TRA psuedogenes within window and fake TRD pseudos outside!
	for group_name, group_data in grouped:
		#find if there are any TRA and TRD genes within the same group
		traFound = any(group_data['alleleName'].apply(lambda x: bool(re.search(r'TRA', x))))
		trdFound = any(group_data['alleleName'].apply(lambda x: bool(re.search(r'TRD', x))))
		trdPseudoFound = False
		traPseudoFound = False
		trdFunctionalFound = False
		traFunctionalFound = False
		#look for specific gene features
		for feature in group_data['feature']:
			if re.search('TRD_V_pseudogene', feature):
				trdPseudoFound = True
			elif re.search('TRA_V_pseudogene', feature):
				traPseudoFound = True
			elif re.search('TRA_V_gene', feature):
				traFunctionalFound = True
			elif re.search('TRD_V_gene', feature):
				trdFunctionalFound = True
		#if both TRA and TRD alleles are found within the same group, determine whether to remove or keep any TRA do TRD functional/pseudogenes
		if traFound and trdFound:
			if traFunctionalFound and trdFunctionalFound:
				print("both TRA and TRD functional genes found, remove TRD")
				group_data = group_data[~group_data['alleleName'].str.contains('TRD')]
				df2 = pd.concat([df2,group_data])
			if traPseudoFound and trdPseudoFound:
				print("both TRA and TRD pseudogenes found, manually inspect these alleles:")
				print(group_data[['feature','start','end','alleleName']])
				df2 = pd.concat([df2,group_data])
			if traFunctionalFound and trdPseudoFound:
				group_data = group_data[~group_data['alleleName'].str.contains('TRD')]
				df2 = pd.concat([df2,group_data])
			if trdFunctionalFound and traPseudoFound:
				group_data = group_data[~group_data['alleleName'].str.contains('TRA')]
				df2 = pd.concat([df2,group_data])
		else:
			df2 = pd.concat([df2,group_data])


	df2 = df2.drop('cluster',axis=1)

	#sort and return the final annotations
	df2.sort_values( ['start', 'end'], ascending=[True, False], inplace=True )

	#copy over header
	header = []
	with open(gffFile, 'r') as readHeader:
		header.append( next(readHeader) )
		header.append( next(readHeader) )
	with open(gffFile, 'w') as printFile:
		for comment in header:
			printFile.write(comment)
		df2.to_csv(printFile, sep='\t', index=False, header=False)
