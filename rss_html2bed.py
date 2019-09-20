#!/usr/bin/env python3

"""
rss_html2bed.py

This script takes html output of RSS predictions from RSSsite
    (https://www.itb.cnr.it/rss/analyze.html) and converts it to BED format
	for use with annotateIgLoci.py.

Usage: rss_html2bed.py INFOLDER OUTFOLDER [ COORDFIXES ]

Options:
    INFOLDER    - Directory with saved html output. There should be at least
                      two files (RSS12 and RSS23) for each locus (IGH, IGK,
					  and IGL); if there are multiple files for a locus/length
				      combination, they will be treated as annotations of
				      different contigs and concatenated.
	OUTFOLDER   - Directory in which to save the bedfile output.
	COORDFIXES  - Optional tab delimited file to fix coordinates if a contig
	                  had to be split for RSSsite to be able to handle it. First
					  column is the name of the split contig (which must be in
					  the output file name) and the second column is the offset
					  to add back to the prediction coordinates.


Created by Chaim A Schramm on 2019-07-24.
Updated and documented by CA Schramm 2019-09-21.

Copyright (c) 2019 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""


from docopt import docopt
import pandas as pd
import glob, re, sys, csv


def main():
	for locus in ("IGH","IGK","IGL"):
		for rss in ("RSS12","RSS23"):

			data = []

			for file in glob.glob(f"{arguments['INFOLDER']}/{locus}*{rss}.html"):

				print(f"Reading {file}... ", end='', flush=True)
				tables = pd.read_html(file,header=1,flavor="html5lib")
				print("Done.")

				df = tables[1]
				df.columns = ['name','startend','sequence','strand','score','passfail']

				pf  = df[ df['passfail'] == "PASS" ]
				pf2 = pf.drop(columns=['startend','passfail'])
				pf2['start'],pf2['end'] = pf.startend.str.split(' - ',1).str
				pf2 = pf2[['name','start','end','sequence','score','strand']]

				pf2.start = pd.to_numeric(pf2.start) - 1

				for c in contigsToFix:
					if re.search(c,file):
						pf2.start = pd.to_numeric(pf2.start) + contigsToFix[c]
						pf2.end   = pd.to_numeric(pf2.end)   + contigsToFix[c]

				data.append( pf2 )

			all_data = pd.concat( data )
			all_data.sort_values(by=["name", "start"], inplace=True)

			all_data.to_csv(f"{arguments['OUTFOLDER']}/{locus}_{rss}.bed",index=False,sep="\t",header=False)


if __name__ == '__main__':

	arguments = docopt(__doc__)

	contigsToFix = dict()
	if arguments['COORDFIXES'] is not None:
		with open( arguments['COORDFIXES'], 'r') as handle:
			reader = csv.reader( handle, delimiter="\t")
			for row in reader:
				contigsToFix[ row[0] ] = int( row[1] )

	main()
