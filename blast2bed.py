"""
blast2bed.py 

added qstart and qend to outfmt by Simone Olubo & CA Schramm 2024-09-26

Copyright (c) 2019-2024 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.
"""

import subprocess
import csv


def blastOnly(blast, subject, query, outfile, outformat="6 sseqid sstart send qseqid bitscore sstrand qstart qend", minPctID='0', evalue="1e-20", maxTarget="100", maxHSP="100"):

	subprocess.call( [blast, "-subject", subject, "-query", query, "-out", outfile, \
						"-outfmt", outformat, "-gapopen", '5', "-gapextend", '2', "-penalty", '-1', \
						"-reward", '1', "-word_size", '7', "-max_target_seqs", maxTarget, '-max_hsps', maxHSP, \
						'-evalue', evalue, '-perc_identity', minPctID ], stderr=subprocess.DEVNULL)

def blast2bed(blast, subject, query, outfile, evalue="1e-20"):

	blastOnly(blast, subject, query, "annoTemp/blast.temp", evalue=evalue)

	hits = list()
	with open("annoTemp/blast.temp", 'r') as input:
		reader = csv.reader(input, delimiter="\t")

		for row	 in reader:
			#make sure start and stop are in correct order for minus strand hits
			if int(row[2]) < int(row[1]):
				temp = row[1]
				row[1] = row[2]
				row[2] = temp
				row[5] = "-"
			else:
				row[5] = "+"

			#also fix coordinates for bed format
			row[1] = int(row[1]) - 1

			hits.append(row)

	with open(outfile, 'w') as output:
		writer = csv.writer(output, delimiter="\t")
		for h in sorted( hits, key=lambda r: (r[0],r[1]) ):
			writer.writerow(h)
