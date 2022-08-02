
from collections import defaultdict

def parseRSS( blastHits, rssPreds, gene ):

	overlap = blastHits.window( rssPreds, w=10 )
	possible = defaultdict( dict )

	for predicted in overlap:

		strpred = "\t".join(predicted[0:6])

		if gene=="J" or gene =="D":
			#looking at 5p end

			if '5p' not in possible[ strpred ]:
				possible[ strpred ][ '5p' ] = []

			if predicted[5] == "+":
				#plus strand so RSS is minus; also check boundaries
				if predicted[11] == "-" and int(predicted[8])+10 >= int(predicted[1]) and int(predicted[7])-10 <= int(predicted[1]):
						#looks good - save it
						possible[ strpred ][ '5p' ].append( predicted[6:12] )
			else:
				#minus strand, reverse everything
				if predicted[11] == "+" and int(predicted[8])+10 >= int(predicted[2]) and int(predicted[7])-10 <= int(predicted[2]):
						#looks good - save it
						possible[ strpred ][ '5p' ].append( predicted[6:12] )

		if gene=="V" or gene=="D":
			#looking at 3p end

			if '3p' not in possible[ strpred ]:
				possible[ strpred ][ '3p' ] = []

			if predicted[5] == "+":
				#plus strand so RSS is also plus; also check boundaries
				if predicted[11] == "+" and int(predicted[8])+10 >= int(predicted[2]) and int(predicted[7])-10 <= int(predicted[2]):
						#looks good - save it
						possible[ strpred ][ '3p' ].append( predicted[6:12] )
			else:
				#minus strand, reverse everything
				if predicted[11] == "-" and int(predicted[8])+10 >= int(predicted[1]) and int(predicted[7])-10 <= int(predicted[1]):
						#looks good - save it
						possible[ strpred ][ '3p' ].append( predicted[6:12] )

	#now go back through the original list, selecting the one with the best score if multiple were predicted
	results = dict()
	for hit in blastHits:
		strhit = "\t".join(hit[0:6])

		try:
			if gene == "V":
				if len(possible[ strhit ][ '3p' ]) > 0:
					#sort by score
					best = sorted( possible[ strhit ][ '3p' ], key=lambda bed: float(bed[4]), reverse=True )[0]
					results[ strhit ] = [best + [ "V_gene_recombination_feature" ]]

			if gene == "J":
				if len(possible[ strhit ][ '5p' ]) > 0:
					#sort by score
					best = sorted( possible[ strhit ][ '5p' ], key=lambda bed: float(bed[4]), reverse=True )[0]
					results[ strhit ] = [best + [ "J_gene_recombination_feature" ]]

			if gene == "D":
				#check both ends
				if len(possible[ strhit ][ '5p' ]) > 0:
					#sort by score
					best5 = sorted( possible[ strhit ][ '5p' ], key=lambda bed: float(bed[4]), reverse=True )[0]
					results[ strhit ] = [ best5+['five_prime_D_recombination_signal_sequence'] ]
				else:
					results[ strhit ] = [ [] ]
				if len(possible[ strhit ][ '3p' ]) > 0:
					#sort by score
					best3 = sorted( possible[ strhit ][ '3p' ], key=lambda bed: float(bed[4]), reverse=True )[0]
					results[ strhit ].append( best3+['three_prime_D_recombination_signal_sequence'] )

		except KeyError:
			#this just means there was no predicted RSSs for this hit, don't worry about it
			pass

	return results
