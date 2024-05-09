#!/bin/bash

# predictRSS.sh
#
# This script uses the DnaGrep application (Merelli et al NAR 2010) to predict the presence
# 	of Recombination Signal Sequences (RSS) using the Recombination Information Content (RIC)
# 	developed by Cowell et al (Genome Biol 2002). It is part of the ALIGaToR immunogenetic
# 	annotation suite.
#
# Usage: predictRSS.sh </path/to/aligator> <genome_to_predict.fasta> <output_prefix>
#
# Output will be in BED format with names OUTPUT_PREFIX.RSS12.bed and OUTPUT_PREFIX.RSS23.bed
#
# Created by Chaim A Schramm on 2024-04-15.
#
# Copyright (c) 2024 Vaccine Research Center, National Institutes of Health, USA.
# All rights reserved.


ALIGATOR_PATH=$1
INPUT_FASTA=$2
OUTPUT_PREFIX=$3

if [[ ! -f $ALIGATOR_PATH/DnaGrep/DnaGrep ]]; then
	git submodule update --init
	#going to assume executable in distro works
fi

$ALIGATOR_PATH/DnaGrep/DnaGrep model -q -m 0 $INPUT_FASTA $ALIGATOR_PATH/DnaGrep/fasta/HS23RSS.fasta $(cat $ALIGATOR_PATH/DnaGrep/model/HS23.model) ca | \
		awk '$6>-58.45 {print $1,$2-1,$3-1,$4,$6,$5}' OFS="\t" > $OUTPUT_PREFIX.RSS23.bed

$ALIGATOR_PATH/DnaGrep/DnaGrep model -q -m 0 $INPUT_FASTA $ALIGATOR_PATH/DnaGrep/fasta/HS12RSS.fasta $(cat $ALIGATOR_PATH/DnaGrep/model/HS12.model) ca | \
		awk '$6>-38.81 {print $1,$2-1,$3-1,$4,$6,$5}' OFS="\t" > $OUTPUT_PREFIX.RSS12.bed
