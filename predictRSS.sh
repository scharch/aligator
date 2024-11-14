#!/bin/bash

# predictRSS.sh
#
# This script uses the DnaGrep application (Merelli et al NAR 2010) to predict the presence
# 	of Recombination Signal Sequences (RSS) using the Recombination Information Content (RIC)
# 	developed by Cowell et al (Genome Biol 2002). It is part of the ALIGaToR immunogenetic
# 	annotation suite.
#
# Usage: predictRSS.sh </path/to/aligator> <genome_to_predict.fasta> <output_prefix> [<species> <12thresh> <23thresh>]
#
# Output will be in BED format with names OUTPUT_PREFIX.RSS12.bed and OUTPUT_PREFIX.RSS23.bed
# Currently only HS and MM are recognized for <species>, with HS as the default.
# <12thresh> and <23thresh> are the minimum RIC for counting a predicted RSS as real. Defaults
#     are the human thresholds from the RSSsite paper, 12thresh=-38.81 and 23thresh=-58.85
#
# Created by Chaim A Schramm on 2024-04-15.
# Tweaked RSS23 threshold slightly (seems to be a typo in the RSSsite 
#     online docs) by CA Schramm 2024-10-09.
# Added option to use mouse model and custom thresholds by CASchramm 2024-11-14
#
# Copyright (c) 2024 Vaccine Research Center, National Institutes of Health, USA.
# All rights reserved.


ALIGATOR_PATH=$1
INPUT_FASTA=$2
OUTPUT_PREFIX=$3
SPECIES=${4:-HS}
MIN12=${5:--38.81}
MIN23=${6:--58.85}

if [[ "$SPECIES" != "HS" && "$SPECIES" != "MM" ]]; then
	echo "Only 'HS' and 'MM' are recognized for \`--species\`..."
	exit 1
fi

if [[ ! -f $ALIGATOR_PATH/DnaGrep/DnaGrep ]]; then
	git submodule update --init
	#going to assume executable in distro works
fi

$ALIGATOR_PATH/DnaGrep/DnaGrep model -q -m 0 $INPUT_FASTA $ALIGATOR_PATH/DnaGrep/fasta/${SPECIES}23RSS.fasta $(cat $ALIGATOR_PATH/DnaGrep/model/${SPECIES}23.model) ca | \
		awk "\$6>$MIN23 {print \$1,\$2-1,\$3-1,\$4,\$6,\$5}" OFS="\t" > $OUTPUT_PREFIX.RSS23.bed

$ALIGATOR_PATH/DnaGrep/DnaGrep model -q -m 0 $INPUT_FASTA $ALIGATOR_PATH/DnaGrep/fasta/${SPECIES}12RSS.fasta $(cat $ALIGATOR_PATH/DnaGrep/model/${SPECIES}12.model) ca | \
		awk "\$6>$MIN12 {print \$1,\$2-1,\$3-1,\$4,\$6,\$5}" OFS="\t" > $OUTPUT_PREFIX.RSS12.bed
