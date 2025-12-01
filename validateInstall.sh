#!/bin/bash

python $ALIGATOR_PATH/aligator extract BK063715 validate.bed

if cmp --silent "$ALIGATOR_PATH/sample_data/IGH/BK063715.bed" "validate.bed"; then
  echo "aligator extract validated successfully"
else
  echo "Possible issue with aligator extract, likely due to IMGT interference. Please try manual verification with the --local option"
fi


python $ALIGATOR_PATH/aligator predict $ALIGATOR_PATH/sample_data/IGH/MF989451.fa validate

if cmp --silent "$ALIGATOR_PATH/sample_data/IGH/MF989451.RSS12.bed" "validate.RSS12.bed" & cmp --silent "$ALIGATOR_PATH/sample_data/IGH/MF989451.RSS23.bed" "validate.RSS23.bed"; then
  echo "aligator predict validated successfully"
else
  echo "possible issue with aligator predict"
fi


python $ALIGATOR_PATH/aligator annotate $ALIGATOR_PATH/sample_data/IGH/MF989451.fa \
										$ALIGATOR_PATH/sample_data/IGH/MF989451.RSS12.bed \
										$ALIGATOR_PATH/sample_data/IGH/MF989451.RSS23.bed \
										IGH $ALIGATOR_PATH/sample_data/IGH/BK063715.fasta \
										$ALIGATOR_PATH/sample_data/IGH/BK063715.bed \
										--outgff validate.gff --debug quiet

#naming is sometimes a little arbtirary/wobbly so just look at the first 8 columns
cut -f1-8 validate.gff > validate_new.gff
cut -f1-8 $ALIGATOR_PATH/sample_data/IGH/predictedGenes.gff > validate_old.gff

if cmp --silent "validate_old.gff" "validate_new.gff"; then
  echo "aligator annotate validated successfully"
else
  echo "possible issue with aligator annotate"
fi

rm validate.bed validate.RSS12.bed validate.RSS23.bed validate.gff validate_old.gff validate_new.gff
