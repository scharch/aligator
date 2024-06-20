<p align="center">
<image src=https://github.com/scharch/aligator/assets/6708960/c6acd3d9-d082-4b0b-9f09-c99c7f8f651a>
</p>

# ALIGaToR - Annotator of Loci for IG and T-cell Receptors
A pipeline for annotating genomic contigs from the IG and TR loci. The pipeline includes:
- Extract: A parsing script that extracts gene, exon, and RSS name and corrdinates from reference annotations of choice of closely related species.
- Predict: A prediction script calls submodule DnaGrep, that predicts RSS sequences based on genomic contigs.
- Annotate: Annotator script that uses the extracted reference genome and genomic information to generate a search databse for blast. Blast hits are matched with predicted RSSs. Other scripts are called to check for start and stop codons, and splice sites.

## Dependencies/Prerequisites 
- Python 3.6 or greater
- Muscle
- Blast+
- BedTools

## Getting Started
Clone the aligator repository:

		git clone https://github.com/scharch/aligator.git

Install required python packages:

		pip install -r aligator/requirements.txt

Set enviromental variable:

		export ALIGATOR_PATH=$(pwd)/aligator

Quick help:

		`aligator help`


## Vignette annotating MF989451 from Ramesh et al Frontiers Immunology 2017:
Data is in `aligator/sample_data`.

First, get reference genome from IMGT:

		#Download BK063715 fasta file from https://imgt.org/ligmdb/view.action?format=FASTA&id=BK063715
		#Then create bedfile with reference annotations
		aligator extract https://imgt.org/ligmdb/view.action?id=BK063715 BK063715

Find possible RSS motifs in the target contig. For MF989451, the output should look the same as `sample_data/MF989451.rss12_pred.bed` and `sample_data/MF989451.rss23_pred.bed`:
   
		aligator predict $ALIGATOR_PATH/sample_data/MF989451.fa MF989451

Finally, annotate the target contig. For MF989451, the actual annotations provided by Ramesh et al are included as `sample_data/MF989451.ground_truth.bed`:
   
		aligator annotate $ALIGATOR_PATH/sample_data/MF989451.fa MF989451.RSS12.bed MF989451.RSS23.bed IGH BK063715.fasta BK063715.bed

