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

Initialize the DnaGrep submodule:

	       git submodule init
	       git submodule update

Set enviromental variable:

		export ALIGATOR_PATH=$(pwd)/aligator

Quick help:

		aligator help


## Example commands and validation results:

Run all validations:

    	bash sample_data/runTests.sh

Compare results to a reference:

		Rscript sample_data/gatorADE <reference annotations> <aligator output>
