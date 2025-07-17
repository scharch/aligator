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
- Muscle v3.8 (currently, muscle5 will cause an error)
- Blast+
- BedTools

## Getting Started
After installing the above dependencies, clone the aligator repository:

		git clone https://github.com/scharch/aligator.git

Install required python packages:

		cd aligator
		pip install -r requirements.txt

Initialize the DnaGrep submodule:

		git submodule init
		git submodule update

Set enviromental variable:

		export ALIGATOR_PATH=$(pwd)

Quick help:

		aligator help

Alternatively, a Singularity container for ALIGaToR is available from Zenodo
with the doi **10.5281/zenodo.16004645**

		singularity exec aligator.sif aligator <extract|predict|annotate> ...

## Example commands and validation results:

Validate install:

		bash validateInstall.sh

Re-run analysis from the paper (requires additional dependencies, see script header):

		bash sample_data/runPaperAnalysis.sh

Compare results to a reference (requires R with the packages bedtoolsr and ggalluvial):

		Rscript sample_data/gatorADE <reference annotations> <aligator output>


## Citing:
If you use ALIGaToR, please cite:
Schramm, Olubo, and Douek, "ALIGaToR: A tool for the automated annotation of immunoglobulin and T cell receptor genomic loci" [preprint](https://www.biorxiv.org/content/10.1101/2025.05.02.651960v1)
