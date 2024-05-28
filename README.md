# ALIGaToR - Annotator of Loci for IG and T-cell Receptors
A pipeline for annotating genomic contigs from the IG and TR loci. The pipeline includes:
	- A parsing script that extracts gene, exon, and RSS name and corrdinates from reference annotations of choice of closely related species.
	- A prediction script is calls submodule DnaGrep, that predicts RSS sequences based on genomic contigs.
	- Annotator script that uses the extracted reference genome and genomic information to generate a search databse for blast. Blast hits are matched with predicted RSSs. Other scripts are called to check for start and stop codons, and splice sites.

# Getting Started
#clone the aligator repository 
git clone https://github.com/scharch/aligator.git

# Dependencies/Prerequisites 
Python
Beautifulsoup 4.12.3
Muscle
Blast+
BedTools
pyBedTools

# Usage
aligator extract 
aligator predict
aligator annotate
aligator --help