# ALIGaToR - Annotator of Loci for IG and T-cell Receptors
A pipeline for annotating genomic contigs from the IG and TR loci. The pipeline includes:
- Extract: A parsing script that extracts gene, exon, and RSS name and corrdinates from reference annotations of choice of closely related species.
- Predict: A prediction script calls submodule DnaGrep, that predicts RSS sequences based on genomic contigs.
- Annotate: Annotator script that uses the extracted reference genome and genomic information to generate a search databse for blast. Blast hits are matched with predicted RSSs. Other scripts are called to check for start and stop codons, and splice sites.

## Getting Started
Clone the aligator repository 
git clone https://github.com/scharch/aligator.git

## Dependencies/Prerequisites 
- Python
- Beautifulsoup 4.12.3
- Muscle
- Blast+
- pyBedTools

## Usage
aligator --help
### Example
		 #Download BK063715 fasta file from IMGT.org
		 #extract IGH annotations from IMGT's rheMac10 
		 aligator extract https://imgt.org/ligmdb/view.action?id=BK063715 BK063715
   
		 #predict RSS for MF989451 and compare to sample data 
		 aligator predict /sample_data /sample_data/MF989451.fa MF989451
   
		 #annotate MF989451 and compare to sample data
		 aligator annotate /sample_data/MF989451.fa /sample_data/MF989451.rss12_pred.bed MF989451.rss23_pred.bed IGH BK063715.fasta BK063715.bed --alleledb coding.fa --outgff annotations.gff --outfasta IgGenes.fa --blast blastn

