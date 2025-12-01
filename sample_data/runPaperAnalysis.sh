#!/bin/bash
#
# These commands can be used to recapitulate the examples
#
# 
# Because IMGT has interfered with programmatic downloads from LIGM-DB,
#     this script is no longer sufficient on its own to fully reproduce 
#     the analysis from the paper. Instead, please first use a browser
#     to save the following URLs as local HTML files:
# https://imgt.org/ligmdb/view.action?id=BK063715 --> IGH/BK063715.html
# https://imgt.org/ligmdb/view.action?id=BK063716 --> IGK/BK063716.html
# https://imgt.org/ligmdb/view.action?id=BK063717 --> IGL/BK063717.html
# https://imgt.org/ligmdb/view.action?id=BK065026 --> TRA/BK065026.html
# https://imgt.org/ligmdb/view.action?id=BK068537 --> TRA/BK068537.html
# https://imgt.org/ligmdb/view.action?id=BK065025 --> TRB/BK065025.html
# https://imgt.org/ligmdb/view.action?id=BK068295 --> TRB/BK068295.html
#
#
# `efetch` is part of the Entrez Direct suite, please see
#     https://www.ncbi.nlm.nih.gov/books/NBK179288/
#     for more information and installation instructions
#
# `extract_refs` and `fix_macaque_gaps` are part of receptor_utils, see:
#     https://williamdlees.github.io/receptor_utils/_build/html/index.html
#
# `digger` can be found at https://williamdlees.github.io/digger/_build/html/index.html
#
# The comparison script `gatorADE` requires the R packages bedtoolsr and ggalluvial,
#     see: https://github.com/PhanstielLab/bedtoolsr and 
#          https://corybrunson.github.io/ggalluvial/


# IGH
mkdir -p IGH
cd IGH/
efetch -db nucleotide -id MF989451 -format fasta | perl -pe 's/\.1//' > MF989451.fa
efetch -db nucleotide -id MF989451 -format gb | perl ../gb2bed.pl > MF989451_geneRef.bed
python ../../aligator extract BK063715 BK063715.bed --local BK063715.html --fasta BK063715.fasta
python ../../aligator predict MF989451.fa MF989451
python ../../aligator annotate MF989451.fa MF989451.RSS12.bed MF989451.RSS23.bed IGH BK063715.fasta BK063715.bed --outgff predictedGenes.gff --outfasta predictedGenes.fasta

#also generate digger annotations for comparison
mkdir -p digger
extract_refs -L IGH "Macaca mulatta"
mv Macaca_mulatta_*fasta digger/
fix_macaque_gaps digger/Macaca_mulatta_IGHV_gapped.fasta digger/Macaca_mulatta_IGHV_gapped_fixed.fasta IGH
cat digger/Macaca_mulatta_IGHV.fasta digger/Macaca_mulatta_IGHD.fasta digger/Macaca_mulatta_IGHJ.fasta > digger/Macaca_mulatta_IGHVDJ.fasta
digger -species rhesus_macaque MF989451.fa -v_ref digger/Macaca_mulatta_IGHV.fasta -d_ref digger/Macaca_mulatta_IGHD.fasta -j_ref digger/Macaca_mulatta_IGHJ.fasta -v_ref_gapped digger/Macaca_mulatta_IGHV_gapped_fixed.fasta -ref imgt,digger/Macaca_mulatta_IGHVDJ.fasta digger/MF989451_digger.csv

Rscript ../gatorADE MF989451_geneRef.bed predictedGenes.gff differences.tsv sankey_plot.png digger/MF989451_digger.csv
cd ..

#IGK
mkdir -p IGK
cd IGK
efetch -db nucleotide -id MF989473 -format fasta | perl -pe 's/\.1//' > MF989473.fa
efetch -db nucleotide -id MF989473 -format gb | perl ../gb2bed.pl > MF989473_geneRef.bed
python ../../aligator extract BK063716 BK063716.bed --local BK063716.html --fasta BK063716.fasta
python ../../aligator predict MF989473.fa MF989473
python ../../aligator annotate MF989473.fa MF989473.RSS12.bed MF989473.RSS23.bed IGK BK063716.fasta BK063716.bed --outgff predictedGenes.gff --outfasta predictedGenes.fasta
Rscript ../gatorADE MF989473_geneRef.bed predictedGenes.gff differences.tsv sankey_plot.png
cd ..

#IGL
mkdir -p IGL
cd IGL
efetch -db nucleotide -id MF989453 -format fasta | perl -pe 's/\.1//' > MF989453.fa
efetch -db nucleotide -id MF989453 -format gb | perl ../gb2bed.pl > MF989453_geneRef.bed
python ../../aligator extract BK063717 BK063717.bed --local BK063717.html --fasta BK063717.fasta
python ../../aligator predict MF989453.fa MF989453
python ../../aligator annotate MF989453.fa MF989453.RSS12.bed MF989453.RSS23.bed IGL BK063717.fasta BK063717.bed --outgff predictedGenes.gff --outfasta predictedGenes.fasta
Rscript ../gatorADE MF989453_geneRef.bed predictedGenes.gff differences.tsv sankey_plot.png
cd ..

#TRA
mkdir -p TRA
cd TRA
python ../../aligator extract BK068537 ferret_TRA_geneRef_withPseudo.bed --local BK068537.html --fasta ferret_TRA.fasta --pseudo
python ../../aligator extract BK065026 dog_TRA.bed --local BK065026.html --fasta dog_TRA.fasta
python ../../aligator predict ferret_TRA.fasta ferret_TRA --species MM
python ../../aligator annotate ferret_TRA.fasta ferret_TRA.RSS12.bed ferret_TRA.RSS23.bed TRA dog_TRA.fasta dog_TRA.bed --outgff predictedGenes.gff --outfasta predictedGenes.fasta
Rscript ../gatorADE ferret_TRA_geneRef_withPseudo.bed predictedGenes.gff dog_differences.tsv dog_sankey_plot.png

#redo predictions with self as a positive control
#don't include pseudogenes in reference for searching
mkdir -p self
cd self
python ../../../aligator extract BK068537 ferret_TRA_geneRef_noPseudo.bed --local ../BK068537.html
python ../../../aligator annotate ../ferret_TRA.fasta ../ferret_TRA.RSS12.bed ../ferret_TRA.RSS23.bed TRA ../ferret_TRA.fasta ferret_TRA_geneRef_noPseudo.bed --outgff self_predictedGenes.gff --outfasta self_predictedGenes.fasta
Rscript ../../gatorADE ../ferret_TRA_geneRef_withPseudo.bed self_predictedGenes.gff self_differences.tsv self_sankey_plot.png
cd ../..

#TRB
mkdir -p TRB
cd TRB
python ../../aligator extract BK068295 ferret_TRB_geneRef_withPseudo.bed --local BK068295.html --fasta ferret_TRB.fasta --pseudo
python ../../aligator extract BK065025 dog_TRB.bed --local BK065025.html --fasta dog_TRB.fasta
python ../../aligator predict ferret_TRB.fasta ferret_TRB --species MM
python ../../aligator annotate ferret_TRB.fasta ferret_TRB.RSS12.bed ferret_TRB.RSS23.bed TRB dog_TRB.fasta dog_TRB.bed --outgff predictedGenes.gff --outfasta predictedGenes.fasta
Rscript ../gatorADE ferret_TRB_geneRef_withPseudo.bed predictedGenes.gff dog_differences.tsv dog_sankey_plot.png

#redo predictions with self as a positive control
#don't include pseudogenes in reference for searching
mkdir -p self
cd self
python ../../../aligator extract BK068295 ferret_TRB_geneRef_noPseudo.bed --local ../BK068295.html
python ../../../aligator annotate ../ferret_TRB.fasta ../ferret_TRB.RSS12.bed ../ferret_TRB.RSS23.bed TRB ../ferret_TRB.fasta ferret_TRB_geneRef_noPseudo.bed --outgff self_predictedGenes.gff --outfasta self_predictedGenes.fasta
Rscript ../../gatorADE ../ferret_TRB_geneRef_withPseudo.bed self_predictedGenes.gff self_differences.tsv self_sankey_plot.png
cd ../..
