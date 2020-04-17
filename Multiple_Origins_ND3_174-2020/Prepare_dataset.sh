#!/bin/bash

#Pipeline script for ND3-174 identification from GenBank records
#Multiple origins of a frameshift insertion in a mitochondrial gene in birds and turtles ; Author: S. Andreu-Sanchez
#Steps: 1. Get Fastas from GenBank 2. Select and Merge fasta records 3. Make MSA  4. Remove alignment positions not seen in 95% of sequences 5. Change gaps in position 175 to position 174 6. Assess the nucleotide found at position 174


####Prepare NCBI sequences#
echo 'ND3 extraction from genebank'
python scripts/Prepare_fasta.py > ND3_genes.fa

####First merge
echo 'Merging refseq complete mitochondria and Genebank records'
python scripts/Merge_0.py

##find taxonomy
echo 'Find taxonomy'
python scripts/find_taxa.py > Taxon_info.txt

####Merging
echo 'Merging with B10K data'
python scripts/Merge.py



###Align
echo 'Submitting MSA job'
bash scripts/MSA.sh ND3_merged.fa


###Cleaning
echo 'Cleaning up MSA'
bash scripts/Clean.sh
python scripts/Adjust_insertion.py Aligned_clean.fa > Corrected_Alignment.fa

###Get position
echo 'Creating table with 174 position'
python scripts/Find_position.py Corrected_Alignment.fa > Position_nucleotide.tsv

###Include taxonomy in position
echo 'Updating table with taxonomical information'
python Incorporate_taxonomy.py




