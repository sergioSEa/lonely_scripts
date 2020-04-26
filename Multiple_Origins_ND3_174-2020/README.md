Scripts from **Multiple origins of a frameshift insertion in a mitochondrial gene in birds and turtles**
Sergio Andreu-Sanchez,, Wanjun Chen, Josefin Stiller, Guojie Zhang.

All scripts here presented were written by S. Andreu-Sanchez.


The scripts are divided in three parts:
1. Preparation of the multiple sequence alignment file and table with infomation on position 174. This can be achieved by running **Prepare_dataset.sh**. This script is to be run using three sets of gbf files, one for Refseq records, other for the whole genbank and a third for records belonging to the B10K project.
* Requirements: Bash, MAFFT, pxclsq, pytho3, biopython

For all subsequent analysis we work with a subset of those records that only contain Diapsida (data available).

2. Obtain a phylogenetic tree for all Diapsiada. This was done using **Make_otol_tree.R** and the R package Rotol. The resulting tree is available with the data. The open tree of life version 11.4, labelled_supertree_ottnames.tre was used as input.
* Requirements: R, tidyverse, Rotl, ape, phytools
3. Analysis of the sequence and phylogenetic patterns of the insertion. This is made by using **Analysis_Diapsida.R**. This script is composed of different functions for each of the analysis from the manuscript. The input data for them is made available. The scripts for generating that data are included in the scripts directory.
* Requirements: R, tidyverse, ape, phytools, ggtree, ggimage, phangorn, ggstance, Biostrings, castor

