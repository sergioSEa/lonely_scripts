Scripts from **Multiple origins of a frameshift insertion in a mitochondrial gene in birds and turtles**
Sergio Andreu-Sanchez,, Wanjun Chen, Josefin Stiller, Guojie Zhang.

All scripts here presented were written by S. Andreu-Sanchez.


The scripts are divided in three parts:
* Preparation of the multiple sequence alignment file and table with infomation on position 174. This can be achieved by running Prepare_dataset.sh. This script is to be run using three sets of gbf files, one for Refseq records, other for the whole genbank and a third for records belonging to the B10K project.

For all subsequent analysis we work with a subset of those records that only contain Diapsiada (data available).

* Obtain a phylogenetic tree for all Diapsiada. This was done using Make_otol_tree.R and the R package Rotol. The resulting tree is available with the data.
* Analysis of the sequence and phylogenetic patterns of the insertion. This is made by using Analysis_Diapsida.R. This script is composed of different functions for each of the analysis from the manuscript. The input data for them is made available. The scripts for generating that data are included in the scripts directory.

