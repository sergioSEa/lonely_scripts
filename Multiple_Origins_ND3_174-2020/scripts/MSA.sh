#!/bin/bash
#SBATCH --partition normal 
#SBATCH --mem-per-cpu 5G #memory per CPU
#SBATCH -c 6 #one core
#SBATCH -t 24:00:00

## Argument 1: Merged ND3, eg: ND3_merged.fa
mafft --localpair --maxiterate 1000 --thread 6 $1  > ND3_msa.fa

