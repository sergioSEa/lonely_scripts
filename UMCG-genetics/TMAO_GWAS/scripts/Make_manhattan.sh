#!/bin/bash
#SBATCH --job-name=manhattan
#SBATCH --output=manhattan.out
#SBATCH --error=manhattan.err
#SBATCH --time=50:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1

module load RPlus/3.5.1-foss-2015b-v19.04.1
Rscript scripts/Data_exploration.R
