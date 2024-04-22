#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=100
#SBATCH -o test

#activate conda
module load anaconda3/2022.10
conda activate compMethods

Rscript subtree_hapauto_mapping.R
