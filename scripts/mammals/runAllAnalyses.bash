#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=100
#SBATCH -o runAll

#activate conda
module load anaconda3/2022.10
conda activate compMethods

#run mammal wide scripts
Rscript multiSAF_mcmc.R
Rscript multiSAF_mapping.R
Rscript hapauto_mcmc_base_params.R
Rscript hapauto_mcmc_adjusted_prior.R
Rscript hapauto_mcmc_final.R
Rscript hapauto_mapping.R
Rscript propSAF_analysis.R
Rscript propSAF_figures.R

#run subtree scripts
Rscript subtree_hapauto_mcmc.R
Rscript subtree_hapauto_mapping.R
Rscript subtree_multiSAF_mcmc.R
Rscript subtree_multiSAF_mapping.R
Rscript subtree_propSAF_analysis.R
Rscript subtree_propSAF_figures.R
Rscript subtree_propSAF_scatter.R



