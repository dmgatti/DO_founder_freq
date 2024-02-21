#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 128G
#SBATCH --time 0-12:00

################################################################################
# Use argyle to get the log2-ratio of intensities and the B-allele frequency.
# These are described in https://pubmed.ncbi.nlm.nih.gov/26684930/.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-12-20 
################################################################################

##### VARIABLES #####

R=/compsci/gedi/DO_founder_freq/containers/argyle.sif

##### MAIN #####

module load singularity

singularity exec -B /compsci,/gedi ${R} Rscript get_lrr_baf.R
