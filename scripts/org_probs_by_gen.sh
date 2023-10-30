#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 32G
#SBATCH --time 0-8:00

################################################################################
# Reorg.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-26
################################################################################

##### VARIABLES #####

set -e -u -o pipefail

# R container with qtl2.
R=/compsci/gedi/DO_founder_freq/containers/r_qtl2.sif

##### MAIN #####

module load singularity

singularity exec -B /compsci ${R} Rscript org_probs_by_gen.R

