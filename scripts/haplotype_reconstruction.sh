#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 64G
#SBATCH --time 0-8:00
#SBATCH --array=1-20

################################################################################
# Reconstruct haplotypes for one project. Run one chromosome per task.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-26
################################################################################

##### VARIABLES #####

set -e -u -o pipefail

PROJECT=DO_Arsenic

CHR=${SLURM_ARRAY_TASK_ID}

if [ ${SLURM_ARRAY_TASK_ID} -eq 20 ]
then
  CHR=X
fi

# R container with qtl2.
R=/compsci/gedi/DO_founder_freq/containers/r_qtl2.sif

##### MAIN #####

module load singularity

singularity exec -B /compsci,/gedi ${R} Rscript haplotype_reconstruction.R ${PROJECT} ${CHR}

