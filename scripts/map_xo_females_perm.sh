#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 256G
#SBATCH --time 0-24:00
#SBATCH --array 1-2

##### VARIABLES #####

R=/compsci/gedi/DO_founder_freq/containers/r_qtl2.sif

##### MAIN #####

cd ${SLURM_SUBMIT_DIR}

module load singularity

singularity exec -B /compsci ${R} Rscript map_xo_females_perm.R ${SLURM_ARRAY_TASK_ID}


