#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 64G
#SBATCH --time 0-24:00

##### VARIABLES #####

R=/compsci/gedi/DO_founder_freq/containers/r_qtl2.sif

##### MAIN #####

cd ${SLURM_SUBMIT_DIR}

module load singularity

singularity exec -B /compsci ${R} Rscript map_xo_females_perm.R

