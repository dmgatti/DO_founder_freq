#!/bin/bash
#SBATCH --qos batch
#SBATCH --nodes=1 # number of nodes
#SBATCH --ntasks=1 # number of cores
#SBATCH --mem=375G # memory pool for all cores
#SBATCH --time=0-8:00 # time (D-HH:MM)
#SBATCH --array=1-20

##### VARIABLES #####

# Current chromosome.
CHR=${SLURM_ARRAY_TASK_ID}

if [ ${CHR} -eq 20 ]
then
  CHR=X

fi

# qtl2 directory with JSON files.
JSON_DIR=/compsci/gedi/DO_founder_freq/data/qtl2

# Current JSON file.
JSON=${JSON_DIR}/chr${CHR}.json

# Current Project (top-level genotype directory)
PROJECT=209_DO_Pack_Sleep 

# R container.
CONTAINER=/compsci/gedi/DO_founder_freq/containers/r_qtl2.sif

# R script for haplotype reconstruction.
RSCRIPT=/compsci/gedi/DO_founder_freq/scripts/haplotype_reconstruction.R


##### MAIN #####

module load singularity

singularity exec -B /compsci ${CONTAINER} Rscript ${RSCRIPT} ${JSON} ${PROJECT} ${CHR}

