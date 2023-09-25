#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 64G
#SBATCH --time 0-8:00

##### VARIABLES #####

# Neogen directories file.
NEOGEN_DIRS=/compsci/gedi/DO_founder_freq/data/neogen_dirs.txt

# Output directory.
OUTPUT_DIR=/compsci/gedi/DO_founder_freq/data

# R container.
R=/compsci/gedi/DO_founder_freq/containers/r_qtl2.sif


##### MAIN #####

module load singularity

singularity exec -B /gedi,/compsci ${R} Rscript --vanilla prescreen_samples.R ${NEOGEN_DIRS} ${OUTPUT_DIR}

