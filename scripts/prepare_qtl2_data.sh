#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 32G
#SBATCH --time 0-8:00


################################################################################
# Gather data and output qtl2 files.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-26
################################################################################

##### VARIABLES #####

set -e -u -o pipefail

# Base project directory.
BASE_DIR=/compsci/gedi/DO_founder_freq

# File containing GigaMUGA project directories.
NEOGEN_FILE=${BASE_DIR}/data/neogen_dirs.txt

# R container.
R=${BASE_DIR}/containers/r_qtl2.sif

##### MAIN #####

module load singularity

PROJECTS=`cat ${NEOGEN_FILE}`
PROJECTS=($PROJECTS)

for P in ${PROJECTS[@]}
do
  echo ${P}
  singularity exec -B /compsci,/gedi ${R} Rscript prepare_qtl2_data.R ${P}
done

