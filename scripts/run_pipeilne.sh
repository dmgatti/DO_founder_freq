#!/bin/bash

##### VARIABLES #####

# File containing the Neogen directories that we will use.
INPUT_FILE=/compsci/gedi/DO_founder_freq/data/neogen_dirs.txt

# Directory from which Neogen genotypes will be read.
SOURCE_DIR=/gedi/resource/MUGA/raw_neogen_data/GigaMUGA

# Directory where genotypes will be written.
DEST_DIR=/compsci/gedi/DO_founder_freq/data/genotypes

R=/compsci/gedi/DO_founder_freq/containers/r_qtl2.sif

##### MAIN #####

module load singularity

# Create founder genotypes.

DIRS=`cat ${INPUT_FILE}`


  CURR_DIR=${DIRS[$INDEX]}

  echo ${CURR_DIR}

  OUT_DIR=${DEST_DIR}/${CURR_DIR}

  mkdir -p ${OUT_DIR}

  LOCAL_DIRS=(`ls -1d ${SOURCE_DIR}/${CURR_DIR}/*/ | xargs echo | sed 's/ /,/g'`)

  # Read in the Neogen files and output genotypes and intensities.
  singularity exec -B /gedi ${R} Rscript read_neogen.R ${LOCAL_DIRS} ${OUT_DIR}
  
  # Infer sex and calculate call rates. 
  
  # Prepare genotypes for qtl2.

  # Estimate Chr Y & M haplotypes.
  
  # Run qtl2 and calculate genoprobs.
  
  # Count & locate crossovers. (Need probs & cross object)
  
  # QC based on crossover count.
  
  # Clean up. (genoprobs, genotypes, intensities)

done
