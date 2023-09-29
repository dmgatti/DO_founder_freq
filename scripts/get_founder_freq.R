################################################################################
# Given a 36-state founder allele probs object, calculate the founder allele
# frequency along the genome.
# Split the samples up by DO generation.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-27
################################################################################

library(qtl2)
library(tidyverse)

##### VARIABLES #####

# Top-level directory for the project.
base_dir = '/compsci/gedi/DO_founder_freq'

# Data directory.
data_dir = file.path(base_dir, 'data')

# Results directory.
results_dir = file.path(base_dir, 'results')

# Metadata directory.
metadata_dir = file.path(data_dir, 'metadata')

# Genoprobs directory.
probs_dir = file.path(results_dir, 'genoprobs')




##### MAIN #####


# Get the founder allele
get_founder_freq = function(meta) {

  # Get a listing of the allele probs objects.
  probs_files = dir(probs_dir, pattern = project,  full.names = TRUE)
  probs_files = probs_files[grep('_alleleprobs_', probs_files)]

  # Extract chromosome from filenames.
  chr = str_replace_all(basename(probs_files), str_c('^', project, '_alleleprobs_chr|\\.rds$'), '')

  for(i in seq_along(probs_files)) {
  
    print(chr[[i]])

    probs = readRDS(probs_files[i])

    

  } # for(i)

} # get_founder_freq()

