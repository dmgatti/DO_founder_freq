################################################################################
# Given a cross object and a genoprobs object, estimate the number of crossovers
# in each mouse in the data set.
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

# Cross object directory.
cross_dir = file.path(results_dir, 'cross_obj')

# Genoprobs directory.
probs_dir = file.path(results_dir, 'genoprobs')


##### MAIN #####


get_crossovers = function(project, markers) {

  # Get a listing of the cross objects.
  cross_files = dir(cross_dir, pattern = project, full.names = TRUE)

  # Get a listing of the probs objects.
  probs_files = dir(probs_dir, pattern = project,  full.names = TRUE)
  probs_files = probs_files[grep('_genoprobs_', probs_files)]

  stopifnot(length(cross_files) == length(probs_files))

  # Extract chromosome from filenames.
  chr = str_replace_all(basename(cross_files), str_c('^', project, '_cross_chr|\\.rds$'), '')

  for(i in seq_along(cross_files)) {

    # Read cross object.
    cross = readRDS(cross_files[i])

    # Read genoprobs object.
    probs = readRDS(probs_files[i])

    maxgt  = maxmarg(probs)
    num_xo = count_xo(maxgt)

  } # for(i)

} # get_crossovers()

