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


get_crossovers = function(project) {

  # Get a listing of the cross objects.
  cross_files = dir(file.path(cross_dir, str_c(project, '*')), full.names = TRUE)

  # Get a listing of the probs objects.
  probs_files = dir(file.path(probs_dir, str_c(project, '*')), full.names = TRUE)

  stopifnot(length(cross_files) == length(probs_files)))

  for(i in seq_along(cross_files)) {

    # Read cross object.
    cross = readRDS(file.path(cross_dir, ''))

    # Read genoprobs object.
    probs = readRDS(file.path(probs_dir, ''))

  } # for(i)

} # get_crossovers()

