################################################################################
# Given a cross object and a genoprobs object, estimate the number of crossovers
# in each mouse in the data set.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-27
################################################################################

library(qtl2convert)
library(qtl2)
library(tidyverse)

##### VARIABLES #####

# Top-level directory for the project.
base_dir = '/compsci/gedi/DO_founder_freq'

# Data directory.
data_dir = file.path(base_dir, 'data')

# Results directory.
results_dir = file.path(base_dir, 'results')

# Genoprobs directory.
probs_dir = file.path(results_dir, 'genoprobs')


##### MAIN #####

# Arguments:
# project: character string containing the top level project directory name.
# markers: data.frame containing marker information for GigaMUGA. 
get_crossovers = function(project, markers) {

  # Get a listing of the allele probs objects.
  probs_files = dir(probs_dir, pattern = project,  full.names = TRUE)
  probs_files = probs_files[grep('_genoprobs_', probs_files)]

  # Extract chromosome from filenames.
  chr = str_replace_all(basename(probs_files), str_c('^', project, '_genoprobs_chr|\\.rds$'), '')

  # Create map object.
  markers = markers %>%
              filter(chr %in% c(1:19, 'X')) %>%
              mutate(pos = bp_grcm39 * 1e-6) %>%
              as.data.frame()
  map     = map_df_to_list(map = markers, pos_column = 'pos')

  num_xo = NULL
  xo_pos = NULL

  # Process one chromosome at a time.
  for(i in seq_along(probs_files)) {

    print(chr[i])

    # Read in genoprobs for this chromosome.
    probs = readRDS(probs_files[i])

    # Get most 
    maxgt  = maxmarg(probs, minprob = 0.95)
    xo_pos[[chr[i]]] = locate_xo(maxgt, map = map)[[1]]
    
    df = data.frame(id  = names(xo_pos[[chr[i]]]), 
                    nxo = sapply(xo_pos[[chr[i]]], length))
    colnames(df)[2] = names(probs)
    
    if(is.null(num_xo)) {
    
      num_xo = df
      
    } else {

      num_xo = full_join(num_xo, df, by = 'id')

    } # else

  } # for(i)
  
  return(list(num_xo = num_xo, xo_pos = xo_pos))

} # get_crossovers()

