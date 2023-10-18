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

# Crossover results directory.
crossover_dir = file.path(results_dir, 'crossovers')

# Genoprobs directory.
probs_dir = file.path(results_dir, 'genoprobs')

# GigaMUGA marker file.
marker_file = file.path(data_dir, 'gm_uwisc_v4.csv')

# Sample metadata file.
metadata_file = file.path(data_dir, 'gigamuga_sample_metadata.csv')

##### FUNCTIONS #####

# Get crossover locations using qtl2 functions. 
# Arguments:
# project: character string containing the top level project directory name.
# markers: data.frame containing marker information for GigaMUGA. 
get_crossovers = function(project, markers) {

  # Get a listing of the genoprobs objects.
  probs_files = dir(probs_dir, pattern = project,  full.names = TRUE)
  probs_files = probs_files[grep('_genoprobs_', probs_files)]
  
  # Remove Chr X for now.
  # TBD: Focus on Chr X at some point.
  probs_files = probs_files[-grep('chrX', probs_files)]

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


# This is an alternate way of getting crossovers which uses matrix subtraction.
# It calls all crossovers and get the founder diplotypes on each side of each
# crossover.
# Arguments:
# project: character string containing the top level project directory name.
# markers: data.frame containing marker information for GigaMUGA. 
get_crossovers_with_diplotypes = function(project, markers) {

  # Get a listing of the genoprobs objects.
  probs_files = dir(probs_dir, pattern = paste0(project, '_genoprobs_'),
                    full.names = TRUE)
  
  # Remove Chr X for now.
  # TBD: Focus on Chr X at some point.
  probs_files = probs_files[-grep('chrX', probs_files)]

  # Extract chromosome from filenames.
  chr = str_replace_all(basename(probs_files), str_c('^', project, '_genoprobs_chr|\\.rds$'), '')

  # Create map object.
  markers = markers %>%
              filter(chr %in% c(1:19, 'X')) %>%
              mutate(pos = bp_grcm39 * 1e-6) %>%
              as.data.frame()
  map     = map_df_to_list(map = markers, pos_column = 'pos')  

  results = NULL

  # Process one chromosome at a time.
  for(i in seq_along(probs_files)) {

    print(chr[i]) 

    # Read in genoprobs for this chromosome.
    probs = readRDS(probs_files[i])[[1]]
    probs = round(2 * probs)
    
    # Synch up markers in probs and map.
    map[[chr[i]]] = map[[chr[i]]][dimnames(probs)[[3]]] 

    xo = apply(probs, 1, get_xo_one_sample, map = map[[chr[i]]])
    
    for(j in 1:nrow(probs)) {
      get_xo_one_sample(probs[j,,], map[[chr[i]]])
    }

  } # for(i)

} # get_crossovers_with_diplotypes()


# Arguments:
# pr: Named, numeric matrix containing the 36-state genoprobs for one sample. 
#     Diplotypes in rows and markers in columns.
# map: Named, numeric vector. qtl2-style marker map for one chromosome.
get_xo_one_sample = function(pr, map) {

  # Get the sum of the rounded probs at each marker.
  pr_sum = colSums(pr)
  pr_max = apply(pr, 2, max)
  pr_min = apply(pr, 2, min)

  # Get the genoprobs difference between adjacent markers.
  # If a row conatins all zeros, then there was no change between the
  # two adjacent markers.
  d = diff(t(pr))
  
  # Get the rows with non-zero values.
  wh = which(rowSums(d != 0) > 0)
  
  # Result should contain the mean crossover position between markers,
  # proximal and distal markers, and proximal and distal diplotypes.
  result = data.frame(prox_mkr   = rep('', length(wh)),
                      dist_mkr   = rep('', length(wh)),
                      position   = rep(0,  length(wh)),
                      prox_diplo = rep('', length(wh)),
                      prox_diplo = rep('', length(wh)))

  # Row index for results data.frame
  idx = 1
  
  for(i in wh) {
  
    prox_mkr = colnames(pr)[i]
    dist_mkr = colnames(pr)[i + 1]
  
    if(pr_max[prox_mkr] == 2 & pr_max[dist_mkr] == 2) {
      pos        = (map[prox_mkr] + map[dist_mkr]) * 0.5
      prox_diplo = rownames(pr)[which.max(pr[,prox_mkr])]
      dist_diplo = rownames(pr)[which.max(pr[,dist_mkr])]
    } else {
      pos = NA
      prox_diplo = NA
      dist_diplo = NA
    } # else
  
    # Add values to results.
    result$prox_mkr[idx]   = prox_mkr 
    result$dist_mkr[idx]   = dist_mkr
    result$position[idx]   = pos
    result$prox_diplo[idx] = prox_diplo
    result$dist_diplo[idx] = dist_diplo
  
    idx = idx + 1

  } # for(i)
  
  result = result[order(result$pos),]
  
  ####################################
  # rle approach
  haps = apply(pr, 1, rle)
  haps = haps[sapply(haps, function(z) { length(z$values) > 1 })]
  
  df = data.frame(gt  = rep(names(haps), sapply(haps, function(z) { length(z$lengths) })),
                  val = unlist(sapply(haps, function(z) { z$values })),
                  idx = unlist(sapply(haps, function(z) { cumsum(z$lengths) })))
  df = df[order(df$idx, df$val),]
  
  tbl = table(df$idx)
  
  result = NULL
  
  for(i in seq_along(tbl)) {
  
    # Current marker index.
    idx = as.numeric(names(tbl)[i])
    # Rows in table with current index.
    wh  = which(df$idx == idx)
  
    # Normal case with one, clear genotype.
    if(all(df$val[wh] %in% c(0, 2))) {
      
      prox_mkr = colnames(pr)[df$idx[wh[1]]]
      dist_mkr = colnames(pr)[df$idx[wh[2]]]
      pos      = (map[prox_mkr] + map[dist_mkr]) * 0.5
      prox_diplo = df$gt[wh[1]]
      dist_diplo = df$gt[wh[2]]
      
    } else {

      # Find the next pair of 0/2 values at the same marker.
      next_idx = which(tbl == 2 & 1:length(tbl) > i)[1]

    } # else
    
    result = rbind(result, data.frame(prox_mkr, dist_mkr, pos, prox_diplo, dist_diplo))
    
  } # for(i)

  return(result)

} # get_xo_one_sample()



##### MAIN #####

# Read in the markers.
markers = read_csv(marker_file, show_col_types = FALSE) %>%
            select(marker:snp) %>%
            filter(chr %in% c(1:19, 'X')) %>%
            mutate(pos = bp_grcm39 * 1e-6)

# Read in sample metadata.
metadata = read_csv(metadata_file, show_col_types = FALSE)

# Get the projects from the genoprobs.
projects = dir(probs_dir, pattern = 'genoprobs')
projects = unique(sub('_genoprobs_chr[0-9|X]+\\.rds$', '', projects))

for(i in seq_along(projects)) {

  print(projects[i])

  # Get the crossovers for the current project.
  cr = get_crossovers(project = projects[i], markers = markers)
  
  # Save crossover counts and positions separately.
  saveRDS(cr$num_xo, file = file.path(crossover_dir, str_c(projects[i], '_num_xo.rds')))
  saveRDS(cr$xo_pos, file = file.path(crossover_dir, str_c(projects[i], '_xo_pos.rds')))

} # for(i) 



