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
probs_dir = file.path(results_dir, 'probs_by_gen')

# Output directory for generation linear model results.
out_dir = file.path(results_dir, 'gen_lm')

##### MAIN #####

# Get the unique generations from the allele probs filenames.
probs_files = dir(probs_dir)
probs_files = str_replace_all(probs_files, 'gen|_alleleprobs_chr[0-9,X]+\\.rds', '')

unique_gen  = sort(as.numeric(unique(probs_files)))

# Just keep G22 and higher for now.
unique_gen = unique_gen[unique_gen >= 22]

# Read in each chomosome and fit a model at each marker of
# probs ~ gen. Save the p-value and R^2.
for(chr in c(1:19, 'X')) {
 
  print(str_c('CHR ', chr))

  # Get the probs filenames for thie chromosome.
  probs_files = dir(probs_dir, pattern = str_c('chr', chr, '\\.rds'), full.names = TRUE)
  probs_files = probs_files[grep(paste0(unique_gen, collapse = '|'), probs_files)]

  # Set the generation as the names of probs_files. 
  names(probs_files) = str_replace_all(probs_files, 
                                       str_c('_alleleprobs_chr', chr,'.rds'), '')
  names(probs_files) = str_sub(names(probs_files), start = -2, end = -1)

  # Read in the allele probs for this chromosome.
  probs = lapply(probs_files, readRDS)
  
  # Get the number of samples in each generation.
  n_samples = sapply(probs, dim)[1,]

  # Create a generation vector for the linear model.
  gen_vec = rep(as.numeric(names(probs)), n_samples)
 # gen     = unique(gen_vec)

  # Get the number of markers.
  markers   = dimnames(probs[[1]])[[3]]
  n_markers = dim(probs[[1]])[3]

  # Create results data structure.
  lm_result = array(0, dim = c(n_markers, 3, 8), 
                    dimnames = list(markers, c('slope', 'rsq', 'p.value'), LETTERS[1:8]))
  
  # TBD: Would tidyverse be easier? more data-wrangling to get data
  # into a format that tidyverse can us.
  
  for(i in 1:nrow(lm_result)) {
  
    if(i %% 100 == 0) print(str_c('   ', i))
  
    # Get the probs into one data structure.
    mkr = markers[i]
    pr  = sapply(probs, function(z) { z[,,mkr] })
    pr  = do.call(rbind, pr)
    #pr  = data.frame(pr)
  
    # Get the mean founder frequency in each generation.
    # mean_freq = t(sapply(split(pr, gen_vec), colMeans))
    
    # Fit a model of probs ~ gen.
 #   mod  = lm(mean_freq ~ gen)
    mod  = lm(pr ~ gen_vec)
    smry = summary(mod)

    # Extract results.
    # slope
    lm_result[i,'slope',]   = sapply(smry, function(z) { z$coefficients[2,1] })
    # adjusted R squared
    lm_result[i,'rsq',]     = sapply(smry, function(z) { z$adj.r.squared })
    # p.value
    lm_result[i,'p.value',] = sapply(smry, function(z) { z$coefficients[2,4] })
  
  } # for(i)

  saveRDS(lm_result, file = file.path(out_dir, str_c('lm_result_chr', chr, '.rds')))  

} # for(chr)


