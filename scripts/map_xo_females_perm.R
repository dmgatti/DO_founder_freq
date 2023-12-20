################################################################################
# Gather the females and map XO as a binary trait.
# Note that you need to use the -B /compsci argument in singularity.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-12-07
################################################################################

library(qtl2convert)
library(qtl2)

##### VARIABLES #####

# Base project directory.
base_dir = '/compsci/gedi/DO_founder_freq'

# Data directory.
data_dir = file.path(base_dir, 'data')

# Results directory.
results_dir = file.path(base_dir, 'results')

# Genoprobs directory.
probs_dir = file.path(results_dir, 'genoprobs')

# Output directory.
out_dir = file.path(results_dir, 'map_xo_females')


##### MAIN #####

# Read in phenotypes.
pheno = readRDS(file = file.path(out_dir, 'pheno_xo.rds'))

# Run permutations on each chromosome.
for(chr in c(1:19, 'X')) {

  print(chr)
  t1 = proc.time()[3]
  
  probs = readRDS(file = file.path(out_dir, paste0('chr', chr, '_alleleprobs_female.rds')))

  perms = scan1perm(probs, pheno, model = 'binary', n_perm = 1000, cores = 20)
  
  saveRDS(perms, file = file.path(out_dir, paste0('perms_chr', chr, '.rds')))

  print(proc.time()[3] - t1)

} # for(chr)

