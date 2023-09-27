################################################################################
# Reconsruct the haplotypes of a data set.
# This will point to a JSON file which will describe the input files for qtl2.
# 
# Daniel Gatti
# Dan.gatti#jax.org
# 2023-09-26
################################################################################

library(qtl2convert)
library(qtl2)

##### VARABLES #####

# Command line arguments.
# args[1]: project name.
# args[2]: chromosome
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 2) {

  print('usage: haploype_reconstruction.R <PROJECT> <CHR>')
  quit(1)

} # if(length(args) != 2)

project = args[1]
chr     = args[2]

# Base directory.
base_dir = '/compsci/gedi/DO_founder_freq'

# Directory containing qtl2 files.
qtl2_dir = file.path(base_dir, 'data', 'qtl2')

# Results directory.
results_dir = file.path(base_dir, 'results')

# Cross object directory.
cross_dir = file.path(results_dir, 'cross_obj')

# Genoprobs direcetory.
probs_dir = file.path(results_dir, 'genoprobs')

##### MAIN #####

# Haplotype reconstruction.

print('Reading cross')
cross = read_cross2(file.path(qtl2_dir, paste0('chr', chr, '.json')))

print(cross)

print('Writing cross object.')
saveRDS(cross, file = file.path(cross_dir, paste0(project, '_cross_chr', chr, '.rds')))

print('Running HR')
probs = calc_genoprob(cross, map = cross$gmap, cores = 1, quiet = FALSE)

print('Writing 36 state genoprobs')
saveRDS(probs, file = file.path(probs_dir, paste0(project, '_genoprobs_chr', chr, '.rds')))

print('Converting 36 state genoprobs to 8 state allele probs')
aprobs = genoprob_to_alleleprob(probs, cores = 1, quiet = FALSE)

print('Writing allele probs')
saveRDS(aprobs, file = file.path(probs_dir, paste0(project, '_alleleprobs_chr', chr, '.rds')))

print('COMPLETED')

