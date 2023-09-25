################################################################################
# Haplotype Reconstruction on GRCm38 forDO Superovulation study.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-12-06
################################################################################

library(qtl2convert)
library(qtl2)

##### VARIABLES #####

args = commandArgs(trailingOnly = TRUE)

json_file = args[1]

project   = args[2]

chr       = args[3]

base_dir    = '/compsci/gedi/DO_founder_freq'

results_dir = file.path(base_dir, 'results')

cross_dir   = file.path(results_dir, 'cross_obj')

probs_dir   = file.path(results_dir, 'probs')


##### MAIN #####

# Haplotype reconstruction.
print('Reading cross')
cross = read_cross2(json_file)

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


