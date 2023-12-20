########################################################################
# Read in the DO variants and filter to retain a somewhat reasonable 
# subset of potential de-novo variants.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2023-12-21
########################################################################

library(tidyverse)

##### VARIABLES #####

# Base project directory.
base_dir = '/compsci/gedi/DO_founder_freq'

# Results directory.
results_dir = file.path(base_dir, 'results', 'do_variants')

# Variant files.
var_files = dir(results_dir, pattern = '^do_snps_chr', full.names = TRUE)

# Corresponding chromosomes.
chrs = gsub('^/compsci/gedi/DO_founder_freq/results/do_variants/do_snps_chr|\\.csv$', '', var_files)

##### MAIN #####

# Read in the variants files.
vars = lapply(var_files, read_csv, col_types = c('ccdcclccidddddddidciiiiiiciiiiiicccccccccccccccccccccccccccccccccccccccccccccccc'))
vars = do.call(bind_rows, vars)
colnames(vars)[1] = 'snp_id'

# How many variants are there?
dim(vars)

# Retain SNPs with depths less than 2X and positive values in 0/0, 0/1, & 1/1.
vars = vars %>%
         filter(DP <= 96 & `0/0` > 0 & `0/1` > 0 & `1/1` > 0)

# How many variants are there now?
dim(vars)

