################################################################################
# Gather the haplotype-level inbreeding degree for each mouse and compile them.
# Then reorganize them based on DO generation.
# Use at least 64G memory.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2023-12-18
################################################################################

library(qtl2)

##### VARIABLES #####

# Base project directory.
base_dir = '/compsci/gedi/DO_founder_freq'

# Results directory.
results_dir = file.path(base_dir, 'results')

# Genoprobs directory.
probs_dir = file.path(results_dir, 'genoprobs')

# Inbreeding output directory.
inbreeding_dir = file.path(results_dir, 'inbreeding')

# Metadata directory.
meta_dir = file.path(base_dir, 'data', 'metadata')


##### MAIN #####

# Read in the sample metadata files.
meta_files = dir(path = meta_dir, pattern = 'csv$', full.names = TRUE)

meta = lapply(meta_files, read.csv)
meta = do.call(rbind, meta)

# Get the genoprobs files.
probs_files = dir(probs_dir, pattern = 'rds$', full.names = TRUE)
 
# Projects.
projects = unique(gsub(paste0('^', probs_dir, '/|_(allele|geno)probs_chr([0-9]+|X)\\.rds$'), '', probs_files))

# Get the 2-letter genotype codes.
# "AA" "AB" "BB" "AC" "BC" "CC" "AD" "BD" "CD" "DD" "AE" "BE" "CE" "DE" "EE"
# "AF" "BF" "CF" "DF" "EF" "FF" "AG" "BG" "CG" "DG" "EG" "FG" "GG" "AH" "BH"
# "CH" "DH" "EH" "FH" "GH" "HH"
# The homozygoutes are c(1, 3, 6, 10, 15, 21, 28, 36)
hom_index = c(1, 3, 6, 10, 15, 21, 28, 36)
result = NULL

for(proj in projects) {

  print(proj)

  # Get the files for this project.
  chr_files = dir(probs_dir, pattern = paste0(proj, '_genoprobs_'), full.names = TRUE)
  chr_files = chr_files[!grepl('chrX', chr_files)]
  
  n_markers = 0
  pr  = readRDS(chr_files[[11]])[[1]]
  hom = data.frame(id  = rownames(pr), 
                   hom = rep(0, nrow(pr)))
  
  # Read in the genoprobs files.
  for(i in seq_along(chr_files)) {
  
    print(paste('  ', i))

    pr = readRDS(chr_files[[i]])

    gt        = maxmarg(probs = pr, minprob = 0.3, return_char = FALSE)[[1]]
    n_markers = n_markers + ncol(gt)
    hom$hom   = hom$hom + rowSums(matrix(gt %in% hom_index, nrow = nrow(gt)))

  } # for(i)
  
  hom$hom = hom$hom / n_markers

  result = rbind(result, hom)
  
  rm(pr, gt)
  gc()

} # for(proj)

# Merge in the generation information.
meta           = meta[,c('Unique.Sample.ID', 'DO.Generation')]
colnames(meta) = c('id', 'gen')
result         = merge(result, meta, by = 'id')

# Save the results.
saveRDS(result, file = file.path(inbreeding_dir, 'do_haplotype_homozygosity.rds'))




