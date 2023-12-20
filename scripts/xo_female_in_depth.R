################################################################################
# After the XO females have been mapped, look at each peak over LOD = 7.
# We are having a hard time estimating model coefficients, so look at the 
# founder allele frequencies in the XO and non-XO groups.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-12-19
################################################################################

##### VARIABLES #####

# Base directory.
base_dir = '/compsci/gedi/DO_founder_freq'

# Data directory.
data_dir = file.path(base_dir, 'data')

# Array intensity directory.
intensity_dir = file.path(data_dir, 'intensities')

# XO mapping results directory.
xo_dir = file.path(base_dir, 'results', 'map_xo_females')

# LOD files from mapping.
lod_files = dir(xo_dir, pattern = '^XO_lod_chr', full.names = TRUE)

# Map files from mapping.
map_files = dir(xo_dir, pattern = '^XO_map_chr', full.names = TRUE)

# Phenotype file.
pheno_file = file.path(xo_dir, 'pheno_xo.rds')

# Allele probs files.
probs_files = dir(xo_dir, pattern = '_alleleprobs_female\\.rds$', full.names = TRUE)

##### MAIN #####

# Load in phenotypes.
pheno = readRDS(pheno_file)

result = NULL

chromosomes = c(1:19, 'X')

for(chr in chromosomes) {

print(chr)

  # Read in the LOD and map for this chromosome.
  lod = readRDS(lod_files[grep(paste0('chr', chr, '\\.rds$'), lod_files)])
  map = readRDS(map_files[grep(paste0('chr', chr, '\\.rds$'), map_files)])

  peaks = find_peaks(lod, map, threshold = 7, peakdrop = 3)
  
  if(nrow(peaks) > 0) {

    probs = readRDS(probs_files[grep(paste0('chr', chr, '_'), probs_files)])[[1]]

    for(p in 1:nrow(peaks)) {
    
      mkr = which(map[[chr]] == peaks$pos[p])
      mkr = names(map[[chr]])[mkr]
      pr    = probs[,,mkr]
    
      for(f in 1:ncol(pr)) {
    
        tt = t.test(x = pr[pheno[,1] == 0, f], y = pr[pheno[,1] == 1, f])
        result = rbind(result,
                       data.frame(chr = chr, pos = peaks$pos[p], 
                                  marker = mkr, lod = peaks$lod[p],
                                  founder = colnames(pr)[f],
                                  non_xo = tt$estimate[1], xo = tt$estimate[2],
                                  non_xo_se = tt$stderr[1], xo_se = tt$stderr[1],
                                  stat = tt$statistic, pv = tt$p.value))
    
      } # for(f)
    } # for(p)

  } # if(nrow(peaks) > 0)

} # for(chr)

write.csv(result, file.path(xo_dir, 'xo_founder_allele_freq.csv'), quote = F, row.names = F)


