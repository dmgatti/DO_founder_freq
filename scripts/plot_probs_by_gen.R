################################################################################
# Plot the founder allele frequency for each chromosome by generation.
# 
# Daniel Gatti
# 2024-02-13
# dan.gatti2jax.org
################################################################################

library(tidyverse)
library(qtl2)


##### VARIABLES #####

base_dir = '/compsci/gedi/DO_founder_freq'

figure_dir = file.path(base_dir, 'figures')

results_dir = file.path(base_dir, 'results')

probs_dir = file.path(results_dir, 'probs_by_gen')

marker_dir = '/compsci/gedi/reference'

marker_file = file.path(marker_dir, 'gm_uwisc_v4.csv')

##### MAIN #####

# Read in the GigaMUGA markers.
markers = read.csv(marker_file) %>%
            select(marker, chr, pos = bp_grcm39) %>%
            mutate(pos = pos * 1e-6)

# Get the probs file names.
probs_files = dir(probs_dir, full.names = TRUE)

# Get the chromosome and generation for each file.
tmp = basename(probs_files)
tmp = sub('\\.rds$', '', tmp)
tmp = strsplit(tmp, '_')

probs_files = data.frame(fn  = probs_files,
                         gen = sapply(tmp, '[', 1),
                         chr = sapply(tmp, '[', 3))
probs_files$gen = sub('^gen', '', probs_files$gen)
probs_files$chr = sub('^chr', '', probs_files$chr)

# Remove generation 14 because it has too few samples.
probs_files = subset(probs_files, gen > 14)

# For each chromosome, plot the founder allele frequencies for each founder.
unique_gen = unique(probs_files$gen)
unique_chr = unique(probs_files$chr)

num_gen = length(unique_gen)
num_chr = length(unique_chr)
gen_colors = colorRampPalette(c('gray30', 'black'))(num_gen)
names(gen_colors) = unique_gen

for(chr in unique_chr) {

  print(chr)

  # Read in the files for the current chromosome.
  curr_files = probs_files[probs_files$chr == chr,]
  pr = lapply(curr_files$fn, readRDS)
  
  # Get the founder allele frequency at each marker.
  pr = lapply(pr, function(z) { apply(z, 2:3, mean) })
  names(pr) = curr_files$gen

  # Add the generation to each founder freq.
  freq = NULL
  for(i in seq_along(pr)) {
  
    pr[[i]] = data.frame(marker = colnames(pr[[i]]), 
                         gen    = names(pr)[i], 
                         t(pr[[i]]))

  } # for(i)

  freq = do.call(rbind, pr)

  # Add the marker positions in Mb.
  freq = freq %>%
           left_join(markers, by = 'marker')

  # Make the plot.
  png(file.path(figure_dir, paste0('founder_freq_chr', chr, '.png')),
      width = 2000, height = 2000)
  freq %>%
    pivot_longer(cols = A:H, names_to = 'founder', values_to = 'freq') %>%
    ggplot(aes(pos, freq, color = gen)) +
      geom_line() +
      facet_wrap(~founder, ncol = 1)
  dev.off()

  saveRDS(freq, file = file.path(results_dir, 'founder_freq', paste0('founder_freq_chr', chr, '.rds')))

  rm(pr)
  gc()

} # for(chr)
