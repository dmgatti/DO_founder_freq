################################################################################
# Read in the allele probs for each project and write them out in one file per
# generation.
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
probs_dir = file.path(results_dir, 'genoprobs')

# Output direcotory.
out_dir = file.path(results_dir, 'probs_by_gen')

##### MAIN #####


# Read in the sample metadata files.
metadata_files = dir(metadata_dir, pattern = 'csv', full.names = TRUE)

metadata = lapply(metadata_files, read_csv, show_col_types = FALSE)
metadata = do.call(rbind, metadata)
metadata = metadata %>%
             select(id        = Unique.Sample.ID,
                    project   = Project, 
                    directory = Directory,
                    sample_id = Sample.ID,
                    sex       = Sex,
                    gen       = DO.Generation) %>%
             filter(gen > 1)

# Split sample metadata by DO generation.
metadata = split(metadata, metadata$gen)

# For each generation, for each chromosome, read in the allele probs from 
# each project and write them out.
for(gen in names(metadata)) {

  # Get the sample metadata for the current generation.
  curr_meta = metadata[[gen]]
  
  print(str_c('Gen: ', gen, '  n = ', nrow(curr_meta)))

  # Get the unique projects.
  projects = unique(curr_meta$project)

  # Get the allele probs filenames for the first project.
  probs_files = dir(probs_dir, pattern = str_c(projects[1], '_alleleprobs_chr'), 
                    full.names = TRUE)

  # Extract chromosome from filenames.
  chr = str_replace_all(basename(probs_files), str_c('^', projects[1], '_alleleprobs_chr|\\.rds$'), '')

  for(c in seq_along(chr)) {
  
    # Read in one chromosome to set up data structures.
    pr = readRDS(probs_files[c])
  
    # Use the sample metadata and the size of this chromosome to craete
    # a data structure for the new data.
    new_probs = array(0, dim = c(nrow(curr_meta), 8, dim(pr[[1]])[[3]]), 
                      dimnames = list(curr_meta$id, LETTERS[1:8], dimnames(pr[[1]])[[3]]))
  
    start_idx = 1
    end_idx   = 0
  
    for(p in projects) {
  
      # Get the probs files for this project.
      probs_files = dir(probs_dir, pattern = str_c(p, '_alleleprobs_chr'), 
                        full.names = TRUE)
      
      # Read in the file for the current chromosome.
      # We're assuming that all of the projects have the same chromosome order.
      pr = readRDS(probs_files[c])[[1]]

      # Subset the probs to contain only the samples in the metadata. 
      pr = pr[curr_meta$id[curr_meta$project == p],,]
    
      end_idx = start_idx + nrow(pr) - 1
      new_probs[start_idx:end_idx,,] = pr
      
      start_idx = end_idx + 1
  
    } # for(p)
    
    # Write out the probs for the current chromosome.
    saveRDS(new_probs, file = file.path(out_dir, str_c('gen', gen, '_alleleprobs_chr', chr[c], '.rds')))

  } # for(c)

} # for(gen)

