################################################################################
# Use argyle to get the log2-ratio and B-allele frequency for each sample
# across the genome.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-12-20
################################################################################

library(argyle)

##### VARIABLES #####

# Base project directory.
base_dir = '/compsci/gedi/DO_founder_freq'

# Data directory.
data_dir = file.path(base_dir, 'data')

# Results directory.
results_dir = file.path(base_dir, 'results')

# Metadata file.
metadata_file = file.path(data_dir, 'gigamuga_sample_metadata.csv')

# Output directory for log2-ratio and B-allele frequency.
out_dir = file.path(results_dir, 'lrr_b_freq')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Neogen top-level directory.
neogen_dir = '/gedi/resource/MUGA/raw_neogen_data/GigaMUGA'

# Marker directory.
marker_dir = '/compsci/gedi/reference'

# Marker file.
marker_file = file.path(marker_dir, 'gm_uwisc_v4.csv')

# UNC Gigamuga cluster file.
cluster_file = file.path(marker_dir, 'clusters.gigamuga.Rdata') 

# Project file.
project_file = file.path(data_dir, 'neogen_dirs.txt')


##### MAIN #####

# Read in the project file.
projects = scan(project_file, what = 'char', sep = '\n')

# Read in the sample metadata file.
metadata = read.csv(metadata_file)

# Read in the markers and reformat for argyle.
markers = read.csv(marker_file)
markers = markers[,c('chr', 'marker', 'cM_cox', 'bp_grcm39', 'snp', 'snp')]
dimnames(markers) = list(markers$marker, c('chr', 'marker', 'cM', 'pos', 'A1', 'A2'))
markers$A1 = substring(markers$A1, 1L, 1L)
markers$A2 = substring(markers$A2, 2L, 2L)

# Read in UNC Gigamuga clusters. (Object name = 'clusters')
load(cluster_file)

# Go through each project and read in all of the genotype data using
# argyle. Extract the log2-ratio (LRR) and B-allele frequency (BAF).
for(project in projects) {

  print(project)
  
  project_dir = file.path(neogen_dir, project)
  
  # Get the Neogen directories in this project.
  neogen_dirs = list.dirs(path = project_dir, full.names = TRUE, recursive = FALSE)
  
  # Genotype object for all samples.
  all_gty = NULL
  
  for(d in neogen_dirs) {
  
    # Read in current Neogen files.
    g = read.beadstudio(prefix  = basename(d),
                        snps    = markers, 
                        in.path = d)
    
    # Add sample prefixes.
    # We have to do surgery on the object. Is there an easier way?
    colnames(g) = paste(basename(d), colnames(g), sep = ';')
    ped         = attr(g, 'ped')
    ped$fid     = colnames(g)
    ped$iid     = colnames(g)
    rownames(ped)  = colnames(g)
    attr(g, 'ped') = ped

    attr(g, 'filter.samples') = colnames(g)
    
    intensities             = attr(g, 'intensity')
    colnames(intensities$x) = colnames(g)
    colnames(intensities$y) = colnames(g)
    attr(g, 'intensity')    = intensities
    
    # We're keeping all markers for CNV analysis.
    # But we are filtering out samples with call rates < 90%.
    g = run.sample.qc(gty = g, max.N = round(0.1 * nrow(g)), apply = TRUE)

    # Combine genotypes with all gentoypes.
    if(is.null(all_gty)) {
      all_gty = g
    } else {
      all_gty = cbind.genotypes(all_gty, g)
    } # else

    # Clean up.
    rm(ped, intensities, g)
    gc()

  } # for(d)
  
  # Normalize the intensities.
  all_gty = tQN(gty = all_gty, clusters = clusters)
  
  # Keep only the samples in the sample metadata.
  all_gty = all_gty[,colnames(all_gty) %in% metadata$Unique.Sample.ID]
  
  # Write out LRR and BAF.
  saveRDS(attr(all_gty, 'lrr'), file = file.path(out_dir, paste0(project, '_lrr.rds')))
  saveRDS(attr(all_gty, 'baf'), file = file.path(out_dir, paste0(project, '_baf.rds')))
  
  # Clean up.
  rm(all_gty)
  gc()
  
} # for(project)


