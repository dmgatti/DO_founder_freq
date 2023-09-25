################################################################################
# Read in all of the FinalReport files and return the call rates and mean 
# intensities. Infer sex and fill in the sample metadata file.
# Run this first to have sex and call rate for downstream processes.
# We need to infer sex with all samples since single-sex studies are hard to
# infer by clustering alone. My goal is to get intensity thresholds to use
# globally.
# After running, it looks like these thresholds will owrk AFTER removing 
# samples with low call rates.
#
# XO:  chrX < 0.45  & chrY < 0.1
# XX:  chrX >= 0.45 & chrY < 0.1
# XY:  chrX < 0.47  & chrY >= 0.1
# XXY: chrX >= 0.47 & chrY > 0.2
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-21
################################################################################


###### LIBRARIES #####

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(qtl2convert))
suppressPackageStartupMessages(library(qtl2))

##### VARIABLES #####

# Get command line arguments.
args = commandArgs(trailingOnly = TRUE)
# args = c('/compsci/gedi/DO_founder_freq/data/neogen_dirs.txt', '/compsci/gedi/DO_founder_freq/data/')

print(args)

if (length(args) != 2) {

  print("usage: read_neogen <NEOGEN_DIR> <OUTPUT_DIR>")
  quit(status = 1)

} # if (length(args) != 2) 

# arg[1]: path to file containing top-level Neogen directories.
neogen_dirs = scan(file = args[1], 
                   what = 'character',
                   sep  = '\n')

# arg[2]: path to output directory. We will write the sample metadata file
# here and use this one from here forward.
output_dir  = args[2]

print(paste('INPUT:', neogen_dirs))

print(paste('OUTPUT:', output_dir))

# Top-level directory for the project.
base_dir = '/compsci/gedi/DO_founder_freq'

# Source directory for Neogen genotypes.
source_dir = '/gedi/resource/MUGA/raw_neogen_data/GigaMUGA'

# Sample metadata file.
metadata_file = file.path(source_dir, 'GigaMUGA.samples.20230918.csv')

# GigaMUGA GRCm39 marker file.
marker_file = '/compsci/gedi/DO_founder_freq/data/gm_uwisc_v4.csv'

neogen_dirs = file.path(source_dir, neogen_dirs)

##### FUNCTIONS #####

source('/compsci/gedi/DO_founder_freq/scripts/read_neogen.R')

# Read in the FinalReport files and return call rates and mean X/Y intensisites.
# This is a pre-screen to filter samples. We need all of the samples together
# to infer sex because we may have single-sex studies which make sex inferrence
# difficult.
# Arguments:
# neogen_dir: character string containing the full path to the top-level
#             project directory.
# metadata: sample metadata file with all samples.
# markers: data.frame of markers from U Wisc. 
# Returns:
# data.frame containing call rates and mean X/Y intensities.
scan_neogen = function(neogen_dir, metadata, markers) {

  results = read_neogen(neogen_dir, markers)
  
  # Filter to retain samples in metadata file.
  samples2keep = results$inten$id[results$inten$id %in% metadata$Unique.Sample.ID]
  
  results$geno  = results$geno[,c(1, which(colnames(results$geno) %in% samples2keep))] 
  results$inten = subset(results$inten, results$inten$id %in% samples2keep)
  
  # Get call rate.
  cr  = get_call_rate(results$geno)
  
  return(merge(cr, results$inten, by = 'id', sort = FALSE))

} # scan_neogen()

source('/compsci/gedi/DO_founder_freq/scripts/get_call_rate.R')
source('/compsci/gedi/DO_founder_freq/scripts/infer_sex.R')

##### MAIN #####

# Read in the sample metadata file.
metadata = read.csv(metadata_file)
  
# Read in the GigaMUGA GRCm39 markers.
markers = read.csv(marker_file)

cr_inten = NULL

# Go through each directory and get the call rates and chr X & Y 
# mean intensities.
for(curr_dir in neogen_dirs) {

  results = scan_neogen(curr_dir, metadata, markers)

  cr_inten = rbind(cr_inten, results)

  saveRDS(cr_inten, file = file.path(output_dir, 'callrate_inten.rds'))  

} # for(curr_dir)

# Save the data one last time.
saveRDS(cr_inten, file = file.path(output_dir, 'callrate_inten.rds'))






