################################################################################
# Gather all of the data required to produce a qtl2 cross object.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-26
################################################################################

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(qtl2convert))
suppressPackageStartupMessages(library(qtl2))


##### VARIABLES #####

# Argument is the name of the project directory.
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 1) {

  print('usage: prepare_qtl2_data.R <PROJECT>')
  quit(1)

} # if(length(args) != 1)

project = args[1]

# Source directory where Neogen files are stored.
source_dir = file.path('', project)

##### MAIN #####

# Read sample metadata.

# Read GigaMUGA markers.

# Read Neogen files.

# Filter samples by call rate.

# Infer sex.

# Get Chr M genotypes.

# Get Chr Y genotypes.

# Save sample metadata.

# Write covar file.

# Write pheno file.

# Write JSON files.


