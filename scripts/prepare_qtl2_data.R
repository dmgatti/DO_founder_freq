################################################################################
# Prepare data for qtl2 to run calc_genoprob().
# We do not rrun calc_genoprob() here becaues of the wide variation in project
# sizes. We perform the following steps:
# 1. Read Neogen data for all directories in a project.
# 2. Calculate call rates and write them to the sample metadata file. Discard 
#    samples below 0.9.
# 3. Infer sex and Chr M & Y haplotypes, and write to sample metadata file.
# 4. Write out physical map, genetic map, founder genotypes, sample genotypes,
#    covariates, cross info, and teh JSON file.
# The next step is to run haplotype reconstruction with enough resources for
# the specific project.
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
# args = c('/gedi/resource/MUGA/raw_neogen_data/GigaMUGA/209_DO_Pack_Sleep', '/compsci/gedi/DO_founder_freq/data/genotypes/209_DO_Pack_Sleep_')

if (length(args) != 2) {

  print("usage: read_neogen <NEOGEN_DIR> <OUTPUT_DIR>")
  quit(status = 1)

} # if (length(args) != 2) 

# arg[1]: top-level project directory for Neogen files..
neogen_dir = args[1]

# arg[2]: path to output directory. This should have a filename prefix
# as well (i.e. /compsci/gedi/DO_founder_freq/genotypes/project1_).
output_dir  = args[2]

print(paste('INPUT:', neogen_dir))

print(paste('OUTPUT:', output_dir))

# Top-level directory for the project.
base_dir = '/compsci/gedi/DO_founder_freq'

# Data directory.
data_dir = file.path(base_dir, 'data')

# Reference data directory.
ref_dir = file.path(data_dir, 'reference')

# Output directory for qtl2 cross files.
qtl2_dir = file.path(data_dir, 'qtl2')

# Source directory for Neogen genotypes.
source_dir = '/gedi/resource/MUGA/raw_neogen_data/GigaMUGA'

# Sample metadata file.
# This is a copy of GigaMUGA.samples.20230918.csv which we modify
# by adding inferred sex, and Chr M & Y genotypes.
metadata_file = file.path(data_dir, 'gigamuga_sample_metadata.csv')

# GigaMUGa GRCm39 marker file.
marker_file = '/compsci/gedi/DO_founder_freq/data/gm_uwisc_v4.csv'

# Allele code file.
allele_file = file.path(ref_dir, 'allele_codes.rds')

##### FUNCTIONS #####

source('/compsci/gedi/DO_founder_freq/scripts/read_neogen.R')
source('/compsci/gedi/DO_founder_freq/scripts/get_call_rate.R')
source('/compsci/gedi/DO_founder_freq/scripts/infer_sex.R')
source('/compsci/gedi/DO_founder_freq/scripts/prepare_sample_geno.R')
source('/compsci/gedi/DO_founder_freq/scripts/chrMY_geno.R')
source('/compsci/gedi/DO_founder_freq/scripts/write_json.R')


##### MAIN #####

# Read in the sample metadata file.
metadata = read.csv(metadata_file)
  
# Subset the metadata to retain the rows for the curernt project.
curr_meta = subset(metadata, metadata$Project == basename(neogen_dir))

# Read in the GigaMUGA GRCm39 markers.
markers = readr::read_csv(marker_file, show_col_types = FALSE) %>%
            select(marker, chr, bp_grcm39, cM_cox, snp) %>%
            mutate(pos = bp_grcm39 * 1e-6)

# Read in the Neogen files and output genotypes and intensities.
result = read_neogen(neogen_dir, markers)
    
# Filter the samples to retain ones with call rates >= 0.9.
cr  = get_call_rate(result$geno)

# Infer sex from the Chr X & Y intensities.
sex = infer_sex(result$inten)

# Fill in metadata.
local_meta = subset(metadata, Project == basename(neogen_dir))
local_meta$Sex = sex$sex[match(local_meta$Unique.Sample.ID, sex$id)]

# Remove the current samples from metadata for now. 
# We will add them back in once we have sex, chr M & Y genotypes.
metadata = subset(metadata, Project != basename(neogen_dir))

# Read in alleles.
alleles = readRDS(allele_file)

# Prepare genotypes for qtl2.
geno = prepare_sample_geno(geno = result$geno, alleles = alleles, markers = markers)

for(i in seq_along(geno)) {

  write_csv(geno[[i]], file = file.path(qtl2_dir, str_c('geno_chr', names(geno)[i], '.csv')))

} # for(i)

# Write covar file.
covar = data.frame(id  = local_meta$Unique.Sample.ID,
                   sex = local_meta$Sex,
                   gen = local_meta$DO.Generation) %>%
         mutate(sex = if_else(sex == 'F',   'female', sex),
                sex = if_else(sex == 'M',   'male',   sex),
                sex = if_else(sex == 'XO',  'male',   sex),
                sex = if_else(sex == 'XXY', 'female', sex))

covar_file = file.path(qtl2_dir, 'covar.csv')
write_csv(covar, file = covar_file)

# Write phenotype file.
pheno = data.frame(id    = covar$id,
                   pheno = rnorm(nrow(covar)))

pheno_file = file.path(qtl2_dir, 'pheno.csv')
write_csv(pheno, file = pheno_file)

# Estimate Chr Y & M haplotypes.
# Chr M
founders = read_csv(file.path(data_dir, 'GigaMUGA_founder_consensus_genotypes_Mt.csv')) %>%
             rownames_to_column(var = 'marker') %>%
             as.data.frame()
samples  = result$geno %>%
             left_join(select(markers, marker, chr), by = 'marker') %>%
             filter(chr == 'M') %>%
             pivot_longer(cols = -marker) %>%
             mutate(value = str_sub(value, start = 1L, end = 1L)) %>%
             pivot_wider(names_from = name, values_from = value) %>%
             select(-chr) %>%
             as.data.frame()
samples = samples[match(founders$marker, samples$marker),]

chrM = chrMY_geno(founders, samples, chr = 'M') %>%
         rename(chrM     = chr_group,
                chrM_cor = cor)

# Chr Y
founders = read_csv(file.path(data_dir, 'GigaMUGA_founder_consensus_genotypes_Y.csv')) %>%
             rownames_to_column(var = 'marker') %>%
             as.data.frame()
samples  = result$geno %>%
             left_join(select(markers, marker, chr), by = 'marker') %>%
             filter(chr == 'Y') %>%
             pivot_longer(cols = -marker) %>%
             mutate(value = str_sub(value, start = 1L, end = 1L)) %>%
             pivot_wider(names_from = name, values_from = value) %>%
             select(-chr) %>%
             as.data.frame()
samples = samples[match(founders$marker, samples$marker),]

chrY = chrMY_geno(founders, samples, chr = 'Y') %>%
         rename(chrY     = chr_group,
                chrY_cor = cor)

# Add the Chr M & Y genotypes to the sample metadata.
local_meta = local_meta %>%
               select(-(chrM:chrY_cor)) %>%
               left_join(chrM, by = c('Unique.Sample.ID' = 'id')) %>%
               left_join(chrY, by = c('Unique.Sample.ID' = 'id'))

if(!'chrM' %in% colnames(metadata)) {

  metadata = cbind(metadata, 
                   data.frame(chrM     = NA_character_,
                              chrM_cor = NA_real_,
                              chrY     = NA_character_,
                              chrY_cor = NA_real_))

} # if(!'chrM' %in% colnames(metadata))

metadata = rbind(metadata, local_meta)
write_csv(metadata, metadata_file)

# Write JSON files.
for(i in c(1:19, 'X')) {

  geno_file    = str_c('geno_chr', i, '.csv')
  founder_file = file.path('reference', str_c('founder_geno_chr', i, '.csv'))
  gmap_file    = file.path('reference', str_c('gmap_chr', i, '.csv'))
  pmap_file    = file.path('reference', str_c('pmap_chr', i, '.csv'))
  
  write_json(out_dir = qtl2_dir, chr = i, geno_file, founder_file, gmap_file, 
             pmap_file, basename(pheno_file), basename(covar_file))

} # for(i in c(1:19, 'X'))


  # Run qtl2 and calculate genoprobs.
  
  # Count & locate crossovers. (Need probs & cross object)
  
  # QC based on crossover count.
  
  # Clean up. (genoprobs, genotypes, intensities)

