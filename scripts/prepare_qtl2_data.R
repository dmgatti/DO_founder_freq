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
# args = 'DO_Arsenic'

if (length(args) != 1) {

  print("usage: read_neogen <NEOGEN_DIR>")
  quit(status = 1)

} # if (length(args) != 1) 

# arg[1]: top-level project directory for Neogen files.
project = args[1]

print(paste('project:', project))

# Top-level directory for the project.
base_dir = '/compsci/gedi/DO_founder_freq'

# Data directory.
data_dir = file.path(base_dir, 'data')

# Reference data directory.
ref_dir = file.path(data_dir, 'reference')

# Output directory for qtl2 cross files.
qtl2_dir = file.path(data_dir, 'qtl2')

# Source directory for Neogen genotypes.
source_dir = file.path('/gedi/resource/MUGA/raw_neogen_data/GigaMUGA', project)

# Output directory for genotypes.
output_geno_dir = file.path(data_dir, 'genotypes')

# Output directory for intensities.
output_inten_dir = file.path(data_dir, 'intensities')

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
metadata = subset(metadata, metadata$Project == project)

# Read in the GigaMUGA GRCm39 markers.
markers = readr::read_csv(marker_file, show_col_types = FALSE) %>%
            select(marker, chr, bp_grcm39, cM_cox, snp) %>%
            mutate(pos = bp_grcm39 * 1e-6)

# Read in the Neogen files and output genotypes and intensities.
neogen = read_neogen(source_dir)

geno    = neogen$geno
inten_x = neogen$inten_x
inten_y = neogen$inten_y

rm(neogen)

#metadata$Unique.Sample.ID[!metadata$Unique.Sample.ID %in% colnames(geno)]

# Only keep samples that are in the metadata.
geno    = geno[,   c('marker', metadata$Unique.Sample.ID)]
inten_x = inten_x[,c('marker', metadata$Unique.Sample.ID)]
inten_y = inten_y[,c('marker', metadata$Unique.Sample.ID)]

# Get the genotype call rates.
cr  = get_call_rate(geno)

# Retain samples with call rates >= 0.9.
cr = cr %>%
       filter(call_rate >= 0.9) %>%
       separate_wider_delim(id, delim = ';', names = c('dir', 'sample'),
                            cols_remove = FALSE)

metadata = metadata[metadata$Unique.Sample.ID %in% cr$id,]
geno     = geno[,   c('marker', cr$id)]
inten_x  = inten_x[,c('marker', cr$id)]
inten_y  = inten_y[,c('marker', cr$id)]

stopifnot(metadata$Unique.Sample.ID == cr$id)
stopifnot(colnames(geno)[-1]        == cr$id)
stopifnot(colnames(inten_x)[-1]     == cr$id)
stopifnot(colnames(inten_y)[-1]     == cr$id)

# Look for duplicate samples. If we find any, retain the one with the
# highest call rate.
dupl = cr$sample[which(duplicated(cr$sample))]
dupl = sort(unique(dupl))

if(length(dupl) > 0) {

  for(s in dupl) {
  
    # Find the sample with the highest call rate.
    wh   = which(cr$sample == s)
    keep = wh[which.max(cr$call_rate[wh])]
    
    # Remove the duplicate samples with lower call rates.
    rem = wh[wh != keep]

    metadata   = metadata[-rem,]
    # marker is in column 1, so add 1 to rem.
    geno    = geno[,-(rem + 1)]
    inten_x = inten_x[,-(rem + 1)]
    inten_y = inten_y[,-(rem + 1)]
    cr      = cr[-rem,]

  } # for(s)

} # if(length(dupl) > 0)

stopifnot(!duplicated(cr$id))

# Save genotypes and intensities.
saveRDS(geno,    file = file.path(output_geno_dir,  str_c(project, '_geno.rds')))
saveRDS(inten_x, file = file.path(output_inten_dir, str_c(project, '_x.rds')))
saveRDS(inten_y, file = file.path(output_inten_dir, str_c(project, '_y.rds')))

# Infer sex from the Chr X & Y intensities.
sex = infer_sex(inten_x, inten_y, markers)

# Fill in metadata.
metadata$Sex       = sex$sex[match(metadata$Unique.Sample.ID,      sex$id)]
metadata$call_rate = cr$call_rate[match(metadata$Unique.Sample.ID, cr$id)]

# Estimate Chr Y & M haplotypes.
# Chr M
founders = read_csv(file.path(data_dir, 'GigaMUGA_founder_consensus_genotypes_Mt.csv')) %>%
             rownames_to_column(var = 'marker') %>%
             as.data.frame()
samples  = geno %>%
             left_join(select(markers, marker, chr), by = 'marker') %>%
             filter(chr == 'M') %>%
             pivot_longer(cols = -marker) %>%
             mutate(value = str_sub(value, start = 1L, end = 1L)) %>%
             pivot_wider(names_from = name, values_from = value) %>%
             select(-chr) %>%
             as.data.frame()
samples = samples[match(founders$marker, samples$marker),]

chrM = chrM_geno(founders, samples)

# Chr Y
founders = read_csv(file.path(data_dir, 'GigaMUGA_founder_consensus_genotypes_Y.csv')) %>%
             rownames_to_column(var = 'marker') %>%
             as.data.frame()
samples  = geno %>%
             left_join(select(markers, marker, chr), by = 'marker') %>%
             filter(chr == 'Y') %>%
             pivot_longer(cols = -marker) %>%
             mutate(value = str_sub(value, start = 1L, end = 1L)) %>%
             pivot_wider(names_from = name, values_from = value) %>%
             select(-chr) %>%
             as.data.frame()
samples = samples[match(founders$marker, samples$marker),]

chrY = chrY_geno(founders, samples, sex = sex)

# Add the Chr M & Y genotypes to the sample metadata.
metadata = metadata %>%
               left_join(chrM, by = c('Unique.Sample.ID' = 'id')) %>%
               left_join(chrY, by = c('Unique.Sample.ID' = 'id'))

# Write out sample metadata.
write_csv(metadata, file = file.path(data_dir, 'metadata', 
          str_c(project, '_metadata.csv')))

# Read in alleles.
alleles = readRDS(allele_file)

# Prepare genotypes for qtl2.
geno = prepare_sample_geno(geno = geno, alleles = alleles, markers = markers)

for(i in seq_along(geno)) {

  write_csv(geno[[i]], file = file.path(qtl2_dir, str_c('geno_chr', names(geno)[i], '.csv')))

} # for(i)

# Write covar file.
covar = data.frame(id  = metadata$Unique.Sample.ID,
                   sex = metadata$Sex,
                   gen = metadata$DO.Generation) %>%
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

