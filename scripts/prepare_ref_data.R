################################################################################
# Gather an prepare DO GigaMUGA reference data.
# Create one file per chromosome for founder genotypes, genetic and physiscal 
# maps.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-23
################################################################################

library(qtl2convert)
library(tidyverse)

##### VARIABLES #####

# Path to reference data.
ref_dir = '/compsci/gedi/DO_founder_freq/data'

# Path to output directory.
out_dir = file.path(ref_dir, 'reference')

# Founder genotypes file, filtered to retain markers with known positions.
geno_file = file.path(ref_dir, 'GigaMUGA_founder_consensus_genotypes.csv')

# Marker file. GRCm39
marker_file = file.path(ref_dir, 'gm_uwisc_v4.csv')


##### MAIN #####

# Read in markers.
markers = readr::read_csv(marker_file, show_col_types = FALSE) %>%
            select(marker, chr, bp_grcm39, cM_cox, snp) %>%
            mutate(pos = bp_grcm39 * 1e-6)

# Read in genotypes.
geno = readr::read_csv(geno_file, show_col_types = FALSE) %>%
         rename(A = `A/J`,
                B = `C57BL/6J`,
                C = `129S1/SvImJ`,
                D = `NOD/ShiLtJ`,
                E = `NZO/HILtJ`,
                F = `CAST/EiJ`,
                G = `PWK/PhJ`,
                H = `WSB/EiJ`)

# Convert genotypes to A/H/B allele codes.
# qtl2convert::encode_geno() requires a character matrix of alleles.
alleles = markers %>%
            select(marker, snp) %>%
            separate_wider_position(snp, widths = c('A' = 1, 'B' = 1)) %>%
            column_to_rownames(var = 'marker') %>%
            as.matrix()

alleles = alleles[geno$marker,]

# Write out alleles.
# Writing out as an RDS because encode_geno() needs a matrix.
saveRDS(alleles, file = file.path(out_dir, 'allele_codes.rds'))

geno = geno %>%
         column_to_rownames(var = 'marker') %>%
         as.matrix()         
geno = qtl2convert::encode_geno(geno, allele_codes = alleles)  

geno = geno %>%
         as.data.frame() %>%
         rownames_to_column(var = 'marker')

# Join markers and genotypes.
geno_mkr = right_join(markers, geno, by = 'marker')

unique_chr = unique(geno_mkr$chr)

for(i in unique_chr) {

  data = filter(geno_mkr, chr == i)
  
  # Write genotypes.
  data %>%
    select(marker, A:H) %>%
    write_csv(file = file.path(out_dir, str_c('founder_geno_chr', i, '.csv')))
    
  # Write genetic map.
  data %>%
    select(marker, chr, pos = cM_cox) %>%
    write_csv(file = file.path(out_dir, str_c('gmap_chr', i, '.csv')))
  
  # Write physical map.
  data %>%
    select(marker, chr, pos) %>%
    write_csv(file = file.path(out_dir, str_c('pmap_chr', i, '.csv')))

} # for(i)

