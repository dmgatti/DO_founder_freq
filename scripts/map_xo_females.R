################################################################################
# Gather the females and map XO as a binary trait.
# Note that you need to use the -B /compsci argument in singularity.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-12-07
################################################################################

library(qtl2convert)
library(qtl2)

##### VARIABLES #####

# Base project directory.
base_dir = '/compsci/gedi/DO_founder_freq'

# Data directory.
data_dir = file.path(base_dir, 'data')

# Sample metadata directory.
meta_dir = file.path(data_dir, 'metadata')

# Results directory.
results_dir = file.path(base_dir, 'results')

# Genoprobs directory.
probs_dir = file.path(results_dir, 'genoprobs')

# Output directory.
out_dir = file.path(results_dir, 'map_xo_females')

# Gigamuga marker file.
gm_marker_file = '/compsci/gedi/reference/gm_uwisc_v4.csv'

# qtl2 SNP query function.
snp_fxn = create_variant_query_func('/compsci/gedi/reference/fv.2021.snps.db3', 
                                    id_field = 'variants_id')

# qtl2 gene query function.
gene_fxn = create_gene_query_func('/compsci/gedi/reference/fv.2021.snps.db3', 
                                  chr_field   = 'chromosome', 
                                  name_field  = 'symbol',
                                  start_field = 'start_position', 
                                  stop_field  = 'end_position')

##### MAIN #####

# Read in the marker file and create a physical map.
markers     = read.csv(gm_marker_file)
markers$pos = markers$bp_grcm39 * 1e-6
map         = map_df_to_list(markers, pos_column = 'pos')

# Get the sample metadata files.
metadata_files = dir(meta_dir, pattern = 'csv$', full.names = TRUE)

# Extract the project names.
projects = gsub(paste0('^', meta_dir, '/|_metadata.csv$'), '', metadata_files)
projects = projects[projects != 'Luanne_Peters']

# Read in all of the metadata and retain the females.
metadata = lapply(metadata_files, read.csv)
metadata = do.call(rbind, metadata)
metadata = subset(metadata, Sex %in% c('F', 'XO'))

# How many females and XOs do we have?
table(metadata$Sex)

# For each project, read in the genoprobs and retain the females.
# Then write out each chromosome.
projects = projects[projects %in% metadata$Project]

# Create the binary phenotype.
pheno = matrix(as.numeric(metadata$Sex == 'XO'), ncol = 1,
               dimnames = list(metadata$Unique.Sample.ID, 'XO'))

saveRDS(pheno, file = file.path(out_dir, 'pheno_xo.rds'))

for(chr in c(1:19, 'X')) {

  print(paste('CHR:', chr))

  # Read in the first probs file to get the probs array dimensions.
  probs_file = dir(probs_dir, pattern = paste0(projects[1], '_alleleprobs_chr', chr, '.rds'), 
                   full.names = TRUE)
  probs = readRDS(probs_file)[[1]]

  all_probs = array(0, dim = c(nrow(metadata), 8, dim(probs)[3]), 
                    dimnames = list(metadata$Unique.Sample.ID, LETTERS[1:8], dimnames(probs)[[3]]))

  for(prj in projects) {

    print(prj)

    # Read in the probs.
    probs_file = dir(probs_dir, pattern = paste0(prj, '_alleleprobs_chr', chr, '.rds'), 
                     full.names = TRUE)
    probs = readRDS(probs_file)[[1]]
    
    # Subset to retain the females (in metadata)
    probs = probs[rownames(probs) %in% metadata$Unique.Sample.ID,,]
    
    # Place the samples in the correct locations in all_probs.
    samples = rownames(probs)
    all_probs[samples,,] = probs

  } # for(prj)

  # Map XO on this chromosome.
  all_probs = list(all_probs)
  names(all_probs) = chr
  class(all_probs) = c('calc_genoprob', 'list')
  attr(all_probs, 'crosstype') = 'do'
  # Setting is_x_chr caused trouble in scan1snps().
  #attr(all_probs, 'is_x_chr')  = FALSE
  attr(all_probs, "alleles")   = LETTERS[1:8]
  attr(all_probs, 'alleleprobs') = TRUE
  
  saveRDS(all_probs, file = file.path(out_dir, paste0('chr', chr, '_alleleprobs_female.rds')))
  
  # LOD
  lod = scan1(genoprobs = all_probs, pheno = pheno, model = 'binary')
  map[[chr]] = map[[chr]][dimnames(all_probs[[1]])[[3]]]
  
  saveRDS(map, file = file.path(out_dir, paste0('XO_map_chr', chr, '.rds')))
  saveRDS(lod, file = file.path(out_dir, paste0('XO_lod_chr', chr, '.rds')))
  
  png(file.path(out_dir, paste0('XO_lod_chr', chr, '.png')), width = 1000,
      height = 800, res = 128)
  plot(lod, map, main = paste('XO, chr', chr))
  dev.off()

  # Founder allele effects.
  eff = scan1coef(genoprobs = all_probs, pheno = pheno, model = 'binary')
  saveRDS(eff, file = file.path(out_dir, paste0('XO_coef_chr', chr, '.rds')))

  png(file.path(out_dir, paste0('XO_coef_chr', chr, '.png')), width = 1000,
      height = 800, res = 128)
  plot_coefCC(eff, map, scan1_output = lod, main = paste('XO, chr', chr))
  dev.off()
  
  # Find peaks above a threshold.
  peaks = find_peaks(lod, map, threshold = 7, prob = 0.95, peakdrop = 3)
  
  # Association Mapping
  if(nrow(peaks) > 0) { 
  
    for(i in 1:nrow(peaks)) {
    
      start = peaks$pos[i] - 2
      end   = peaks$pos[i] + 2
  
      assoc = scan1snps(genoprobs = all_probs, map = map, pheno = pheno, 
                        model = 'binary', query_func = snp_fxn, chr = chr,
                        start = start, end = end, keep_all_snps = TRUE)
                    
      genes = gene_fxn(chr, start = start, end = end)
                    
      png(file.path(out_dir, paste0('XO_assoc_chr', chr, '_', round(start), '_', 
                    round(end), '.png')), 
          width = 1000, height = 800, res = 128)
      colnames(assoc$snpinfo)[1] = 'snp_id'
      plot_snpasso(scan1output = assoc$lod, snpinfo = assoc$snpinfo, 
                   genes = genes, drop_hilit = 1, minlod = 1, sdp_panel = TRUE, 
                   main = paste('XO, chr', chr))
      dev.off()

    } # for(i)

  } # if(nrow(peaks > 0))
  
} # for(chr)


# Permutations.
pheno = readRDS(file = file.path(out_dir, 'pheno_xo.rds'))

for(chr in c(1:19, 'X')) {

  print(chr)
  t1 = proc.time()[3]
  
  probs = readRDS(file = file.path(out_dir, paste0('chr', chr, '_alleleprobs_female.rds')))

  perms = scan1perm(probs, pheno, model = 'binary', n_perm = 10, cores = 5)
  
  saveRDS(perms, file = file.path(out_dir, paste0('perms_chr', chr, '.rds')))

  print(proc.time()[3] - t1)

} # for(chr)


