################################################################################
# Read in the variant calling in the 48 DO samples that were parrt of the
# lcGBS SSIF. Remove Indels and filter out SNPs that are in the Sanger VCF.
# Then look at coverage and allele frequencies in the remaining samples.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-11-30 
################################################################################

##### LIBRARIES #####

library(VariantAnnotation)

#### VARIABLES #####

# Base project directory.
base_dir = '/compsci/gedi/DO_founder_freq'

# Results directory.
results_dir = file.path(base_dir, 'results')

# DO Variants directory.
var_dir = file.path(results_dir, 'do_variants')

# Sanger SNP file.
sanger_vcf_file = '/compsci/gedi/reference/mgp_REL2021_snps.vcf.gz'

# DO called variants from combined VCF.
comb_vcf_files = dir(var_dir, pattern = '_combined.vcf.gz$',
                full.names = TRUE)

# DO called variants from uncombined VCFs.
sample_vcf_files = dir(var_dir, pattern = '^do_variants_chr([0-9]+|X).vcf.gz$',
                       full.names = TRUE)

##### MAIN #####

# Read in each combined VCF and remove Indels and SNPs that are in the 
# Sanger VCF.

# Get the DO founder strain names.
hdr     = scanVcfHeader(sanger_vcf_file)
strains = samples(hdr)[c(4, 2, 34, 37, 19, 40, 51)]

for(i in seq_along(comb_vcf_files)) {

  # Get the chromosome for the current file.
  chr   = sub('chr', '', strsplit(basename(comb_vcf_files[i]), '_')[[1]][3])

  print(chr)

  # Read in the Sanger VCF.
  param   = ScanVcfParam(samples = strains, 
                         which   = GRanges(chr, IRanges(0, 200e6)))
  sanger_vcf = readVcf(sanger_vcf_file, genome = 'GRCm39', param = param)
  gt = geno(sanger_vcf)$GT
  # Retain polymorphic SNPs of any kind.
  gt = cbind(C57BL_6J = '0/0', gt)

  sanger_vcf = subset(sanger_vcf, rowMeans(gt == '0/0') < 1.0)

  # Read in the combined VCF.
  param    = ScanVcfParam(which   = GRanges(chr, IRanges(0, 200e6)))
  comb_vcf = readVcf(comb_vcf_files[i], genome = 'GRCm39', param = param) 
  
  # Remove Indels.
  comb_vcf = subset(comb_vcf, info(comb_vcf)$INDEL == FALSE)

  # Subset the de novo VCF to retain only SNPs that are NOT in the Sanger VCF.
  ol       = findOverlaps(comb_vcf, sanger_vcf)
  keep     = setdiff(1:nrow(comb_vcf), queryHits(ol))
  comb_vcf = comb_vcf[keep,]
  
  # Get the genotypes.
  gt              = geno(comb_vcf)$GT
  colnames(gt)[1] = 'geno'

  # Reformat the strand-specific read depths.
  dp4           = do.call(rbind, info(comb_vcf)$DP4)
  colnames(dp4) = c('ref_for', 'ref_rev', 'alt_for', 'alt_rev')
  
  # Assemble results data.frame.
  result    = data.frame(info(comb_vcf)[,-16], dp4, gt)
  result$AC = unlist(sapply(result$AC, paste0, collapse = ';'))
  
  pos_info = matrix(unlist(strsplit(rownames(result), '\\:|_|/')), ncol = 4, byrow = TRUE)
  result   = data.frame(chr = pos_info[,1], pos = as.numeric(pos_info[,2]), 
                        ref = pos_info[,3], alt = pos_info[,4], result)
  rm(pos_info)

  # Read in the sample VCF.
  param      = ScanVcfParam(which   = GRanges(chr, IRanges(0, 200e6)))
  sample_vcf = readVcf(sample_vcf_files[i], genome = 'GRCm39', param = param) 
  
  # Remove Indels.
  sample_vcf = subset(sample_vcf, info(sample_vcf)$INDEL == FALSE)
  
  # Filter the sample VCF to retain only SNPs in the filtered combined VCF.
  sample_vcf = subsetByOverlaps(sample_vcf, comb_vcf)
  comb_vcf   = subsetByOverlaps(comb_vcf,   sample_vcf)
  
  stopifnot(length(sample_vcf) == length(comb_vcf))
  stopifnot(rownames(rowRanges(sample_vcf)) == rownames(rowRanges(comb_vcf)))
  
  # Get the genotypes.
  gt = geno(sample_vcf)$GT
  colnames(gt) = gsub('^/compsci/gedi/DO_founder_freq/data/bams/plexWell-|_GT23-[0-9]+_[ACGT]+-[ACGT]+_S[0-9]+_L[0-9]+_sorted.bam$', '', colnames(gt))
  
  # Count the allele frequency across samples.
  unique_alleles = sort(unique(as.vector(gt)))
  af             = lapply(data.frame(t(gt)), factor, levels = unique_alleles)
  af             = t(sapply(af, table))
  rownames(af)   = rownames(gt)
  
  # Merge in allele frequencies.
  result = merge(result, af, by = 'row.names')
  rownames(result) = result$Row.names
  result = result[,-1]
  result = merge(result, gt, by = 'row.names')
  rownames(result) = result$Row.names
  result = result[,-1]
  
  # Write out the results.
  write.csv(result, file = file.path(results_dir, paste0('do_snps_chr', chr, '.csv')),
            quote = FALSE)

} # for(i)


