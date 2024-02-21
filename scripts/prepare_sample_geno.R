################################################################################
# Prepare sample genotypes by filtering and coverting to A/H/B alleles.
# GigaMUGA version
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-23
################################################################################

# Arguments:
# samples: data.frame containing genotypes. Markers in a column called 'marker'.
#       Markers in rows, samples in columns.
# founders: data.frame containing founder genotypes. Markers in a column called 
#           'marker'. Markers in rows, samples in columns.
# alleles: character matrix containing allele codes for each marker.
# markers: data.frame containing marker positions.
# Returns:
#   list containing one element per chromosome. Each element is a data.frame
#     containing genotypes for one chromosome, encoded as A/H/B.
prepare_sample_geno = function(samples, founders, alleles, markers) {

  # Combine founders and geno.
  all_geno = full_join(founders, samples, by = 'marker')

  # Convert genoypes to matrix for encode_geno
  all_geno = all_geno %>%
               column_to_rownames(var = 'marker') %>%
               as.matrix()

  # Synch up marker order in genotypes and alleles.
  all_geno = all_geno[rownames(alleles),]

  stopifnot(rownames(all_geno) == rownames(alleles))

  # Encode genotypes as A/H/B.
  all_geno = qtl2convert::encode_geno(geno = all_geno, allele_codes = alleles)
  markers  = subset(markers, marker %in% rownames(all_geno))
  all_geno = all_geno[markers$marker,]

  # Separate founders and samples.
  f = all_geno[,1:(ncol(founders) - 1)]
  s = all_geno[,-(1:(ncol(founders) - 1))]

  # Organize founders by chromosome.
  # Select marker, chr, and samples.
  f = f %>%
        as.data.frame() %>%
        rownames_to_column(var = 'marker') %>%
        left_join(markers, by = 'marker') %>%
        select(-(bp_grcm39:pos))
  
  # Split samples by chromosome.
  f_chr = factor(f$chr)
  f     = f[,colnames(f) != 'chr']
  f     = split(f, f_chr)

  # Organize samples by chromosome.
  # Select marker, chr, and samples.
  s = s %>%
        as.data.frame() %>%
        rownames_to_column(var = 'marker') %>%
        left_join(markers, by = 'marker') %>%
        select(-(bp_grcm39:pos))
  
  # Split samples by chromosome.
  s_chr = factor(s$chr)
  s     = s[,colnames(s) != 'chr']
  s     = split(s, s_chr)

  return(list(founders = f, samples = s))

} # prepare_sample_genotypes()
