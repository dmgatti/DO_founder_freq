################################################################################
# Prepare sample genotypes by filtering and coverting to A/H/B alleles.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-23
################################################################################

# Arguments:
# geno: data.frame containing genotypes. Markers in a column called 'marker'.
#       Markers in rows, samples in columns.
# alleles: character matrix containing allele codes for each marker.
# markers: data.frame containing marker positions.
# Returns:
#   list containing one element per chromosome. Each element is a data.frame
#     containing genotypes for one chromosome, encoded as A/H/B.
prepare_sample_geno = function(geno, alleles, markers) {

  # Convert genoypes to matrix for encode_geno
  geno = geno %>%
           column_to_rownames(var = 'marker') %>%
           as.matrix()

  # Synch up marker order in genotypes and alleles.
  geno = geno[rownames(alleles),]

  stopifnot(rownames(geno) == rownames(alleles))

  # Encode genotypes as A/H/B.
  geno = qtl2convert::encode_geno(geno = geno, allele_codes = alleles)

  # Select marker, chr, and samples.
  geno = geno %>%
           as.data.frame() %>%
           rownames_to_column(var = 'marker') %>%
           left_join(markers, by = 'marker') %>%
           select(-(bp_grcm39:pos))

  # Split samples by chromosome.
  geno_chr = factor(geno$chr)
  geno = geno[,colnames(geno) != 'chr']
  
  return(split(geno, geno_chr))

} # prepare_sample_genotypes()
