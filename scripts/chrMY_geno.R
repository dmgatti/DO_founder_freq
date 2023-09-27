################################################################################
# Given the CC/DO founder and sample genotypes, assign each sample to one of the
# founder mitochondrialor chr Y genotypes.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-02-23
################################################################################

# Assign mitochondrial genotypes.
# Genotypes can be 1-letter or 2-letter.
# Arguments:
# founders: character matrix containing genotypes. Markers in rows, samples
#           in columns, named as A:H. First column should be 'marker'.
# samples:  character matrix containing genotypes. Markers in rows, samples
#           in columns. First column should be 'marker'.
# We expect five groups in chr M:
# ABCD
# E
# F
# G
# H
#
# We expect six groups in Chr Y:
# A
# BCE
# D
# F
# G
# H
#
# Returns: data.frame containing three columns.
#   id: the sample ID provided in the samples argument.
#   mito: character string indicating the most probable mitochondrial 
#         genotype.
#   cor: the correlation of the most probable mitochondrial genotype with
#        the corresponding founder.
#
chrM_geno = function(founders, samples) {

  # Verify that the data contains the same markers. We'll trust that
  # they are mitochondrial, but perhaps adding markers to the arguments 
  # would be a good idea.
  stopifnot(nrow(founders) == nrow(samples))
  stopifnot(founders$marker == samples$marker)

  # Make substitution vector for mitochondrial groups.
  chr_groups = setNames(rep(c('ABCD', 'E', 'F', 'G', 'H'), c(4, 1, 1, 1, 1)),
                        LETTERS[1:8])

  # Convert founder data into a long numeric matrix with 8 columns,
  # one for each founder.
  f = founders
  rownames(f) = f$marker
  f  = data.frame(t(f[,-1]))
  dn = dimnames(f)
  f  = lapply(f, factor)
  f  = matrix(as.numeric(unlist(f)), ncol = length(f[[1]]), byrow = TRUE, 
              dimnames = list(NULL, dn[[1]]))

  s = samples
  rownames(s) = s$marker
  s  = data.frame(t(s[,-1]))
  dn = dimnames(s)
  s  = lapply(s, factor)
  s  = matrix(as.numeric(unlist(s)), ncol = length(s[[1]]), byrow = TRUE,
              dimnames = list(NULL, dn[[1]]))

  fg_cor  = cor(f, s, use = 'pairwise')
  max_cor = apply(fg_cor, 2, max)
  max_grp = apply(fg_cor, 2, which.max)
  max_grp = setNames(colnames(f)[max_grp], names(max_grp))

  stopifnot(names(max_cor) == names(max_grp))

  return(data.frame(id       = names(max_cor),
                    chrM     = max_grp,
                    chrM_cor = max_cor))

} # chrM_geno()



# Arguments:
# founders: character matrix containing genotypes. Markers in rows, samples
#           in columns, named as A:H. First column should be 'marker'.
# samples:  character matrix containing genotypes. Markers in rows, samples
#           in columns. First column should be 'marker'.
# sex:      data.frame containing 'id' column and 'sex' column.
chrY_geno = function(founders, samples, sex) {

  # Verify that the data contains the same markers. We'll trust that
  # they are mitochondrial, but perhaps adding markers to the arguments 
  # would be a good idea.
  stopifnot(nrow(founders) == nrow(samples))
  stopifnot(founders$marker == samples$marker)

  # Make substitution vector for Chr Y groups.
  chr_groups = setNames(rep(c('A', 'BCE', 'D', 'F', 'G', 'H'), c(1, 3, 1, 1, 1, 1)),
                        LETTERS[1:8])

  # Convert founder data into a long numeric matrix with 8 columns,
  # one for each founder.
  f = founders
  rownames(f) = f$marker
  f  = data.frame(t(f[,-1]))
  dn = dimnames(f)
  f  = lapply(f, factor)
  f  = matrix(as.numeric(unlist(f)), ncol = length(f[[1]]), byrow = TRUE, 
              dimnames = list(NULL, dn[[1]]))

  s = samples
  rownames(s) = s$marker
  s  = data.frame(t(s[,-1]))
  dn = dimnames(s)
  s  = lapply(s, factor)
  s  = matrix(as.numeric(unlist(s)), ncol = length(s[[1]]), byrow = TRUE,
              dimnames = list(NULL, dn[[1]]))

  fg_cor  = cor(f, s, use = 'pairwise')
  # Replace NAs with 0 so that max_grp isn't turned into a list.
  # This may make females look like they match with 'ABCD'.
  fg_cor[is.na(fg_cor)] = 0
  max_cor = apply(fg_cor, 2, max)
  max_grp = apply(fg_cor, 2, which.max)
  max_grp = setNames(colnames(f)[max_grp], names(max_grp))  

  stopifnot(names(max_cor) == names(max_grp))
  
  # Set females to NA.
  retval = data.frame(id       = names(max_cor),
                      chrY     = max_grp,
                      chrY_cor = max_cor) %>%
             full_join(select(sex, id, sex), by = 'id') %>%
             mutate(chrY     = if_else(sex == 'F' | sex == 'XO', NA, chrY),
                    chrY_cor = if_else(sex == 'F' | sex == 'XO', NA, chrY_cor)) %>%
             select(-sex)

  return(retval)

} # chrY_geno()

