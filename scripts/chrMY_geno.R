################################################################################
# Given the CC/DO founder and sample genotypes, assign each sample to one of the
# founder mitochondrial or chr Y genotypes.
# GigaMUGA version.
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

  # List the Chr M haplogroups.
  chrM_haplogroups = c('ABCD', 'E', 'F', 'G', 'H')

  # Make substitution vector for mitochondrial groups.
  chr_groups = setNames(rep(chrM_haplogroups, c(4, 1, 1, 1, 1)),
                        LETTERS[1:8])
                        
  # Create a founder data.frame contining only the haplogroups.
  f = founders
  colnames(f)[-1] = chr_groups
  f = f[,unique(colnames(f))]
  
  # Combine the founder and sample data and convert genotypes to numbers.
  comb = merge(f, samples, by = 'marker', sort = FALSE)

  # Transpose that data to put markers in columns.
  rownames(comb) = comb$marker
  comb  = data.frame(t(comb[,-1]))
  
  # Get the dimnames to put on the combined matrix.
  dn    = dimnames(comb)
  
  # Factor the genotypes and convert to numbers. 
  comb  = lapply(comb, factor)
  comb  = lapply(comb, as.numeric)
  comb  = matrix(unlist(comb), nrow = length(comb), ncol = length(comb[[1]]), 
                 byrow = TRUE, dimnames = list(dn[[2]], dn[[1]]))

  # Split up the founder haplogroups and samples.
  f = comb[,1:length(chrM_haplogroups)]
  s = comb[,-(1:length(chrM_haplogroups))]
  rm(comb) 

  stopifnot(rownames(f) == rownames(s))

  # Get the correlation between all haplogrups and all samples.
  fg_cor  = cor(f, s, use = 'pairwise')
  
  # Select the founder haplogroup with the highest correlation to each sample.
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
  
  # If we have only females, we can return now. 
  if(all(unique(sex$sex) %in% c('F', 'XO'))) {
    retval = data.frame(id       = sex$id,
                        chrY     = rep(NA, nrow(sex)),
                        chrY_cor = rep(NA, nrow(sex)))

    return(retval)
  
  } # if(all(unique(sex$sex) %in% c('F', 'XO')))
  
  # List the Chr Y haplogroups.
  chrY_haplogroups = c('A', 'BCE', 'D', 'F', 'G', 'H')

  # Make substitution vector for Chr Y groups.
  chr_groups = setNames(c('A', 'BCE', 'BCE','D', 'BCE','F', 'G', 'H'),
                        LETTERS[1:8])

  # Create a founder data.frame contining only the haplogroups.
  f = founders
  colnames(f)[-1] = chr_groups
  f = f[,unique(colnames(f))]
  
  # Combine the founder and sample data and convert genotypes to numbers.
  comb = merge(f, samples, by = 'marker', sort = FALSE)

  # Transpose that data to put markers in columns.
  rownames(comb) = comb$marker
  comb  = data.frame(t(comb[,-1]))
  
  # Get the dimnames to put on the combined matrix.
  dn    = dimnames(comb)
  
  # Factor the genotypes and convert to numbers. 
  comb  = lapply(comb, factor)
  comb  = lapply(comb, as.numeric)
  comb  = matrix(unlist(comb), nrow = length(comb), ncol = length(comb[[1]]), 
                 byrow = TRUE, dimnames = list(dn[[2]], dn[[1]]))

  # Split up the founder haplogroups and samples.
  f = comb[,1:length(chrY_haplogroups)]
  s = comb[,-(1:length(chrY_haplogroups))]
  rm(comb) 

  stopifnot(rownames(f) == rownames(s))

  # Get the correlation between all haplogroups and all samples.
  fg_cor  = cor(f, s, use = 'pairwise')
  fg_cor[is.na(fg_cor)] = 0
  
  # Select the founder haplogroup with the highest correlation to each sample.
  max_cor = apply(fg_cor, 2, max)
  max_grp = apply(fg_cor, 2, which.max, simplify = TRUE)
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

