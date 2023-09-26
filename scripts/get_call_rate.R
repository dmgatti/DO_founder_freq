################################################################################
# Given a data.frame of genotypes, with no-call coded as '--', calculate the
# call-rate.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-26 
################################################################################

# Arguments:
# geno: data.frame containing genotypes a two-letter codes. First column is
#       'marker'. Markers in rows, samples in columns.
get_call_rate = functionn(geno) {

  geno = as.matrix(geno[,-1])
  cr   = colMeans(geno == '--', na.rm = TRUE) 

  return(data.frame(id = names(cr), call_rate = cr))
  
} # get_call_rate()
