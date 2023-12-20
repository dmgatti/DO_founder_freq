# Get the call rate for each sample.
# Arguments:
# geno: data.frame containing genotypes. Markers in rows, samples in columns.
# Returns:
# data.frame containing sample ID & call rate (0 to 1).
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-21
################################################################################
get_call_rate = function(geno) {

  geno = as.matrix(geno[,-1])

  cr   = colMeans(geno != '--', na.rm = TRUE)

  return(data.frame(id        = names(cr), 
                    call_rate = cr))

} # get_call_rate()

