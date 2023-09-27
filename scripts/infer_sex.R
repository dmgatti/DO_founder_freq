# Infer the sex of the mice from the sample intensities.
# Arguments:
# inten: data.frame containing sample id, X & Y intensities.
# Returns:
# data.frame containing sample id and inferred sex.
# XO:  chrX < 0.45  & chrY < 0.1
# XX:  chrX >= 0.45 & chrY < 0.1
# XY:  chrX < 0.47  & chrY >= 0.1
# XXY: chrX >= 0.47 & chrY > 0.2
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-21
################################################################################
infer_sex = function(inten) {

  retval = data.frame(id  = inten$id,
                      sex = NA)

  # Females: high chrX, low chrY.
  retval$sex[inten$chrX >= 0.45 & inten$chrY < 0.1] = 'F'

  # Males: low chrX, high chrY.
  retval$sex[inten$chrX < 0.47  & inten$chrY >= 0.1] = 'M'

  # XO females: low chrX, low chrY.
  retval$sex[inten$chrX < 0.45  & inten$chrY < 0.1] = 'XO'

  # XXY: high chrX, high chrY.
  retval$sex[inten$chrX >= 0.47 & inten$chrY > 0.2] = 'XXY'

  return(retval)

} # infer_sex()
