# Infer the sex of the mice from the sample intensities.
# Arguments:
# x: data.frame containing X channel intensities. Markers in rows, samples
#    in columns. 'marker' in first column.
# y: data.frame containing Y channel intensities. Markers in rows, samples
#    in columns. 'marker' in first column.
# markers: data.frame containing marker positions. 
# Returns:
# data.frame containing sample id and inferred sex.
# XO:  chrX <  0.9 & chrY < 0.3
# XX:  chrX >= 0.9 & chrY < 0.3
# XY:  chrX <  1.0 & chrY >= 0.3
# XXY: chrX >= 1.0 & chrY > 0.3
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-21
################################################################################
infer_sex = function(x, y, markers) {

  # Get Chr X & Y markers.
  chrx = markers %>%
           select(marker, chr) %>%
           filter(chr == 'X') %>%
           pull(marker)

  chry = markers %>%
           select(marker, chr) %>%
           filter(chr == 'Y') %>%
           pull(marker)

  # Convert intensities into matrices.
  x = x %>%
        column_to_rownames(var = 'marker') %>%
        as.matrix()
  y = y %>%
        column_to_rownames(var = 'marker') %>%
        as.matrix()

  # Get Chr X and Y mean intensities for each sample.
  inten = data.frame(id   = colnames(x),
                     chrx = 0,
                     chry = 0,
                     sex  = '')

  inten$chrx = colMeans(x[chrx,] + y[chrx,], na.rm = TRUE)
  inten$chry = colMeans(x[chry,] + y[chry,], na.rm = TRUE)
  
  # Females: high chrX, low chrY.
  inten$sex[inten$chrx >= 0.9 & inten$chry < 0.3] = 'F'

  # Males: low chrX, high chrY.
  inten$sex[inten$chrx < 1.0  & inten$chry >= 0.3] = 'M'

  # XO females: low chrX, low chrY.
  inten$sex[inten$chrx < 0.9  & inten$chry < 0.3] = 'XO'

  # XXY: high chrX, high chrY.
  inten$sex[inten$chrx >= 1.0 & inten$chry > 0.3] = 'XXY'

  return(inten)

} # infer_sex()
