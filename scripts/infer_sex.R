################################################################################
#
################################################################################

# female: chrX > 0.47  & chrY < 0.15
# male:   chrX <= 0.47 & chrY >= 0.15
# XO:     chrX < 0.47  & chrY < 0.15
# XXY:    chrX > 0.47  & chrY > 0.25

# Arguments:
# intensities: data.frame containing three columns: 
#              id: mouse ID,
#              chrX: mean X chromosome intensity,
#              chrY: mean Y chromosome intensity.
# Returns:
# data.frame containing mouse ID and sex call 
infer_sex = function(intensities) {

  retval = data.frame(id  = intensities$id, 
                      sex = '')

  retval[intensities$chrX > 0.47  & intensities$chrY < 0.15] = 'female'
  retval[intensities$chrX <= 0.47 & intensities$chrY >= 0.15] = 'male'
  retval[intensities$chrX < 0.47  & intensities$chrY < 0.15] = 'XO'
  retval[intensities$chrX > 0.47  & intensities$chrY > 0.15] = 'XXY'

  return(retval)

} # infer_sex()

