################################################################################
# Given a list of allele or genoprobs filenames, combine them into a single
# qtl2-style genoprobs object. 
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-10-13
################################################################################

require(qtl2)

# Gather multiple allele or genoprobs objects together in one qtl2-style object.
# You should combine either allele probs OR genoprobs, not both at the same 
# time.
# Arguments:
# file_list: string vector containing the full path to each of the allele or
#            genoprobs files. 
# Returns:
# qtl2-style genoprobs object containing one list element per chromosome. Each
# list element is a 3-dimensional array with samples in rows, founders in 
# columns, and markers in slices. The list is named by chromosome.
combine_genoprobs = function(file_list) {

  # Find the chromosome names in the filenames.
  chr_pos = regexpr(pattern = 'chr', text = file_list)

  chr = substr(file_list, chr_pos, chr_pos + 4)
  chr = gsub('^chr|\\.$', '', chr)
  # This generates a warning for Chr X. OK to ignore.
  chr = chr[order(as.numeric(chr))]
   
  probs = setNames(vector('list', length(chr)), nm = chr)
  
  for(curr_chr in chr) {
  
    # Find the file which matches the current chromsome.
    # Need to add some text after the chromosome so that 'chr1' doesn't
    # return chr1 & chr 10.
    f = file_list[grep(paste0('chr', curr_chr, '.rds'), file_list)]
  
    # Read probs object.
    p = readRDS(f)
  
    # Assign chromosome probs to probs object.
    probs[[curr_chr]] = p[[1]]
  
  } # for(curr_chr)

  # Set the attributes for the probs object.
  attr(probs, 'crosstype')    = 'do'
  attr(probs, 'is_x_chr')     = setNames(c(rep(FALSE, 19), TRUE), names(probs))
  attr(probs, 'alleles')      = LETTERS[1:8]
  attr(probs, 'alleleprobs')  = grepl('allele', file_list[1])
  class(probs) = c('calc_genoprob', 'list')

  return(probs)

} # combine_genoprobs()

