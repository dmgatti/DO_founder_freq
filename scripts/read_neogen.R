################################################################################
# Read a set of Neogen FinalReport files for one project.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-26
################################################################################

suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(tidyverse))

# Read all of the FinalReport files for one project.  Get the subdirectories
# and read in each FinalReport file. Separate into genotypes and intensities
# and return these in a list.
#
# Arguments:
# top_dir: path to the top-level project directory. The Neogen directories
#          containing the FinalReport files should be in this directory.
# Returns:
# list containing two elements: 
#    geno: data.frame containing genotypes. Markers in rows, samples in columns.
#          'marker' in first column.
#    inten_x: data.frame containing X channel intensities.
#    inten_y: data.frame containing Y channel intensities. 
read_neogen = function(top_dir) {

  # Get the neogen directories.
  neogen_dirs = list.dirs(path = top_dir, full.names = TRUE, recursive = FALSE)

  # Read in the FinalReport files.
  geno    = NULL
  inten_x = NULL
  inten_y = NULL

  for(i in 1:length(neogen_dirs)) {

    print(neogen_dirs[i])

    # Look for FinalReport file. We want it zipped.
    fr_file = dir(neogen_dirs[i], pattern = 'FinalReport.zip', full.names = TRUE)

    # If FinalReport file was not found, abort.
    if(length(fr_file) == 0) {
    
      print(paste('No FinalReport file in', neogen_dirs[i]))
      quit(1)    
    
    } # if(length(fr_file) == 0)

    #####
    # Read in FinalReport file. Select marker, id, genotype and intensities columns.
    # Append Neogen directory to sample ID.
    fr = fread(cmd = paste0('gunzip -cq ', fr_file), skip = 9,
               sep = '\t', header = TRUE, colClasses = rep(c('character', 'numeric'), c(8, 3))) %>%
           select(`SNP Name`:`SNP Name`:`Allele2 - Forward`, X:Y) %>%
           rename(marker = `SNP Name`, id = `Sample ID`) %>%
           unite('geno', `Allele1 - Forward`, `Allele2 - Forward`, sep = '') %>%
           mutate(id = str_c(basename(neogen_dirs[i]), id, sep = ';'))

    #####
    # Get genotypes as two-letters.
    g = fr %>%
          select(marker, id, geno) %>%
          pivot_wider(names_from = id, values_from = geno)

    if(is.null(geno)) {
    
      geno = g

    } else {

      geno = full_join(geno, g, by = 'marker')

    } # else

    #####
    # Get intensities on X and Y channels.
    ii = fr %>%
           select(marker, id, X, Y) 
           
    ix = ii %>%
           select(-Y) %>%
           pivot_wider(names_from = id, values_from = X)

    iy = ii %>%
           select(-X) %>%
           pivot_wider(names_from = id, values_from = Y)

    if(is.null(inten_x)) {

      inten_x = ix
      inten_y = iy

    } else {

      inten_x = full_join(inten_x, ix, by = 'marker')
      inten_y = full_join(inten_y, iy, by = 'marker')

    } # else

    rm(g, fr, ii, ix, iy)
    gc()

  } # for(i)

  return(list(geno = geno, inten_x = inten_x, inten_y = inten_y))

} # read_neogen()
