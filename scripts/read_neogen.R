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
#    inten: data.frame containing mean Chr X & Y intensities. 
read_neogen = function(top_dir, markers) {

  # Get the neogen directories.
  neogen_dirs = list.dirs(path = top_dir, full.names = TRUE, recursive = FALSE)

  # Read in the FinalReport files.
  geno  = NULL
  inten = NULL

  for(i in 1:length(neogen_dirs)) {

    print(neogen_dirs[i])

    fr_file = dir(neogen_dirs[i], pattern = 'FinalReport', full.names = TRUE)

    if(length(fr_file) == 0) {
    
      print(paste('No FinalReport file in', neogen_dirs[i]))
      quit(1)    
    
    } # if(length(fr_file) == 0)

    fr = fread(cmd = paste0('gunzip -cq ', fr_file), skip = 9,
               sep = '\t', header = TRUE, colClasses = rep(c('character', 'numeric'), c(8, 3))) %>%
           select(`SNP Name`:`SNP Name`:`Allele2 - Forward`, X:Y) %>%
           rename(marker = `SNP Name`, id = `Sample ID`) %>%
           unite('geno', `Allele1 - Forward`, `Allele2 - Forward`, sep = '')

    g = fr %>%
          select(marker, id, geno) %>%
          pivot_wider(names_from = id, values_from = geno)

    if(is.null(geno)) {
      geno = g
    } else {
      geno = full_join(geno, g, by = 'marker')
    } # else

    ii = fr %>%
           select(marker, id, X, Y) %>%
           left_join(markers, by = 'marker') %>%
           filter(chr %in% c('X', 'Y')) %>%
           pivot_longer(cols = X:Y, names_to = 'channel', values_to = 'intensity') %>%
           group_by(chr, id) %>%
           summarize(mean_x = mean(intensity, na.rm = TRUE),
                     .groups = 'drop') %>%
           pivot_wider(names_from = chr, values_from = mean_x) %>%
           mutate(id = as.character(id))

    rm(g, fr)
    gc()

    inten = bind_rows(inten, ii)

  } # for(i)

  return(list(geno = geno, inten = inten))

} # read_neogen()
