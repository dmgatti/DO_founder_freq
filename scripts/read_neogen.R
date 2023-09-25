require(data.table)
require(tidyverse)

# Read in the Neogen files and return genotypes and intensities.
# Arguments:
# neogen_dir: character string containing the full path to the top-level
#             project directory.
# markers: data.frame of markers from U Wisc. 
# Returns:
# list containing two elements: "geno" which contains two-letter genotypes,
# and "inten" which contains contains mean X & Y intensities.
read_neogen = function(neogen_dir, markers) {

  # Get the directories.
  neogen_dirs = list.dirs(path = neogen_dir, full.names = TRUE, recursive = FALSE)

  # Read in the FinalReport files.
  geno  = NULL
  inten = NULL

  for(i in 1:length(neogen_dirs)) {

    print(paste('Reading', neogen_dirs[i]))

    fr_file = dir(neogen_dirs[i], pattern = 'FinalReport.zip',
                  full.names = TRUE)

    if(length(fr_file) == 1) {

      fr = fread(cmd = paste0('gunzip -cq ', fr_file), skip = 9,
                 sep = '\t', header = TRUE, colClasses = rep(c('character', 'numeric'), c(8, 3))) %>%
             select(`SNP Name`:`Allele2 - Forward`, X:Y) %>%
             rename(marker = `SNP Name`, 
                    id     = `Sample ID`) %>%
             mutate(id = str_c(basename(neogen_dirs[i]), id, sep = ';')) %>%
             unite('geno', `Allele1 - Forward`, `Allele2 - Forward`, sep = '')

      # Get genotypes.
      g = fr %>%
            select(marker, id, geno) %>%
            pivot_wider(names_from = id, values_from = geno)

      if(is.null(geno)) {
        geno = g
      } else {
        geno = full_join(geno, g, by = 'marker')
      } # else

      # Get mean X & Y chromosome intensities.
      ii = fr %>%
             select(marker, id, X, Y) %>%
             left_join(select(markers, marker, chr), by = 'marker') %>%
             filter(chr %in% c('X', 'Y')) %>%
             mutate(chr = str_c('chr', chr)) %>%
             pivot_longer(cols = X:Y, names_to = 'channel', values_to = 'intensity') %>%
             group_by(chr, id) %>%
             summarize(mean_x = mean(intensity, na.rm = TRUE),
                       .groups = 'drop') %>%
             pivot_wider(names_from = chr, values_from = mean_x) %>%
             mutate(id = as.character(id))

        rm(g, fr)
        gc()

        inten = bind_rows(inten, ii)

    } else {

       print(paste('Warning: No FinalReport in', neogen_dirs[i]))
    
    } # else
  } # for(i)
  
  return(list(geno = geno, inten = inten))

} # read_neogen()

