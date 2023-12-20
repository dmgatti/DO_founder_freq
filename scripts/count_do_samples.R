################################################################################
# Count the number of DO samples in the MUGA, MegaMUGA, & GigaMUGA directories.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2023-11-17
################################################################################

library(tidyverse)

##### VARIABLES #####

# Base directory for Gedi data.
base_dir = '/gedi/resource/MUGA/raw_neogen_data'

# MUGA directory.
muga_dir = file.path(base_dir, 'MUGA')

# MUGA inventory file.
muga_inventory_file = file.path(muga_dir, 'MUGA_Samples_20170905.csv')

# MegaMUGA directory.
mm_dir = file.path(base_dir, 'MegaMUGA')

# MegaMUGA inventory file.
mm_inventory_file = file.path(mm_dir, 'MegaMUGA.samples.20160412.csv')

# GigaMUGA directory.
gm_dir = file.path(base_dir, 'GigaMUGA')

# GigaMUGA inventory file.
gm_inventory_file = file.path(gm_dir, 'GigaMUGA.samples.20231106.csv')

# Output directory.
out_dir = '/compsci/gedi/DO_founder_freq/data/sample_inventory'

##### MAIN #####

# MUGA

# Read in the inventory file.
muga_inv = read_csv(muga_inventory_file)

# Retain samples with generation data, which should be DO.
muga_inv = subset(muga_inv, !is.na(muga_inv$`DO Generation`))
muga_inv = subset(muga_inv, !muga_inv$`DO Generation` %in% c('f1' ,'founder'))

# Write out generation counts.
count(muga_inv, `DO Generation`) %>%
  write_csv(file = file.path(out_dir, 'muga_sample_table.csv'))


# MegaMUGA

# Read in the inventory file.
mm_inv = read_csv(mm_inventory_file)

# Retain samples with generation data, which should be DO.
mm_inv = subset(mm_inv, !is.na(mm_inv$`DO Generation`))
mm_inv = subset(mm_inv, !is.na(as.numeric(mm_inv$`DO Generation`)))

# Write out generation counts.
count(mm_inv, `DO Generation`) %>%
  write_csv(file = file.path(out_dir, 'megamuga_sample_table.csv'))


# GigaMUGA

# Read in the inventory file.
gm_inv = read_csv(gm_inventory_file)

# Retain samples with generation data, which should be DO.
gm_inv = subset(gm_inv, !is.na(gm_inv$`DO Generation`))
gm_inv = subset(gm_inv, !is.na(as.numeric(gm_inv$`DO Generation`)))

# Write out generation counts.
count(gm_inv, `DO Generation`) %>%
  write_csv(file = file.path(out_dir, 'gigamuga_sample_table.csv'))




