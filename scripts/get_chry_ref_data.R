################################################################################
# Gather the Chr Y refrence data for the DO founders.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2024-02-07
################################################################################

library(qtl2convert)

##### VARIABLES #####

base_dir = '/compsci/gedi'

ref_dir = file.path(base_dir, 'reference')

src_dir = '/gedi/resource/MUGA/raw_neogen_data/GigaMUGA/Controls'

geno_file = file.path(src_dir, 'gigamuga_geno.txt')

marker_file = file.path(ref_dir, 'gm_uwisc_v4.csv')


##### MAIN #####

# Read in the markers
markers = read.csv(marker_file)

# Subset to retain Chr Y markers.
markers = subset(markers, chr == 'Y')

# Read in the genotypes.
geno = read.delim(geno_file)

# Get the founders.
geno = geno[,c(1, grep('^(AA|BB|CC|DD|EE|FF|GG|HH)_', colnames(geno)))]
colnames(geno)[1] = 'marker'
rownames(geno) = geno$marker
geno = as.matrix(geno[,-1])

# Subset to retain markers on Chr Y.
geno = geno[markers$marker,]

# Keep only males.
cn    = colnames(geno)
cn    = sub('_[0-9]+_R[0-9]+C[0-9]+$', '', cn)
males = grep('M$', cn)
geno  = geno[,males]

# Get the strain groups.
grps = factor(substring(colnames(geno), 1, 2))

cons_geno = sapply(levels(grps), function(z) { find_consensus_geno(genotypes = geno[,grps == z]) })

# Go back and do BL6 manually.
bb    = geno[,grep('^BB', colnames(geno))]
bb_gt = apply(bb, 1, function(z) { unique(z[z %in% c('AA', 'CC', 'GG', 'TT')]) })
bb_gt[sapply(bb_gt, length) == 0] = '--'

cons_geno[,'BB'] = unlist(bb_gt)
cons_geno[is.na(cons_geno)] = '--'

# Reformat and write out the Chr Y genotypes.
cons_geno           = data.frame(marker = rownames(cons_geno), cons_geno)
colnames(cons_geno) = c('marker', LETTERS[1:8])

# PWK has a "TG" at marker JAX00725095. This may mean that there is some
# homology with a region on Chr X. But I don't want to meet up the calling
# algorithm, so I'm setting it to "--".
cons_geno['JAX00725095', 'G'] = '--'

write.csv(cons_geno, file = file.path(ref_dir, 'gm_founder_geno_chry.csv'), 
          quote = FALSE, row.names = FALSE)


