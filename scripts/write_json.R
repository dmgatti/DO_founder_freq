################################################################################
# Given the filenames, write out the JSON control file for one chromosome.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-09-23
################################################################################

# Arguments:
# out_dir: full path to the directory where the JSON will be written.
# chr: character indicating the current chromosome for these files.
# geno_file: full path to the sample genotype file for the current chromosome.
# founder_file: full path to the founder genotype file for the current chromosome.
# gmap_file: full path to the genetic map file for the current chromosome.
# pmap_file: full path to the physical map file for the current chromosome.
# pheno_file: full path to the phenotypes file for the current chromosome.
# covar_file: full path to the covariates file for the current chromosome.
write_json = function(out_dir, chr, geno_file, founder_file, gmap_file, 
                      pmap_file, pheno_file, covar_file) {

  # JSON control file

  json = paste0('{
    "description": "GigaMUGA",
    "crosstype": "do",
    "sep": ",",
    "na.strings": ["-", "NA"],
    "comment.char": "#",
    "geno": "', geno_file, '",
    "founder_geno": "', founder_file, '",
    "gmap": "', gmap_file, '",
    "pmap": "', pmap_file, '",
    "pheno": "', pheno_file, '",
    "covar": "', covar_file, '",
    "alleles": ["A", "B", "C", "D", "E", "F", "G", "H"],
    "x_chr": "X",
    "genotypes": {
      "A": 1,
      "H": 2,
      "B": 3
    },
    "geno_transposed": true,
    "founder_geno_transposed": true,
    "sex": {
      "covar": "sex",
      "female": "female",
      "male": "male"
    },
    "cross_info": {
      "covar": "gen"
    }
  }')
  
  writeLines(text = json, con = file.path(out_dir, paste0('chr', chr, '.json')))

} # write_json()


