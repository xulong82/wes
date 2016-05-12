#!/bin/bash
#PBS -l mem=64gb,nodes=1:ppn=1,walltime=30:00:00

module load R
# module load R/3.0.3
# module load R/3.1.1

# w/o parameter
# Rscript --no-restore --quiet ~/Dropbox/GitHub/wes/genotypes/geno.R
# Rscript --no-restore --quiet ~/Dropbox/GitHub/wes/genotypes/qc.R

# w/ parameter
  Rscript --no-restore --quiet ~/Dropbox/GitHub/wes/genotypes/golden.R ${index}
# Rscript --no-restore --quiet ~/Dropbox/GitHub/wes/optimizing.R ${index}
# Rscript --no-restore --quiet ~/Dropbox/GitHub/wes/fisher.R ${index}

