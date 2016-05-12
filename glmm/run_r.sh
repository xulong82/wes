#!/bin/bash
#PBS -l mem=64gb,nodes=1:ppn=1,walltime=12:00:00

module load R

# w/o parameter
# Rscript --no-restore --quiet ~/Dropbox/GitHub/wes/glmm/glmm.R
Rscript --no-restore --quiet ~/Dropbox/GitHub/wes/glmm/gmmat.R

# w/ parameter
# Rscript --no-restore --quiet ~/Dropbox/GitHub/wes/optimizing.R ${index}

