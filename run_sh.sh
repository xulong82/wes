#!/bin/bash
#PBS -l mem=64gb,nodes=1:ppn=1,walltime=12:00:00

# sleep 60m

# ~/Dropbox/GitHub/wes/genotypes/myVcf.sh

# w/ parameter
~/Dropbox/GitHub/wes/genotypes/plink.sh ${index}
# ~/Dropbox/GitHub/Adsp/geno/kinship.sh ${index}

