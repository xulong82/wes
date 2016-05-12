#!/bin/bash

#------- Note: merge split tped

for idx in `seq 1 22`; do
  echo $idx
  cat chr${idx}.tped >> wes.tped
done

#------- New pipeline

hpc="/data/xwang/wes2/plink"
loc="/home/xwang/Dropbox/GitHub/wes/plink"

#------- Note: tped to bed

cd $hpc

~/plink_linux_x86_64/plink --tped wes.tped -tfam ${loc}/wes.tfam --make-bed --out wes

#------- Note: sample clustering

~/plink_linux_x86_64/plink --bfile wes --genome  --out ibs
~/plink_linux_x86_64/plink --bfile wes --read-genome ibs.genome --out mds --cluster --mds-plot 4

#------- Note: kinship by plink

~/plink_linux_x86_64/plink --bfile wes --cluster --matrix --out wes_kin

#------- Note: PCA by plink

~/plink_linux_x86_64/plink --bfile wes --pca --out ${loc}/wes_pca

#------- Note: GWAS

~/plink_linux_x86_64/plink \
  --bfile wes \
  --logistic --out ${loc}/plink \
  --covar ${loc}/wes_cov.txt \
  --covar-number 1-4

# remove "Enriched"

~/plink_linux_x86_64/plink \
  --bfile wes \
  --logistic --out ${loc}/plink_remove_enriched \
  --remove-fam ${loc}/enriched.txt \
  --covar ${loc}/wes_cov.txt \
  --covar-number 1-2

# only "Enriched"

~/plink_linux_x86_64/plink \
  --bfile wes \
  --logistic --out "$loc"/plink_enriched \
  --keep-fam ${loc}/enriched.txt \
  --covar ${loc}/wes_cov.txt \
  --covar-number 1-2

