#!/bin/bash

cd /data/xwang/wes/plink2
out="/home/xwang/Dropbox/GitHub/wes/plink"

~/plink_linux_x86_64/plink \
  --bfile autosome \
  --logistic --out "$out"/plink \
  --covar wes_cov.txt \
  --covar-number 1-2

~/plink_linux_x86_64/plink \
  --bfile autosome \
  --logistic --out "$out"/plink_pca1 \
  --covar wes_cov.txt \
  --covar-number 1-3

~/plink_linux_x86_64/plink \
  --bfile autosome \
  --logistic --out "$out"/plink_broad \
  --keep-fam /home/xwang/Dropbox/GitHub/wes/docs/remove.txt \
  --covar wes_cov.txt \
  --covar-number 1-2

~/plink_linux_x86_64/plink \
  --bfile autosome \
  --logistic --out "$out"/plink_remove_enriched \
  --remove-fam /home/xwang/Dropbox/GitHub/wes/docs/remove_enriched.txt \
  --covar wes_cov.txt \
  --covar-number 1-2

~/snptest_v2.5.2_linux_x86_64_static/snptest_v2.5.2 \
  -data autosome.bed wes.sample \
  -frequentist 1 \
  -pheno bin1 \
  -method newml \
  -o "$out"/snptest.out \
  -cov_names cov1 cov2

cd /data/xwang/Adsp/Plink
out="/data/xwang/Adsp"

~/plink_linux_x86_64/plink \
  --bfile autosome \
  --linear --covar adsp_cov.txt \
  --out "$out"/plink

~/snptest_v2.5.2_linux_x86_64_static/snptest_v2.5.2 \
  -data autosome.bed wgs.sample \
  -frequentist 1 \
  -pheno pheno1 \
  -method score \
  -o "$out"/snptest.out \
  -cov_names cov1 cov2

