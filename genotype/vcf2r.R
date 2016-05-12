library(VariantAnnotation)

setwd("/data/xwang/wes2/merged2")

for (chr in 1:22) {
  cat(paste0("chr", chr), "\n")
  fl = paste0("chr", chr, ".merged.filtered.recode.vcf")
  geno = readGT(fl)
  geno[ geno == "0/0" | geno == "." ] = 0
  geno[ geno == "0/1" ] = 1
  geno[ geno == "1/1" ] = 2
  class(geno) = "numeric"
  save(geno, file = paste0("../R/chr", chr, ".rdt"))
}

# test

fl <- "/data/xwang/wes2/test.vcf"
geno = readGT(fl)

geno[ geno == "0/0" | geno == "." ] = 0
geno[ geno == "0/1" ] = 1
geno[ geno == "1/1" ] = 2

class(geno) = "numeric"

