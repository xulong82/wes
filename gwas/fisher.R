library(dplyr)

rm(list = ls())

arg <- commandArgs(TRUE) 
chr <- paste("chr", arg, sep = "")

load("~/Dropbox/GitHub/wes/data/pheno.rdt")

load(paste0("/data/xwang/wes/golden/geno_", chr, ".rdt"))

pheno <- pheno[match(colnames(geno), pheno$SRR), ]

Ad <- rep(0, nrow(pheno))
Ad[pheno$status == "case"] = 1

chisq = apply(geno, 1, function(x) { tbl = table(Ad, x); chisq.test(tbl)$p.value })
# fisher = apply(geno, 1, function(x) { tbl = table(Ad, x); fisher.test(tbl)$p.value })

save(chisq, file = paste0("/data/xwang/wes/fisher/", chr, ".rdt"))

