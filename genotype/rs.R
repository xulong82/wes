chr <- commandArgs(TRUE)  # 1-22, X, Y

load("~/Dropbox/GitHub/wes/data/pheno.rdt") # pheno
load(paste0("/data/xwang/wes/geno/geno_chr", chr, ".rdt"))  # geno
load(paste0("/data/xwang/wes/geno/gInf_chr", chr, ".rdt"))  # gMix

geno <- geno[names(geno) %in% pheno$SRR] 
n.sample <- length(geno); id.sample <- names(geno)

# rs variants only

rsInf <- gInf[grepl("rs", gInf$ID), ]
rsInf <- rsInf[order(rsInf$POS), ]

duprs <- rsInf$UID[duplicated(rsInf$UID)]
rsInf <- rsInf[! rsInf$UID %in% duprs, ]

n.var <- nrow(rsInf); id.var <- rsInf$UID

rs <- matrix(0, n.var, n.sample, dimnames = list(id.var, id.sample))

for (i in 1:n.sample) {
  geno1 <- geno[[i]]
  geno1 <- geno1[geno1$UID %in% id.var, ]
  rs[match(geno1$UID, id.var), i] <- geno1[, "GTN"]
}

rsInf$SNP <- "N"
NC <- c("A", "T", "C", "G")
rsInf$SNP[rsInf$REF %in% NC & rsInf$ALT %in% NC] <- "Y"

maf <- rowSums(rs) / n.sample / 2
rsInf$MA <- rsInf$ALT
rsInf$MA[maf > .5] <- rsInf$REF[maf > .5] 
maf[maf > .5] <- 1 - maf[maf > .5]

# Note: ALT allele, not necessarily the minor alle, was tested

rsInf$MAF <- maf
rsInf$HET <- apply(rs, 1, function (x) sum(x == 1)) / n.sample
rsInf$N.SRR <- apply(rs, 1, function (x) sum(x > 0))

index <- which(rsInf$MAF > 0.01 & rsInf$HET < 0.99)

geno <- rs[index, ]
meta <- rsInf[index, ]

save(meta, file = paste0("/data/xwang/wes/rs/meta_chr", chr, ".rdt"))
save(geno, file = paste0("/data/xwang/wes/rs/geno_chr", chr, ".rdt"))

