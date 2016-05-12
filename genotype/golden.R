chr <- commandArgs(TRUE)

load("~/Dropbox/GitHub/wes/data/mdata.rdt") # pheno

load("/data/xwang/wes/genoAll.rdt") # geno
load("/data/xwang/wes/qcAll.rdt") # gMix
load("/data/xwang/wes/gMix.rdt") # qc
 
all(names(qc) %in% mdata$SRR)
all(names(geno) %in% mdata$SRR)
names(qc) <- gsub(".sorted.bam", "", names(qc))
names(geno) <- gsub(".sorted.bam", "", names(geno))

qc <- qc[mdata$SRR]
geno <- geno[mdata$SRR]
 
n.sample <- length(geno); id.sample <- names(geno)

meta <- rbind(gMix$bSnpInf, gMix$bIndInf)

meta <- meta[meta$CHR == chr, ]
meta <- meta[order(meta$POS), ]

n.var <- nrow(meta); id.var <- meta$UID
biGeno <- matrix(0, n.var, n.sample, dimnames = list(id.var, id.sample))

for (i in 1:n.sample) {
  geno1 = geno[[i]]
  geno1 = geno1[geno1$UID %in% id.var, ]
  biGeno[match(geno1$UID, id.var), i] <- geno1[, "GTN"]
}

meta$SNP <- "N"
NC <- c("A", "T", "C", "G")
meta$SNP[meta$REF %in% NC & meta$ALT %in% NC] <- "Y"

maf <- rowSums(biGeno) / n.sample / 2
meta$MA <- meta$ALT
meta$MA[maf > .5] <- meta$REF[maf > .5] 
maf[maf > .5] <- 1 - maf[maf > .5]

meta$MAF <- maf
meta$HET <- apply(biGeno, 1, function (x) sum(x == 1)) / n.sample
meta$N.SRR <- apply(biGeno, 1, function (x) sum(x > 0))

idx <- which(meta$MAF > 0.01 & meta$HET < 0.99)

id.var <- id.var[idx]; n.var = length(id.var)

meta$PR <- 2 # default

pr <- rep(0, n.var)
for(i in 1:n.var) {
  qc1 = qc[sapply(qc, function(x) id.var[i] %in% x$UID)]
  qc1 = sapply(qc1, function(x) x[match(id.var[i], x$UID),"QUAL"])
  pr[i] = mean(qc1 == "PASS")
}

meta$PR[idx] = pr
save(meta, file = paste0("/data/xwang/wes/metaAll/meta_", chr, ".rdt"))

geno <- biGeno[meta$MAF > 0.01 & meta$HET < 0.99 & meta$PR > 0.2, ]
meta <- meta[meta$MAF > 0.01 & meta$HET < 0.99 & meta$PR > 0.2, ]

save(meta, file = paste0("/data/xwang/wes/golden/meta_", chr, ".rdt"))
save(geno, file = paste0("/data/xwang/wes/golden/geno_", chr, ".rdt"))

f1 = paste0("/data/xwang/wes/golden/meta_", chr, ".txt")
f2 = paste0("/data/xwang/wes/golden/geno_", chr, ".txt")

write.table(meta, file = f1, row.names = F, col.names = F, quote = F)
write.table(geno, file = f2, row.names = F, col.names = F, quote = F)
 
