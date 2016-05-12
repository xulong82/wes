rm(list = ls())

load("~/Dropbox/GitHub/wes/data/pheno.rdt")
load("/data/xwang/wes/golden/geno_chr22.rdt")

pheno <- pheno[match(colnames(geno), pheno$SRR), ]

M1 <- pheno$SRR
N1 <- length(M1)

KS <- matrix(0, N1, N1, dimnames = list(M1, M1))

kin0 <- read.delim("/data/xwang/wes/kinship/autosome.kin0")
kin0 <- as.matrix(kin0)

for (j in 1:nrow(kin0)) { 
  if(j %% 1e6 == 0) cat(j, nrow(kin0), "\n") 
  index = sort(kin0[j, c("ID1", "ID2")])
  KS[index[2], index[1]] <- kin0[j, "Kinship"]
}

diag(KS) <- 1
KS[upper.tri(KS)] <- t(KS)[upper.tri(KS)]

save(KS, file = "~/Dropbox/GitHub/wes/data/kinship.rdt")

# ---

for (x in paste0("no_chr", 1:22)) { cat(x, "\n")
  kin <- read.delim(paste(x, "kin", sep = "."))
  kin0 <- read.delim(paste(x, "kin0", sep = "."))
  kin230 <- c(kin$Kinship, kin0$Kinship)
  kin23 <- cbind(kin23, kin230)
   
  KS <- KS0
  for (j in 1:nrow(idx)) KS[idx[j, 2], idx[j, 1]] <- kin230[j]
  diag(KS) <- 1
  KS[upper.tri(KS)] <- t(KS)[upper.tri(KS)]
  kinship[[x]] <- KS
}

colnames(kin23) <- c("autosome", paste0("no_chr", 1:22))
kinship$kin23 <- kin23

