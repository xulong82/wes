gInf <- NULL
geno <- list()

map <- data.frame(row.names = c("0/0", "0/1", "1/1", "1/2"), g = c(0, 1, 2, 3))
name <- list.files(path = "/data/xwang/wes/vcf", pattern = "*.vcf")
 
for (i in 1:length(name)) { if (i %% 1e2 == 0) cat(i, "\n")
  filepath <- file.path("/data/xwang/wes/vcf", name[i])
  geno1 <- read.delim(filepath, sep = "", stringsAsFactors = F)
 
  geno1$UID <- paste(geno1$CHR, geno1$POS, sep = "-")
  geno1$FID <- apply(geno1[c(7, 4:5)], 1, function (x) paste(x, collapse = "-"))
  geno1$GTN <- map[geno1$GT, "g"]
  rownames(geno1) <- NULL
  geno[[i]] <- geno1[c("UID", "FID", "GTN")]
 
  gInf1 <- geno1[! geno1$FID %in% gInf$FID, ]
  gInf <- rbind(gInf, gInf1[c("UID", "CHR", "POS", "ID", "REF", "ALT", "FID")])
}
 
names(geno) <- gsub(".vcf", "", name)
rownames(gInf) <- NULL

base4 <- c("A", "T", "C", "G")
myidx <- gInf$REF %in% base4 & gInf$ALT %in% base4
snpInf <- gInf[myidx, ]
mixInf <- gInf[! myidx, ]
polyInf <- mixInf[grepl("," , mixInf$ALT), ]

snpInf <- snpInf[! snpInf$UID %in% mixInf$UID, ]
pSnpId <- snpInf$UID[duplicated(snpInf$UID)]
bSnpInf <- snpInf[! snpInf$UID %in% pSnpId, ]  # bi-SNP

myidx <- polyInf$REF %in% base4 & nchar(polyInf$ALT) == 3
pSnpId <- unique(c(pSnpId, polyInf$UID[myidx]))
pSnpInf <- gInf[gInf$UID %in% pSnpId, ]  # poly-SNP

indInf <- mixInf[! mixInf$UID %in% polyInf$UID, ]
pIndId <- indInf$UID[duplicated(indInf$UID)]
bIndInf <- indInf[! indInf$UID %in% pIndId, ]  # bi-INDEL

pIndId <- unique(c(pIndId, polyInf$UID[! nchar(polyInf$ALT) == 3]))
pIndInf <- gInf[gInf$UID %in% pIndId, ]  # poly-INDEL

# length(unique(gInf$UID)) 
# length(intersect(pSnpId, pIndId))  # loci with polyINDEL and polySNP
# nrow(bSnpInf) + length(unique(pSnpInf$UID)) + nrow(bIndInf) + length(unique(pIndInf$UID))
 
gMix <- list()
gMix$bSnpInf <- bSnpInf
gMix$pSnpInf <- pSnpInf
gMix$bIndInf <- bIndInf
gMix$pIndInf <- pIndInf
   
save(geno, file = "/data/xwang/wes/genoAll.rdt")
save(gInf, file = "/data/xwang/wes/gInf.rdt")
save(gMix, file = "/data/xwang/wes/gMix.rdt")
   
sapply(gMix, nrow) # SAVED

