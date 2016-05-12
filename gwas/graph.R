rm(list = ls())
setwd("~/Dropbox/GitHub/wes")

chrlen <- read.delim("~/Dropbox/GitHub/X/genomes/human.hg19.genome", header = F)
chrlen <- chrlen[match(paste0("chr", 1:22), chrlen$V1), ]
chrlen <- cumsum(as.numeric(chrlen$V2)) * 1e-6; names(chrlen) <- c(1:22)
chrmid <- diff(c(0, chrlen)) * 0.5 + c(0, chrlen[-length(chrlen)])

options(stringsAsFactors = F)
plink <- read.table("./gwas/plink.assoc.logistic", header = T)
plink <- read.table("./gwas/plink_greg.assoc.logistic", header = T)

plink <- read.table("./gwas/plink_broad.assoc.logistic", header = T)
plink <- read.table("./gwas/plink_other.assoc.logistic", header = T)

plink <- plink[plink$TEST == "ADD", ]
plink$POS <- c(0, chrlen)[plink$CHR] + plink$BP * 1e-6

broad.na = plink$SNP[is.na(plink$P)]
other.na = plink$SNP[is.na(plink$P)]

plink = plink[! plink$SNP %in% c(broad.na, other.na), ]

vId <- plink$SNP
chr <- as.numeric(gsub("-.*", "", vId))
pos <- as.numeric(gsub(".*-", "", vId)) * 1e-6 # Mb
pos <- c(0, chrlen)[chr] + pos

manhattan <- data.frame(uid = vId, chr = chr, pos = pos, pval = -log10(plink$P))
manhattan$col <- rep("o", nrow(manhattan)) 
manhattan$col[chr %% 2 == 1] <- "e"
manhattan$peak = "N"
# manhattan$peak[manhattan$uid %in% group$SNP] = "Y"

png("./gwas/manhattan.png", width = 2e3, height = 1e3, res = 200)
png("./gwas/manhattan_broad.png", width = 2e3, height = 1e3, res = 200)
png("./gwas/manhattan_other.png", width = 2e3, height = 1e3, res = 200)
png("./gwas/manhattan_greg.png", width = 2e3, height = 1e3, res = 200)

ggplot(manhattan, aes(x = pos, y = pval, color = col, shape = peak)) + 
  geom_point(alpha = 0.9) + # geom_hline(yintercept = 7.3, color = "black") +
  scale_x_continuous(breaks = chrmid, labels = names(chrlen), limits = c(0, max(pos))) + 
  scale_color_manual(values = c("dodgerblue3", "firebrick1")) + 
  theme_bw() + xlab("") + ylab("-log10(P)") + guides(shape = F, color = F) +
  theme(legend.key = element_blank())

dev.off()

(y = manhattan[which(manhattan$pval > 45), ])
(y = manhattan[which(manhattan$pval > 4), ])

y = rownames(x)
vep_input = data.frame(chr = gsub("-.*", "", y), start = gsub("^.*-", "", y))

vep_input = data.frame(chr = y$chr, start = gsub("^.*-", "", y$uid))
vep_input$end = vep_input$start
vep_input$allele = "A/T"
vep_input$strand = "+"

write.table(vep_input, file = "./gwas/vep_input.txt", quote = F, sep = "\t", row.names = F)

y = x[order(x$P), ]

vep = read.table("./gwas/vep_greg.txt", header = T)
(genes = unique(vep$SYMBOL))

y$uid = gsub("-", "_", y$uid)
vep$SNP = gsub("_A/T", "", vep$Uploaded_variation)
(vep = vep[, c("SNP", "SYMBOL", "Consequence")])
vep = vep[! duplicated(paste0(vep$SNP, vep$SYMBOL, vep$Consequence)), ]
unique(vep$SNP)
vep$CHR = gsub("_.*", "", vep$SNP)
vep$BP = gsub("^.*_", "", vep$SNP)

vep$SNP = gsub("_", "-", vep$SNP)

greg = cbind(vep, plink[match(vep$SNP, plink$SNP), c("OR", "P")])

x = plink[plink$P < 5e-8, ]

group <- list()
group[[1]] <- x[1, ]
group_idx <- 1

for (i in 2:nrow(x)) {
  chromosome = x$CHR[i] == x$CHR[i-1]
  position = x$BP[i] - x$BP[i-1] < 1e6
  
  if ( all(chromosome, position) )
    group[[group_idx]] = rbind(group[[group_idx]], x[i, ])
  else {
    group_idx = group_idx + 1
    group[[group_idx]] = x[i, ]
  }
}

group <- lapply(group, function(x) x[which.min(x$P), ])
group <- do.call(rbind, group)


