qc <- list()

name <- list.files(path = "/data/xwang/wes/qc", pattern = "*.vcf")
 
for (i in 1:length(name)) { if (i %% 1e2 == 0) cat(i, "\n")
  filepath <- file.path("/data/xwang/wes/qc", name[i])
  qc1 <- read.delim(filepath, sep = "", stringsAsFactors = F)
 
  rownames(qc1) <- NULL
  qc1$UID <- paste(qc1$CHR, qc1$POS, sep = "-")
  qc[[i]] <- qc1[c("UID", "QUAL")]
}
 
names(qc) <- gsub(".vcf", "", name)

save(qc, file = "/data/xwang/wes/qcAll.rdt")

chr = c(2:5, 7:10, 12:18, 20:22)

cd /data/xwang/wes/golden 

metaAll = lapply(chr, function(x) {
  load(paste0("meta_", x, ".rdt")); meta
})

metaAll = do.call(rbind, metaAll)

table(metaAll$PR > 0.8)




   
