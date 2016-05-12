# KINSHIP IN R

library(matrixcalc)

setwd("/data/xwang/wes2/plink")

date()
kin = read.table("wes_kin.mibs")
date()

id = read.table("wes_kin.mibs.id")

rownames(kin) = colnames(kin) = id$V1

is.positive.definite(as.matrix(kin))

save(kin, file = "/data/xwang/wes2/plink/kin.rdt")

L = t(chol(kin))
save(L, file = "/data/xwang/wes2/plink/chol_kin.rdt")

# VEP INPUT

vep.input = read.table("/data/xwang/wes2/plink/wes.bim")
vep.input$V7 = paste(vep.input$V5, vep.input$V6, sep = "/")
vep.input = vep.input[c("V1", "V4", "V4", "V7")]
vep.input$V8 = "+"

write.table(vep.input, file = "/data/xwang/wes2/vep_input.txt", quote = F, sep = "\t", row.names = F)

vep = read.table("/data/xwang/wes2/vep_output.txt", stringsAsFactors = F, header = T)

