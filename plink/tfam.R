#------- Note: make tfam file

# --- Plink TFAM FILE ---
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype (0=missing; 1=unaffected; 2=affected)

rm(list = ls())
setwd("~/Dropbox/GitHub/wes")

options(stringsAsFactors = F)
mdata <- read.delim("./docs/ADSP_wes_samplesheet_Jan2016_Updated_01262016.txt")
mdata = mdata[mdata$SampleSelection != "SeqControl", ]
table(duplicated(mdata$ADSP_SM_ID))

tfam <- read.table("./plink/wes.tfam")

all(tfam$V1 %in% mdata$ADSP_SM_ID)
mdata = mdata[match(tfam$V1, mdata$ADSP_SM_ID), ]
all(tfam$V1 == mdata$ADSP_SM_ID)

tfam$V5 = mdata$Gender + 1
tfam$V6 <- rep(1, nrow(mdata))
tfam$V6[mdata$status == "case"] = 2

write.table(tfam, file = "wes.tfam", row.names = F, col.names = F, quote = F)

#------- Note: make covariates file

cov = tfam; cov$V3 = cov$V4 = cov$V5 = cov$V6 = NULL
cov$Age = as.numeric(gsub("\\+", "", mdata$Age))
cov$Sex <- mdata$Gender

eigenvec <- read.table("./plink/wes_pca.eigenvec")[-c(1:2)]
names(eigenvec) <- paste0("PC", 1:20)

cov = cbind(cov, eigenvec[1:3])

write.table(cov, file = "./plink/wes_cov.txt", row.names = F, col.names = F, quote = F)

#------- Note: Enriched samples

enriched = mdata$ADSP_SM_ID[mdata$SampleSelection == "Enriched"]
write(enriched, ncolumns = 1, file = "./plink/enriched.txt")

#------- Note: make SNPTEST file

# Snptest sample file
# ID_1 ID_2 missing cov_1 cov_2 cov_3 cov_4 mdata1 bin1
# 0 0 0 D D C C P B
# 1 1 0.007 1 2 0.0019 -0.008 1.233 1
# 2 2 0.009 1 2 0.0022 -0.001 6.234 0
# 3 3 0.005 1 2 0.0025 0.0028 6.121 1
# 4 4 0.007 2 1 0.0017 -0.011 3.234 1
# 5 5 0.004 3 2 -0.012 0.0236 2.786 0

wes.sample <- data.frame(ID1 = 1:nrow(mdata), ID2 = mdata$SRR)
wes.sample$missing <- rep(0, nrow(mdata))
wes.sample$cov_1 <- mdata$Age
wes.sample$cov_2 <- as.numeric(mdata$Sex) - 1
wes.sample$Ad <- rep(0, nrow(mdata))
wes.sample$Ad[mdata$status == "case"] = 1

setwd("/data/xwang/wes/plink")
write.table(wes.sample, file = "wes.sample", quote = F, row.names = F, col.names = T)
