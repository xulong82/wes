library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())

arg <- commandArgs(TRUE) 
chr <- paste("chr", arg, sep = "")

load("~/Dropbox/GitHub/wes/data/pheno.rdt")
load(paste0("/data/xwang/wes/golden/geno_", chr, ".rdt"))

model <- stan_model("~/Dropbox/GitHub/wes/glmm/noK.stan")

all(colnames(geno) == pheno$SRR)

cov <- pheno[c("Sex")]
cov$Sex <- as.integer(cov$Sex) - 1

Ad <- rep(0, nrow(pheno))
Ad[pheno$status == "case"] = 1

g0 <- rep(0, length(Ad))
dat = list(N = length(Ad), D = 1, cov = cov, g = g0, Ad = Ad)
init = list(a = 0.15, p = 0, beta = -0.06)

myGWAS_optimizing <- function(geno) {

  N1 = nrow(geno)
  vId = rownames(geno)
  pId = c("a", "p", "beta[1]")

  y1 = matrix(nrow = N1, ncol = 7, dimnames = list(vId, c(pId, paste0("se_", pId), "lp")))

  for (i in 1:N1) {
    dat1 = within(dat, {g = geno[i, ]})
    fit1 = optimizing(model, verbose = FALSE, hessian = TRUE, algorithm = "LBFGS", data = dat1)
#   fit1 = optimizing(model, verbose = FALSE, hessian = TRUE, algorithm = "LBFGS", init = init, data = dat1)
    se1 = tryCatch(sqrt(diag(solve(-fit1$hessian))), error=function(e) NULL)
    lk1 = sum(fit1$par[paste0("lp[", 1:length(Ad), "]")])
    y1[i, ] = c(fit1$par[pId], se1, lk1)
  }
  return(y1)
}

n.core <- 20 
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1 
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
genoList <- mclapply(1:n.core, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)

y <- mclapply(genoList, myGWAS_optimizing, mc.cores = n.core)
fit <- do.call(rbind, y)

save(fit, file = paste0("/data/xwang/wes/optimizing/noK/", chr, ".rdt"))

