library(rstan)
library(parallel)
library(dplyr)

#--- Setup

rm(list = ls())

load("~/Dropbox/GitHub/wes/mdata.rdt")
load("/data/xwang/wes2/plink/kinship.rdt")

cov <- mdata[c("Age", "Sex")]
dat <- list(N = 9949, D = 2, Ad = mdata$AD, cov = cov, L = L)

glmm <- stan_model("~/Dropbox/GitHub/wes/glmm/glmm2.stan")
dat0 <- within(dat, {g = rep(0, 9949)})


sam0 = sampling(glmm, data = dat0, pars = c("p", "beta"), chain = 1, iter = 400, warmup = 200)

save(sam0, file = "/data/xwang/wes2/gwas/sam0.rdt")
# summary(fit, pars = c("a", "p", "beta"))$summary

#--- GWAS

# chr <- commandArgs(TRUE) 
# load(paste0("/data/xwang/wes2/geno/chr_", chr, ".rdt"))
# 
# fit0 = optimizing(glmm, data = dat0)
# init = list(p = 0, beta = fit$par[paste0("beta[", 1:2, "]")], c = fit$par[paste0("c[", 1:3, "]")])
# init$Sigma = fit$par["Sigma"]
# init$z = fit$par[paste0("z[", 1:576, "]")]
# save(init, file = "~/Dropbox/GitHub/glmm/data/init_glmm.rdt")
# load("~/Dropbox/GitHub/glmm/data/init_glmm.rdt")
# 
# n.core <- 20 
# u.core <- round(nrow(geno) / n.core)
# idx1 <- ((1:n.core) - 1) * u.core + 1 
# idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
# genoList <- mclapply(1:n.core, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)
# 
# y <- mclapply(genoList, myGWAS_optimizing, mc.cores = n.core)
# fit <- do.call(rbind, y)
# 
# save(fit, file = paste0("/data/xwang/wes2/glmm/chr_", chr, ".rdt"))

