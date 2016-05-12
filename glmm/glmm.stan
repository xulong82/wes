data {
  int<lower=1> N;  // Sample
  int<lower=1> D;  // Covariate
  row_vector[D] cov[N];  // Covariate
  cholesky_factor_cov[N] L;  // Covariance
  vector[N] g;  // genotype
  int<lower=0,upper=1> Ad[N];  // phenotype
}

parameters {
  real p;  // variant effect
  vector[D] beta;  // covariate effect
  real<lower=0.001> sigma;
  vector[N] u;
} 

model {

  vector[N] mu; 
  for (n in 1:N) mu[n] <- 0;

  p ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigma ~ inv_gamma(2, 1);
  u ~ multi_normal_cholesky(mu, sigma * L);

  for (n in 1:N)
    Ad[n] ~ bernoulli_logit(p * g[n] + cov[n] * beta + u[n]);
}

