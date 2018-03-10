
## design matrix, download from https://github.com/randel/mvMISE_simulation
load('X.rdata')
source('https://raw.githubusercontent.com/randel/mvMISE_simulation/master/simulation%20functions.r')

library(parallel)
cl = makeCluster(4)
library(doParallel)
registerDoParallel(cl)
getDoParWorkers()


########################################################### mvMISE_b

## Table 1 & 3: MSE (and bias) & power

# exponential missing
res_b1 = list()
k = 0
probit = F
sigma0 = 1
for(K in c(5, 20, 100)) {
  k = k + 1
  b0 = 10
  res_b1[[k]] = sim_mvmise(K, beta = c(b0, 1.3/sqrt(K)), phi = c(ifelse(probit, .7, 0), -0.015), sigma0, sigma1 = sqrt(2), nsim = 1000, X = X,
                           probit = probit, corr = c("b"), maxIter = 1000)
}

# probit missing
res_b2 = list()
k = 0
probit = T
sigma0 = 1
for(K in c(5, 20, 100)) {
  k = k + 1
  b0 = 10
  res_b2[[k]] = sim_mvmise(K, beta = c(b0, 1.3/sqrt(K)), phi = c(ifelse(probit, .7, 0), -0.015), sigma0, sigma1 = sqrt(2), nsim = 1000, X = X,
                           probit = probit, corr = c("b"), maxIter = 1000)
}




## Table 3: type I error rate

# exponential missing, sigma0 = 1
res_b01 = list()
k = 0
probit = F
sigma0 = 1
for(K in c(5, 20, 100)) {
  k = k + 1
  b0 = 10
  res_b01[[k]] = sim_mvmise(K, beta = c(b0, 0), phi = c(ifelse(probit, .7, 0), -0.015), sigma0, sigma1 = sqrt(2), nsim = 1000, X = X,
                            probit = probit, corr = c("b"), maxIter = 1000)
}

# probit missing
res_b02 = list()
k = 0
probit = T
sigma0 = 1
for(K in c(5, 20, 100)) {
  k = k + 1
  b0 = 10
  res_b02[[k]] = sim_mvmise(K, beta = c(b0, 0), phi = c(ifelse(probit, .7, 0), -0.015), sigma0, sigma1 = sqrt(2), nsim = 1000, X = X,
                            probit = probit, corr = c("b"), maxIter = 1000)
}




################################################################ mvMISE_e

## Table 2: MSE (and bias)

res_e = list()
k = 0
for(probit in c(FALSE, TRUE)) {
  for(K in c(5, 20, 100)) {
    k = k + 1
    b0 = 10
    Sigma = sim_Sigma_power_law(K)
    res_e[[k]] = sim_mvmise(K, beta = c(b0, 1.4/sqrt(K)), phi = c(ifelse(probit, 1,0), -0.015), sigma0 = .1, sigma1 = sqrt(2), nsim = 1000, X = X, 
                            Sigma = Sigma, D = matrix(1), probit = probit, corr = c("e"), maxIter = 1000, lambda = sqrt(log(K)/(36*4)))
  }
}



## Table 4: type I error rate based on permutation

res_e0 = list()
k = 0
for(probit in c(FALSE, TRUE)) {
  for(K in c(5, 20, 100)) {
    k = k + 1
    b0 = 10
    Sigma = sim_Sigma_power_law(K)
    res_e0[[k]] = sim_mvmise_e_fisher(K, beta = c(b0, 0), phi = c(ifelse(probit, 1,0), -0.015), sigma0 = .1, sigma1 = sqrt(2), nsim = 1000, X = X, 
                                      Sigma = Sigma, probit = probit, maxIter = 1000, lambda = sqrt(log(K)/(36*4)), beta_prop = 1, beta_norm = F)
  }
}



## Table 4: power based on permutation

# all equal signals, correlated errors
res_e_power = list()
k = 0
for(probit in c(FALSE, TRUE)) {
  for(K in c(5, 20, 100)) {
    k = k + 1
    b0 = 10
    Sigma = sim_Sigma_power_law(K)
    res_e_power[[k]] = sim_mvmise_e_fisher(K, beta = c(b0, .6/log10(K)), phi = c(ifelse(probit, 1,0), -0.015), sigma0 = .1, sigma1 = sqrt(2), nsim = 1000, X = X,
                                           Sigma = Sigma, probit = probit, maxIter = 1000, lambda = sqrt(log(K)/(36*4)), beta_prop = 1, beta_norm = F)
  }
}


# all equal signals, indpendent errors
res_e_power_ind = list()
k = 0
for(probit in c(FALSE, TRUE)) {
  for(K in c(5, 20, 100)) {
    k = k + 1
    b0 = 10
    Sigma = diag(K)
    res_e_power_ind[[k]] = sim_mvmise_e_fisher(K, beta = c(b0, .6/log10(K)), phi = c(ifelse(probit, 1,0), -0.015), sigma0 = .1, sigma1 = sqrt(2), nsim = 1000, X = X,
                                               Sigma = Sigma, probit = probit, maxIter = 1000, lambda = sqrt(log(K)/(36*4)), beta_prop = 1, beta_norm = F)
  }
}


# half signals, correlated errors
res_e_power_prop = list()
k = 0
for(probit in c(FALSE, TRUE)) {
  for(K in c(5, 20, 100)) {
    k = k + 1
    b0 = 10
    Sigma = sim_Sigma_power_law(K)
    res_e_power_prop[[k]] = sim_mvmise_e_fisher(K, beta = c(b0, .8/log10(K)), phi = c(ifelse(probit, 1,0), -0.015), sigma0 = .1, sigma1 = sqrt(2), nsim = 1000, X = X,
                                                Sigma = Sigma, probit = probit, maxIter = 1000, lambda = sqrt(log(K)/(36*4)), beta_prop = .5, beta_norm = F)
  }
}


# random signals, correlated errors
res_e_power_norm = list()
k = 0
for(probit in c(FALSE, TRUE)) {
  for(K in c(5, 20, 100)) {
    k = k + 1
    b0 = 10
    Sigma = sim_Sigma_power_law(K)
    res_e_power_norm[[k]] = sim_mvmise_e_fisher(K, beta = c(b0, 1/log10(K)), phi = c(ifelse(probit, 1,0), -0.015), sigma0 = .1, sigma1 = sqrt(2), nsim = 1000, X = X, 
                                                Sigma = Sigma, probit = probit, maxIter = 1000, lambda = sqrt(log(K)/(36*4)), beta_prop = 1, beta_norm = T)
  }
}
