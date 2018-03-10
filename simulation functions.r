
#################### simulate data
# nsim: number of replications
# X_real: design matrix X
# K: number of multivariate outcomes
# beta: true fixed effects
# phi: parameters in the missing-data model, covariates in the order: 1, y, x
# probit: probit missing-data mechanism or exponential if FALSE 
# error sd: sigma0, sigma1
# cov of RE: D
# Sigma: covariance matrix for the errors
# D_exp: var of experiment-level RE
# irt: IRT variance component for RE

# beta_prop: proportion of genes with nonzero beta of interest when correlatedError = TRUE
# beta_norm: simulate betas from a uniform distribution

sim_dat = function(nsim = 1000, X_real, K = 4, beta = c(10, -0.7, 0.7), 
                   phi = c(0, -0.025, 0), sigma0 = 1, sigma1 = sqrt(2), D, D_exp = 0, 
                   Sigma, probit = F, remove_allMiss = T, irt = F, correlatedError = F, 
                   seed = 1, beta_prop = 1, beta_norm = F) {
  set.seed(seed)
  nexp = nrow(X_real)
  ni = dim(X_real)[3]  # samples within cluster/experiment
  
  ## multivariate
  library(mvtnorm)
  
  dat_all = list()
  mrate = rep(NA, nsim)
  
  for (k in 1:nsim) {
    
    id = rep(1:(nexp), rep(ni * K, nexp))
    Y = Y_full = array(NA, dim = c(nexp, ni, K))
    y = X = NULL
    
    ## protein-specific beta
    beta_vec = rbinom(K, 1, beta_prop) * beta[length(beta)]
    # negative + positive effects
    if (beta_norm) 
      if (beta[length(beta)] != 0) 
        beta_vec = runif(K, -beta[length(beta)], beta[length(beta)])
    
    for (i in 1:nexp) {
      Yi = matrix(NA, ni, K)
      if (irt) 
        b = rnorm(1) * D else b = rmvnorm(n = 1, sigma = D)
        
        # ni x k
        Xi = t(X_real[i, , ])
        
        b_exp = rnorm(n = 1, sd = sqrt(D_exp))
        
        if (correlatedError) {
          error = matrix(rmvnorm(n = 1, sigma = kronecker(Sigma, diag(c(sigma0, rep(sigma1, ni - 1))^2))), nrow = ni)
          b = rep(rnorm(n = 1, sd = sqrt(D)), K)  # single random intercept
        }
        
        # selection model
        for (j in 1:K) {
          
          if (correlatedError) {
            beta[length(beta)] = beta_vec[j]
            Yi[, j] = Xi %*% beta + b[j] + error[, j] + b_exp
          } else {
            error = c(rnorm(n = 1, sd = sigma0), rnorm(n = ni - 1, sd = sigma1))
            Yi[, j] = Xi %*% beta + b[j] + error + b_exp
          }
          
          # check sensitivity of missing-data model: probit
          p = ifelse(probit, pnorm(sum(phi * c(1, sum(Yi[, j]), sum(Xi[, ncol(Xi)])))), min(exp(sum(phi * c(1, sum(Yi[, j]), sum(Xi[, ncol(Xi)])))), 1))  
          # peptide missing probability: exp()
          
          Y_full[i, , j] = Yi[, j]
          if (rbinom(n = 1, size = 1, p = p) == 1) 
            Yi[, j] = NA
          Y[i, , j] = Yi[, j]
        }
        
        y = c(y, as.vector(Yi))
        X = rbind(X, Xi[rep(1:ni, K), ])
        
    }
    
    # option to remove experiments with no observations
    if (remove_allMiss) {
      id0 = NULL
      for (j in 1:(nexp)) {
        if (mean(is.na(y[id == j])) == 1) 
          id0 = c(id0, j)
      }
      if (length(id0) > 0) {
        idx = which(id %in% id0)
        id = rep(1:(nexp - length(id0)), rep(ni * K, nexp - length(id0)))
        y = y[-idx]
        X = X[-idx, ]
      }
    }
    
    mrate[k] = mean(is.na(y))
    
    dat = list(y = y, X = X, id = id, K = K, Y = Y, Y_full = Y_full)
    
    dat_all[[k]] = dat
  }
  
  # print(summary(mrate))
  return(list(dat = dat_all, mrate = mean(mrate)))
}


#################### analyze simulated data
### parameters
# K: the number of outcomes
# beta: fixed effects, first element is the intercept, the second is the effect of interest 
# sigma0, sigma1: standard deviations of errors
# icc: determines tau (sd of RE) based on sigma0 if tau == NULL
# lambda: tuning parameter for ADMM algorithm
# corr: 'b' if random effects (b) are correlated, 'e' if error terms are correlated

sim_mvmise = function(K, beta = c(10, 0.7), phi = c(0, -0.025), sigma0 = 1, sigma1 = sqrt(2), icc = 0.7, nsim = 1000, X, seed = 1, tau = NULL, 
                      permutation = F, probit = F, Sigma = NULL, D = NULL, lambda = 0.5, corr = c("b", "e"), maxIter = 1000) {
  
  # simulate data
  Z = kronecker(diag(K), rep(1, 4))
  
  if(length(phi) == 2) phi0 = c(phi, 0) else phi0 = phi
  set.seed(seed)  # only work for mvmise_b
  
  if (corr == "b") {
    ratio = (icc/(1 - icc))  # 7/3
    if (is.null(tau)) 
      tau = sqrt(runif(K, max(sigma0^2 - 1, 0), sigma0^2 + 1) * ratio)
    
    dat = sim_dat(nsim = nsim, X_real = X, K, beta = beta, phi = phi0, 
                  sigma0 = sigma0, sigma1 = sigma1, D = tau, D_exp = 0, Sigma, 
                  probit = probit, remove_allMiss = F, irt = T, correlatedError = F)
    # average multivariate data to univariate data
    dat1 = array(NA, dim = c(dim(dat$dat[[1]]$Y)[-3], nsim)) # nexp, ni, nsim
    for (i in 1:nsim) dat1[, , i] = apply(dat$dat[[i]]$Y, 1:2, mean, na.rm = T)
  }
  
  if (corr == "e") {
    dat0 <- foreach(i = 1:nsim, .errorhandling = "pass") %dopar% {
      source('https://raw.githubusercontent.com/randel/mvMISE_simulation/master/simulation%20functions.r')
      sim_dat(nsim = 1, X_real = X, K = K, beta = beta, phi = phi0, 
              sigma0 = sigma0, sigma1 = sigma1, D = D, D_exp = 0, Sigma = Sigma, 
              probit = probit, remove_allMiss = F, correlatedError = T, seed = i)
    }
    
    (err <- sapply(dat0, function(x) is.null(x$dat[[1]]$Y)))
    nsim = sum(!err)
    if(sum(err) > 0) print(c(nsim, dat0[err][[1]]))
    
    dat = array(NA, dim = c(dim(dat0[!err][[1]]$dat[[1]]$Y), nsim)) # nexp, ni, K, nsim
    for (i in 1:nsim) dat[, , , i] = dat0[!err][[i]]$dat[[1]]$Y
    
    dat1 = array(NA, dim = c(dim(dat[, , , 1])[-3], nsim)) # nexp, ni, nsim
    for (i in 1:nsim) dat1[, , i] = apply(dat[, , , i], 1:2, mean, na.rm = T)
  }
  
  # print('multi/uni-variate missing rates')
  if (corr == "b") print(c(K = K, sigma0 = sigma0, probit = probit, b0 = b0, mrate = dat$mrate, mrate1 = mean(is.na(dat1))))
  if (corr == "e") print(c(K = K, sigma0 = sigma0, probit = probit, b0 = b0, mrate = mean(is.na(dat)), mrate1 = mean(is.na(dat1))))
  
  
  #################################################################################### mvMISE
  
  X0 = apply(X, 2, as.vector)
  
  # Y an outcome matrix. Each row is a sample, and each column is an outcome variable; 
  # X a covariate matrix. Each row is a sample, and each column is a covariate
  
  if (corr == "b") {
    t_mvmise = system.time(res_mvmise <- foreach(i = 1:nsim, .errorhandling = "pass") %dopar% 
    {
      library(lme4)
      library(MASS)
      Y = apply(dat$dat[[i]]$Y, 3, as.vector)
      
      source("https://raw.githubusercontent.com/randel/mvMISE/master/R/mvmise_b.r")
      mvMISE_b(Y, X0, id = rep(1:36, 4), sigma_diff = T, maxIter = maxIter)
    })
  }
  
  if (corr == "e") {
    t_mvmise = res_mvmise <- foreach(i = 1:nsim, .errorhandling = "pass") %dopar% 
    {
      library(lme4)
      library(MASS)
      Y = apply(dat[, , , i], 3, as.vector)
      source("https://raw.githubusercontent.com/randel/mvMISE/master/R/mvmise_e.r")
      mvMISE_e(Y, X0, id = rep(1:36, 4), sigma_diff = T, lambda = lambda, maxIter = maxIter)
    }
  }
  
  (err <- sapply(res_mvmise, function(x) length(x) == 2))
  
  if (sum(err) > 0) {
    print("mvMISE algorithm error #")
    print(sum(err))
    print(res_mvmise[[which(err)[1]]])
  }
  
  ######################################################################################## mixEMM 
  # estimated with mvMISE (either mvMISE_b or mvMISE_e) with K=1
  
  t_mixemm = system.time(res_mixemm_ave <- foreach(i = 1:nsim, .errorhandling = "pass") %dopar% 
  {
    library(lme4)
    library(MASS)
    Y = apply(dat1[, , i, drop = F], 3, as.vector)
    
    source("https://raw.githubusercontent.com/randel/mvMISE/master/R/mvmise_b.r")
    mvMISE_b(Y, X0, id = rep(1:36, 4), sigma_diff = T, maxIter = maxIter)
  })
  
  table(err <- sapply(res_mixemm_ave, function(x) length(x) == 2))
  if (sum(err) > 0) {
    print("mixEMM algorithm error #")
    print(sum(err))
  }
  

  
  ######################################################################################## lm
  
  t_lm = system.time(res_lm_ave <- foreach(i = 1:nsim, .errorhandling = "pass") %dopar% 
  {
    y0 = as.vector(apply(dat1[, , i], 1, function(x) x[-1] - x[1]))  # log-transformed / normal data
    X0 = NULL
    for (j in 1:nrow(dat1[, , i])) X0 = rbind(X0, t(X[j, , -1]))
    lm(y0 ~ X0 - 1)
  })
  
  #################################################### print all results
  
  res_list = list(res_mvmise = res_mvmise, res_mixemm_ave = res_mixemm_ave, res_lm_ave = res_lm_ave)
  res_summary = summary_res(res_list, beta, corr)
  
  
  ################################################################### permute each sample for mvMISE
  # only for testing
  
  if (permutation) {
    
    stat_mvMISE0 <- foreach(i = 1:nsim, .errorhandling = "pass") %dopar% 
    {
      library(lme4)
      library(MASS)
      
      Y = apply(dat[sample(1:nrow(dat)), , , i], 3, as.vector)
      source("https://raw.githubusercontent.com/randel/mvMISE/master/R/mvmise_e.r")
      res = mvMISE_e(Y, X0, id = rep(1:36, 4), sigma_diff = T, lambda = lambda, maxIter = maxIter)
      res$beta[2]/res$se[2]
    }
    
    stat_mixemm0 <- foreach(i = 1:nsim, .errorhandling = "pass") %dopar% 
    {
      library(lme4)
      library(MASS)
      Y = apply(dat1[sample(1:nrow(X)), , i, drop = F], 3, as.vector)
      
      source("https://raw.githubusercontent.com/randel/mvMISE/master/R/mvmise_b.r")
      res = mvMISE_b(Y, X0, id = rep(1:36, 4), sigma_diff = T, maxIter = maxIter)
      
      res$beta[2]/res$se[2]
    }
    
    stat_lm0 <- foreach(i = 1:nsim, .errorhandling = "pass") %dopar% 
    {
      y0 = as.vector(apply(dat1[sample(1:nrow(X)), , i], 1, function(x) x[-1] - x[1]))  # log-transformed / normal data
      X0 = NULL
      for (j in 1:nrow(dat1[, , i])) X0 = rbind(X0, t(X[j, , -1]))
      lmfit = lm(y0 ~ X0 - 1)
      coef(lmfit)[2]/sqrt(vcov(lmfit)[2, 2])
    }
    
    pval_perm = rep(NA, 2)
    pval_perm[1] <- mean(sapply(stat_lm, function(x, stat0) mean(abs(stat0) >= abs(x)), stat0 = unlist(stat_lm0)) <= 0.05, na.rm = T)
    pval_perm[2] <- mean(sapply(stat_mixemm, function(x, stat0) mean(abs(stat0) >= abs(x)), 
                                stat0 = unlist(stat_mixemm0[sapply(stat_mixemm0, length) == 1])) <= 0.05, na.rm = T)
    pval_perm[3] <- mean(sapply(stat_mvMISE, function(x, stat0) mean(abs(stat0) >= abs(x)), 
                                stat0 = unlist(stat_mvMISE0[sapply(stat_mvMISE0, length) == 1])) <= 0.05, na.rm = T)
    
    print("permutation")
    print(pval_perm)
    
    # simulated data (dat) are too large for more than 1k replications
    return(list(dat = dat, res_mvmise = res_mvmise, res_mixemm_ave = res_mixemm_ave, 
                stat_lm0 = stat_lm0, stat_mvMISE0 = stat_mvMISE0, res_summary = res_summary,
                stat_mixemm0 = stat_mixemm0, res_lm_ave = res_lm_ave, pval_perm = pval_perm))
    
  } else return(list(dat = dat, res_mvmise = res_mvmise, res_mixemm_ave = res_mixemm_ave, 
                     res_lm_ave = res_lm_ave, res_summary = res_summary, time = c(t_lm[3], t_mixemm[3], t_mvmise[3]), 
                     nodes = getDoParWorkers()))
}




# input lists generated by sim_mvmise(), return summary results
summary_res = function(list_res, beta, corr) {
  
  res_mvmise = list_res$res_mvmise
  res_mixemm_ave = list_res$res_mixemm_ave
  res_lm_ave = list_res$res_lm_ave
  
  # mvMISE
  # check possible error
  (err <- sapply(res_mvmise, function(x) length(x) == 2))
  
  par_all = t(sapply(res_mvmise[which(!err)], function(x) c(x$iter, x$beta, x$phi, x$sigma2, as.numeric(x$tau))))
  
  if(corr=='b') {
    pval_all = sapply(res_mvmise[which(!err)], function(x) 2 * pnorm(-abs((x$beta[2] - beta[2])/sqrt(x$var[2, 2]))))
    
    pval_all_power = sapply(res_mvmise[which(!err)], function(x) 2 * pnorm(-abs((x$beta[2])/sqrt(x$var[2, 2]))))
    
    # stat_mvMISE = sapply(res_mvmise[which(!err)], function(x) x$beta[2]/sqrt(x$var[2, 2]))
  }
  
  if(corr=='e') {
    pval_all = sapply(res_mvmise[which(!err)], function(x) 2 * pnorm(-abs((x$beta[2] - beta[2])/(x$beta[2]/x$stat[2]))))
    
    pval_all_power = sapply(res_mvmise[which(!err)], function(x) 2 * pnorm(-abs((x$stat[2]))))
    
    # stat_mvMISE = sapply(res_mvmise[which(!err)], function(x) x$stat[2])
  }
  
  # mixEMM
  # check possible error
  (err <- sapply(res_mixemm_ave, function(x) length(x) == 2))
  par_all1 = t(sapply(res_mixemm_ave[which(!err)], function(x) c(x$iter, x$beta, x$phi, x$sigma2, as.numeric(x$tau))))
  pval_all1 = sapply(res_mixemm_ave[which(!err)], function(x) 2 * pnorm(-abs((x$beta[2] - beta[2])/sqrt(x$var[2, 2]))))
  stat_mixemm = sapply(res_mixemm_ave[which(!err)], function(x) x$beta[2]/sqrt(x$var[2, 2]))
  pval_all_power1 = sapply(res_mixemm_ave[which(!err)], function(x) 2 * 
                             pnorm(-abs((x$beta[2])/sqrt(x$var[2, 2]))))
  
  
  ## lm
  beta_lm = sapply(res_lm_ave, function(lmfit) coef(lmfit)[2])
  pval_lm = sapply(res_lm_ave, function(lmfit) 2 * pnorm(-abs(coef(lmfit)[2]/sqrt(vcov(lmfit)[2, 2]))))
  pval0_lm = sapply(res_lm_ave, function(lmfit) 2 * pnorm(-abs((coef(lmfit)[2] - beta[2])/sqrt(vcov(lmfit)[2, 2]))))
  stat_lm = sapply(res_lm_ave, function(lmfit) coef(lmfit)[2]/sqrt(vcov(lmfit)[2, 2]))
  
  # print(c('ave estimate', 'MSE', 'type I', 'power'))
  print(round(res <- c(biasI = c(mean(par_all1[, 2]) - beta[1], mean(par_all[, 2]) - beta[1]),
                       biasII = c(mean(beta_lm) - beta[2], mean(par_all1[, 3]) - beta[2], mean(par_all[, 3]) - beta[2]),
                       MSEI = c(mean((par_all1[, 2] - beta[1])^2), mean((par_all[, 2] - beta[1])^2)), 
                       MSEII = c(mean((beta_lm - beta[2])^2), mean((par_all1[, 3] - beta[2])^2), 
                                 mean((par_all[, 3] - beta[2])^2)), 
                       typeI = c(mean(pval0_lm <= 0.05), mean(pval_all1 <= 0.05), mean(pval_all <= 0.05)),
                       power = c(mean(pval_lm <= 0.05), mean(pval_all_power1 <= 0.05), mean(pval_all_power <= 0.05))), 3))
  
  return(res)
}



Fisher.test <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
  return(c(Xsq = Xsq, p.value = p.val))
}


# pathway analysis: permute each sample for mvMISE_e and then use Fisher's method

sim_mvmise_e_fisher = function(K, beta=c(10, .7), phi = c(0, -0.025), sigma0 = 1, sigma1 = 1, D =  matrix(1), nsim = 1000, X, 
                               Sigma, lambda = .5, beta_prop = 1, beta_norm = F, probit = F, maxIter = 1000) {
  
  ############################################################## simulate data
  
  Z = kronecker(diag(1),rep(1,4))
  
  if(length(phi)==2) phi0 = c(phi, 0) else phi0 = phi
  
  dat0 <- foreach(i = 1:nsim, .errorhandling = "pass") %dopar% {
    source('https://raw.githubusercontent.com/randel/mvMISE_simulation/master/simulation%20functions.r')
    sim_dat(nsim = 1, X_real = X, K = K, beta = beta, phi = phi0, 
            sigma0 = sigma0, sigma1 = sigma1, D = D, D_exp = 0, Sigma = Sigma, 
            probit = probit, remove_allMiss = F, correlatedError = T, seed = i, beta_prop=beta_prop, beta_norm=beta_norm)
  }
  
  dat = array(NA, dim = c(dim(dat0[[1]]$dat[[1]]$Y), nsim)) # nexp, ni, K, nsim
  for (i in 1:nsim) dat[, , , i] = dat0[[i]]$dat[[1]]$Y
  
  dat1 = array(NA, dim = c(dim(dat[, , , 1])[-3], nsim)) # nexp, ni, nsim
  for (i in 1:nsim) dat1[, , i] = apply(dat[, , , i], 1:2, mean, na.rm = T)
  
  X0 = apply(X, 2, as.vector)
  
  ############################################################## mvMISE_e - Fisher
  
  t_mvmise = system.time(res_mvmise <- foreach(i=1:nsim, .errorhandling = 'pass') %dopar% {
    library(lme4)
    library(MASS)
    source('https://raw.githubusercontent.com/randel/mvMISE_simulation/master/simulation%20functions.r')
    
    Y = apply(dat[, , , i], 3, as.vector)
    # stack X across outcomes
    X_mat = X0[rep(1:nrow(X0), K), ]
    # Y_ind is the indicator matrix corresponding to different outcomes
    Y_ind = kronecker(diag(K), rep(1, nrow(Y)))
    # generate outcome-specific covariates
    cidx = 2 # the index for the covariate with outcome-specific coefficient
    X_mat = cbind(1, X_mat[, cidx] * Y_ind)
    source("https://raw.githubusercontent.com/randel/mvMISE/master/R/mvmise_e.r")
    res0 = mvMISE_e(Y, X_mat, id = rep(1:36, 4), sigma_diff = T, lambda = lambda, maxIter = maxIter)
    pval0 = 2*pnorm(-abs(res0$stat[-1]))
    pval = Fisher.test(pval0)
    
    # permutation
    set.seed(i+1)
    Y_perm = apply(dat[sample(1:nrow(dat)), , , i], 3, as.vector)
    res_perm = mvMISE_e(Y_perm, X_mat, id = rep(1:36, 4), sigma_diff = T, lambda = lambda, maxIter = maxIter)
    pval0 = 2*pnorm(-abs(res_perm$stat[-1]))
    stat_perm = Fisher.test(pval0)[1]
    
    return(list(res=res0, stat=pval[1], pval=pval[2], stat_perm = stat_perm))
  })
  
  err <- sapply(res_mvmise, function(x) length(x)==2)
  if(sum(err)>0) {print(c('mvmise', sum(err))); print(res_mvmise[[which(err)[1]]])}
  stat_mvmise = sapply(res_mvmise[which(!err)], function(x) x$stat)
  stat_perm_mvmise = sapply(res_mvmise[which(!err)], function(x) x$stat_perm)
  
  ############################################################## mixEMM-Fisher
  
  t_mixemm = system.time(res_mixemm <- foreach(i=1:nsim, .errorhandling = 'pass') %dopar% {
    
    library(lme4)
    library(MASS)
    source("https://raw.githubusercontent.com/randel/mvMISE/master/R/mvmise_b.r")
    source('https://raw.githubusercontent.com/randel/mvMISE_simulation/master/simulation%20functions.r')
    res = list()
    pval0 = pval0_perm = rep(NA, K)
    set.seed(i+1)
    sample_id = sample(1:nrow(dat))
    for(k in 1:K) {
      res[[k]] = mvMISE_b(apply(dat[,,,i][,,k,drop=F], 3, as.vector), X0, id = rep(1:36, 4), sigma_diff = T, maxIter = maxIter)
      pval0[k] = res[[k]]$pval[-1]
      
      # permutation
      tmp = mvMISE_b(apply(dat[sample_id,,,i][,,k,drop=F], 3, as.vector), X0, id = rep(1:36, 4), sigma_diff = T, maxIter = maxIter)
      pval0_perm[k] = tmp$pval[-1]
    }
    pval = Fisher.test(pval0)
    stat_perm = Fisher.test(pval0_perm)[1]
    
    return(list(res=res, stat=pval[1], pval=pval[2], stat_perm = stat_perm))
  })
  
  err <- sapply(res_mixemm, function(x) length(x)==2)
  if(sum(err)>0) {print(c('mixemm', sum(err)))} #; print(res_mixemm[[which(err)[1]]])}
  
  stat_mixemm = sapply(res_mixemm[which(!err)], function(x) x$stat)
  stat_perm_mixemm = sapply(res_mixemm[which(!err)], function(x) x$stat_perm)
  
  ############################################################## lm-Fisher
  
  X0 = NULL
  for(j in 1:nrow(X)) X0 = rbind(X0, t(X[j,,-1]))
  
  t_lm = system.time(res_lm <- foreach(i=1:nsim, .errorhandling = 'pass') %dopar% {
    
    source('https://raw.githubusercontent.com/randel/mvMISE_simulation/master/simulation%20functions.r')
    
    pval0 = pval0_perm = rep(NA, K)
    set.seed(i+1)
    sample_id = sample(1:nrow(dat))
    
    for(j in 1:K) {
      y0 = as.vector(apply(dat[,,j,i], 1, function(x) x[-1]-x[1])) # log-transformed / normal data
      lmfit = lm(y0 ~ X0-1)
      pval0[j] = 2*pnorm(-abs(coef(lmfit)[2] / sqrt(vcov(lmfit)[2,2])))
      
      # permutation
      y0 = as.vector(apply(dat[sample_id,,j,i], 1, function(x) x[-1]-x[1])) # log-transformed / normal data
      lmfit = lm(y0 ~ X0-1)
      pval0_perm[j] = 2*pnorm(-abs(coef(lmfit)[2] / sqrt(vcov(lmfit)[2,2])))
    }
    
    return(list(pval = Fisher.test(pval0)[2], stat = Fisher.test(pval0)[1], stat_perm = Fisher.test(pval0_perm)[1]))
  })
  
  err <- sapply(res_lm, function(x) length(x)==2)
  if(sum(err)>0) {print(c('lm', sum(err))); print(res_lm[[which(err)[1]]])}
  
  stat_perm_lm = sapply(res_lm[which(!err)], function(x) x$stat_perm) # check errors
  stat_lm = sapply(res_lm[(!err)], function(x) x$stat)
  sig_perm_lm = mean(sapply(stat_lm, function(x) mean(abs(stat_perm_lm)>=abs(x)))<=.05, na.rm=T)
  
  sig_perm_mixemm <- mean(sapply(stat_mixemm, function(x) mean(abs(stat_perm_mixemm)>=abs(x)))<=.05, na.rm=T)
  
  sig_perm_mvmise <- mean(sapply(stat_mvmise, function(x) mean(abs(stat_perm_mvmise)>=abs(x)))<=.05, na.rm=T)
  
  print(c(K=K, probit = probit, lm=sig_perm_lm, mixemm=sig_perm_mixemm, mvmise=sig_perm_mvmise))
  
  return(list(sig_rate = c(sig_perm_lm, sig_perm_mixemm, sig_perm_mvmise), res_mvmise=res_mvmise, stat_perm_mvmise=stat_perm_mvmise, dat = dat,
              res_lm = res_lm, res_mixemm = res_mixemm, time = c(t_lm[3], t_mixemm[3], t_mvmise[3])))
}



## Power-law degree distribution
# http://igraph.org/r/doc/sample_degseq.html
## Note that we correct the degree sequence if its sum is odd

sim_Sigma_power_law = function(K, PLOT=F) {
  library(igraph)
  library(MASS)
  
  set.seed(2)
  degs <- sample(1:K, K, replace=TRUE, prob=(1:K)^(-2.3))
  if (sum(degs) %% 2 != 0) { degs[1] <- degs[1] + 1 }
  g5 <- degree.sequence.game(degs, method="vl")
  # all(degree(g5) == degs)
  if(PLOT) plot(g5)
  
  # Peng et al. (2009)
  A = as.matrix(get.adjacency(g5))
  diag(A) = 1
  
  for(i in 1:K) {
    for(j in 1:K) {
      if(i!=j & A[i,j]!=0) A[i,j] = runif(1, .1, .4)*sample(c(-1,1), 1)
    }
  }
  
  for(i in 1:K) A[i,-i] = A[i,-i] / (2*sum(abs(A[i,-i]))) # force positive definiteness
  A = (A + t(A))/2
  A_inv = ginv(A)
  
  Sigma = diag(K)
  for(i in 1:K) {
    for(j in 1:K) {
      if(i!=j) Sigma[i,j] = A_inv[i,j] / sqrt(A_inv[i,i]*A_inv[j,j]) # make sure all diag A_inv positive
    }
  }
  
  return(Sigma)
}
