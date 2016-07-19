###############################################################
## Functions to run Gibbs sampling with Matrices             ##
## Model :                                                   ##
## - Y_i ~ Bin(N_i, zeta_i*theta_i^n + (1-zeta_i)*theta_i^t) ##
## where :                                                   ##
## - zeta_i ~DP(alpha, Beta(1,1))                            ##
## - theta_i^n ~ Dirichlet(1,1,1), theta in {0, 0.5, 1}      ##
###############################################################


#Metropolis-Hastings algorithm for Binomial-Beta
mh_algoZetaMat <- function(dataMat, zeta, theta, ds){
  zetaCurr    = zeta
  alpha_base  = 0
  beta_base   = 0.5
  mhloop      = 1
  
  dataMat_up = cbind(dataMat, theta_n = theta[dataMat[,"ucl_n"]], 
                     theta_t = theta[dataMat[,"ucl_t"]])
  
  ratio = numeric(length(zetaCurr))
  params = numeric(length(zetaCurr))
  
  for(t in 1:mhloop){
    zetaCand = runif(length(zetaCurr),alpha_base, beta_base)
    
    #Calculating the acceptance probability 
    for (i in 1:length(zetaCurr)){
      mydatta  = dataMat_up[which(dataMat_up[, "cluster"] == i),]
      
      if (class(mydatta) == "numeric")
        mydatta = as.matrix(t(mydatta))
      
      if (ds == 1){
        alpha_base = 0
        beta_base  = 0.5
        
        zetaCand = runif(length(zetaCurr),alpha_base, beta_base)
        
        n_vec = sum(dbinom(mydatta[,"mixDat"],
                           mydatta[,"trials"],
                           zetaCand[i]*mydatta[,"theta_n"] + 
                             (1-zetaCand[i])*mydatta[,"theta_t"], 
                           log = TRUE), 
                    dunif(zetaCand[i], 
                          alpha_base,
                          beta_base,
                          log = TRUE))
        
        d_vec = sum(dbinom(mydatta[,"mixDat"],
                           mydatta[,"trials"],
                           zetaCurr[i]*mydatta[,"theta_n"]+ 
                             (1-zetaCurr[i])*mydatta[,"theta_t"],
                           log = TRUE), 
                    dunif(zetaCurr[i],
                          alpha_base,
                          beta_base,
                          log = TRUE)) 
        
        
      } else if (ds == 2){
        alpha_base = 0.5
        beta_base  = 1
        
        zetaCand = runif(length(zetaCurr),alpha_base, beta_base)
        
        n_vec = sum(dbinom(mydatta[,"mixDat"],
                           mydatta[,"trials"],
                           zetaCand[i]*mydatta[,"theta_n"] + 
                             (1-zetaCand[i])*mydatta[,"theta_t"], 
                           log = TRUE), 
                    dunif(zetaCand[i], 
                          alpha_base,
                          beta_base,
                          log = TRUE))
        
        d_vec = sum(dbinom(mydatta[,"mixDat"],
                           mydatta[,"trials"],
                           zetaCurr[i]*mydatta[,"theta_n"]+ 
                             (1-zetaCurr[i])*mydatta[,"theta_t"],
                           log = TRUE), 
                    dunif(zetaCurr[i],
                          alpha_base,
                          beta_base,
                          log = TRUE)) 
      }
      
      n_vec[n_vec == -Inf] <- -100000
      d_vec[d_vec == -Inf] <- -100000
      
      ratio = n_vec-d_vec
      
      params[i] = if(log(runif(1)) < ratio) zetaCand[i] else zetaCurr[i]
      
    }
    zetaCurr = params 
  }
  return(zetaCurr)
}

############# Function for cluster assignments ################
mh_cluster2 = function(datanum1, datanum2, theta, alpha_0){
  K = 9
  k = 3
  dataMat1 = datanum1$data
  dataMat2 = datanum2$data
  zeta1    = datanum1$zeta[dataMat1[, "cluster"]]
  zeta2    = datanum2$zeta[dataMat2[, "cluster"]]
  theta_n  = theta[dataMat1[,"ucl_n"]]
  theta_t  = theta[dataMat1[,"ucl_t"]]
  num_obs  = nrow(dataMat1)
  
  counts_n = sapply(1:k, function(r) sum(dataMat1[, "ucl_n"] == r))
  counts_t = sapply(1:k, function(r) sum(dataMat1[, "ucl_t"] == r))
  
  cluster_n = dataMat1[,"ucl_n"]
  cluster_t = dataMat1[,"ucl_t"]
  
  thetaRep_n = rep(theta, times = c(3,3,3))
  thetaRep_t = rep(theta, 3)
  
  for (j in 1:num_obs){
    
    counts_n[cluster_n[j]] = counts_n[cluster_n[j]] - 1
    counts_t[cluster_t[j]] = counts_t[cluster_t[j]] - 1
    
    #Posterior dirichlet distribution
    alphaPrime_n = alpha_0 + counts_n
    alphaPrime_t = alpha_0 + counts_t
    
    alphaRep_n = rep(alphaPrime_n, times = c(3,3,3)) 
    alphaRep_t = rep(alphaPrime_t, 3)
    
    #Compute cluster probabilities
    prob_clust = log(alphaRep_n) + log(alphaRep_t) + 
        dbinom(dataMat1[j, "mixDat"], dataMat1[j, "trials"], 
               (zeta1[j])*thetaRep_n + 
                 (1-zeta1[j])*thetaRep_t, log = TRUE)+  
        dbinom(dataMat2[j,"mixDat"], dataMat2[j,"trials"], 
               (zeta2[j])*thetaRep_n + 
                 (1-zeta2[j])*thetaRep_t, log = TRUE) 
    
    prob_clust[prob_clust == -Inf] <- -100000
    
    #Includes normalizing 
    prob_norm = prob_clust - logSumExp(prob_clust)
    
    #Sample new cluster 
    new_c = sample(1:K, 1, FALSE, exp(prob_norm))
    
    #Update cluster
    cluster_n[j] = ceiling(new_c/3)
    cluster_t[j] = if(new_c %% 3 == 0) 3 else new_c %% 3

    #Update counts
    counts_n[cluster_n[j]] = counts_n[cluster_n[j]] + 1
    counts_t[cluster_t[j]] = counts_t[cluster_t[j]] + 1
    
  }
  
    dataMat1[, "ucl_n"] = dataMat2[, "ucl_n"] = cluster_n
    dataMat1[, "ucl_t"] = dataMat2[, "ucl_t"] = cluster_t
    
    return(list(dataMat1 = dataMat1, dataMat2 = dataMat2))
}

rho_calc = function(X, theta){
  
  a = X$zeta[X$data[,"cluster"]]*theta[X$data[,"ucl_n"]] + 
    (1-X$zeta[X$data[,"cluster"]])*theta[X$data[,"ucl_t"]]
  return(a)
}

rho_calc_ori = function(X, theta){
  a = X$zeta*theta_n + (1-X$zeta)*theta_t
  return(a)
}

