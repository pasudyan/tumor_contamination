###############################################################
## Script to run on Merged Data with Matrices                ##
## Model :                                                   ##
## - Y_i ~ Bin(N_i, zeta_i*theta_i^n + (1-zeta_i)*theta_i^t) ##
## where :                                                   ##
## - zeta_i ~DP(alpha, Beta(1,1))                            ##
## - theta_i^n ~ Dirichlet(1,1,1), theta in {0, 0.5, 1}      ##
###############################################################

rm(list = ls())

source("./mh_algoZetaMat_Sim.R")
source("./dataGen.R")
source("./library.R")

# set.seed(1234)

################ Use simulated Data #################
kk = 3
num_obs  = 1500
zeta_sym = 0.1
diffzeta = c(0.1, 0.9)
diffN_i  = 100
theta    = c(0, 0.5, 1)
datta    = dataGen(kk, num_obs, theta, diffN_i, diffzeta)

data1 = datta$oridata1[c('trials', 'mixDat')]
data2 = datta$oridata2[c('trials', 'mixDat')]

############### Parameter initialization ##############
alpha     = 100

alpha_base = 0 #parameters of base distribution
beta_base  = 1
m          = 2 #number of auxiliary params
numIter    = 8000 
burnIn     = 3000 
telliter   = 500 
xtraSpace  = 2*num_obs
alpha_0    = c(1, 1, 1)

#initialize number of clusters
K = num_obs

#initialize cluster assignment for zeta 
data1$cluster = sample(c(1:K, sample(1:K, num_obs-K, replace=TRUE)), 
                       num_obs, replace=FALSE)
data2$cluster = sample(c(1:K, sample(1:K, num_obs-K, replace=TRUE)),
                       num_obs, replace=FALSE)

# data1$cluster = sapply(1:num_obs, function(x){
#   if (datta$oridata1$zeta[x] == diffzeta[1])
#     1
#   else 
#     2
# })
# 
# data2$cluster = sapply(1:num_obs, function(x){
#   if (datta$oridata2$zeta[x] == diffzeta[2])
#     1
#   else 
#     2
# })

# data1$cluster = sample(
#   c(1:K, sample(
#     1:K, num_obs-K, replace=TRUE)
#     ), 
#   num_obs, replace=FALSE 
#   )
# data2$cluster = sample(
#   c(1:K, sample(
#     1:K, num_obs-K, replace=TRUE)
#     ),
#   num_obs, replace=FALSE
#   )

#initialize cluster assignment for theta
theta = c(0, 0.5, 1)

# u_n = sapply(1:num_obs, function(x){
#   if (datta$oridata1$theta_n[x] == 0)
#     1
#   else if (datta$oridata1$theta_n[x] == 0.5)
#     2
#   else
#     3
# })
# 
# u_t = sapply(1:num_obs, function(x){
#   if (datta$oridata1$theta_t[x] == 0)
#     1
#   else if (datta$oridata1$theta_t[x] == 0.5)
#     2
#   else
#     3
# })

data1$ucl_n = sample(1:3, num_obs, replace = TRUE)
data1$ucl_t = sample(1:3, num_obs, replace = TRUE)

data2$ucl_n = data1$ucl_n 
data2$ucl_t = data1$ucl_t 

rho = list()
rho[[1]] = matrix(NA, nrow = num_obs, ncol = numIter)
rho[[2]] = matrix(NA, nrow = num_obs, ncol = numIter)

#initialize zeta for observations in each cluster
zeta1 = runif(K, alpha_base, beta_base)
zeta2 = runif(K, alpha_base, beta_base)

datanum = list()
datanum[[1]] = list(data   = as.matrix(data1), 
                    zeta   = zeta1)
datanum[[2]] = list(data   = as.matrix(data2),
                    zeta   = zeta2)

store_zeta  = list(zeta1 = matrix(NA, 
                                  nrow = num_obs, 
                                  ncol = numIter),
                   zeta2 = matrix(NA, 
                                  nrow = num_obs, 
                                  ncol = numIter))
store_theta_n = matrix(NA, nrow = num_obs, ncol = numIter)
store_theta_t = matrix(NA, nrow = num_obs, ncol = numIter)

#Rprof("iteration.out")
system.time(
for (ii in 1:numIter){
  
  #For loop to update zeta
  for (d in 1:2){
    dataMat   = datanum[[d]]$data
    zeta      = numeric(xtraSpace)
    K         = length(datanum[[d]]$zeta)
    zeta[1:K] = datanum[[d]]$zeta
    theta_n = theta[dataMat[,"ucl_n"]]
    theta_t = theta[dataMat[,"ucl_t"]]
    
    counts  = numeric(xtraSpace)
    counts[1:K]  = sapply(1:K, function(r) 
      sum(dataMat[,"cluster"] == r))
    prob_clust   = numeric(xtraSpace)
    rand_samples = sample(1:num_obs,
                          num_obs,
                          replace=FALSE)
    
#     if (d == 1){
#       alpha_base = 0
#       beta_base  = 1
#     } else {
#       alpha_base = 0
#       beta_base  = 1
#     }
#     
    #iteration for each observation
    for (b in 1:num_obs){
      
      j = rand_samples[b]
      h = K + m
      
      cl = dataMat[,"cluster"]
      
      #removing the observation from the cluster
      counts[cl[j]] = counts[cl[j]]-1
      
      #reevaluating cluster assignments
      if (counts[cl[j]]!=0){
        
        #drawing values for the aux params from G_0
        zeta[(K+1):h] = runif(m, alpha_base, beta_base) 
        
      } else if (counts[cl[j]]==0 & counts[cl[j]]!=K){
        
        #renumbering the clusters
        counts[cl[j]] = counts[K]
        cl[cl == K] = cl[j]
        counts[K] = 0
        
        #renumbering the zeta
        temp = zeta[cl[j]]
        zeta[cl[j]] = zeta[K]
        zeta[K] = temp
        
        #reducing the number of clusters and aux clusters
        K = K-1
        h = h-1
        
        #drawing values for the aux params from G_0
        zeta[(K+2):h] = runif(m-1, alpha_base, beta_base)
        
      } else if (counts[cl[j]]==0 & cl[j]==K) {
        
        #reducing the number of clusters and aux clusters
        K = K-1
        h = h-1
        
        #drawing values for the aux params from G_0
        zeta[(K+2):h] = runif((m-1), alpha_base, beta_base)
        
      }
      
      #prob of choosing existing cluster 
      prob_clust[1:K] = log(counts[1:K]) +
        dbinom(dataMat[j, "mixDat"], 
               dataMat[j, "trials"],
               zeta[1:K]*theta_n[j] + (1-zeta[1:K])*theta_t[j],
               log=TRUE)
      
      #prob of choosing new cluster
      prob_clust[(K+1):h] = log(alpha/m) + 
        dbinom(dataMat[j, "mixDat"],
               dataMat[j, "trials"],
               zeta[(K+1):h]*theta_n[j] + (1-zeta[(K+1):h])*theta_t[j],
               log=TRUE)
      
      # prob_clust = ifelse(prob_clust == -Inf, -100000, prob_clust)
      prob_clust[prob_clust == -Inf] <- -100000
      
      #normalizing constant
      prob_norm = prob_clust[1:h] -
        logSumExp(prob_clust[1:h])
      
      #sampling new cluster assignments
      new_cl  = sample(1:h,
                       1,
                       replace=FALSE,
                       exp(prob_norm))
      
      #new table addition
      if (new_cl > K){
        cl[j] = K+1
        zeta[K+1]   = zeta[new_cl]
        counts[K+1] = 1
        K = K+1
      } else {
        cl[j] = new_cl 
        counts[new_cl] = counts[new_cl]+1
      }
      
      zeta[(K+1):xtraSpace]  = rep(0,xtraSpace-K)
      counts[(K+1):xtraSpace]= rep(0,xtraSpace-K)
      dataMat[,"cluster"] = cl
      
    } #end of observations loop 
    
    #sampling the new parameters using the MH Algorithm 
    zeta[1:K] = mh_algoZetaMat(dataMat, zeta[1:K], theta, d)
    
    datanum[[d]] = list(data   = dataMat, 
                        zeta   = zeta[1:K])

    #storing params value for observation 
    store_zeta[[d]][,ii] = zeta[dataMat[,"cluster"]]
    
  } #end of data loop
  
  #For loop to update theta
  u_Update = mh_cluster2(datanum[[1]], datanum[[2]], theta,
                         alpha_0)
  datanum[[1]]$data = u_Update$dataMat1
  datanum[[2]]$data = u_Update$dataMat2
  
  store_theta_n[,ii] = u_Update$dataMat1[,"ucl_n"]
  store_theta_t[,ii] = u_Update$dataMat1[,"ucl_t"]
  
  rho[[1]][,ii] = rho_calc(datanum[[1]], theta)
  rho[[2]][,ii] = rho_calc(datanum[[2]], theta)
  
} #end of iteration loop
)

# Rprof(NULL)
# rprof_res = summaryRprof("iteration.out")$by.total
# print(rprof_res)

# source("./postPlots_Sim.R")
