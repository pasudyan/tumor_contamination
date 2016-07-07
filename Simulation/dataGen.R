############### Generating mix data #####################
## Y_i ~ Bin(N_i, zeta_i*theta_n + (1-zeta_i)*theta_t) ##
## Input Arguments:                                    ##
## - K = number of theta cluster                       ##
## - num_obs = number of observations                  ##
## - theta   = vector of unique theta                  ##
## - diffN_i  = vector of unique number of trials (N_i)##
## - diffzeta = vector of unique zeta                  ##
##                                                     ##
## Output Arguments:                                   ##
## - List containing data set 1 and 2                  ##
#########################################################

dataGen = function(K, num_obs, theta, diffN_i, diffzeta){
  zeta1   = rep(diffzeta[1], num_obs)
  zeta2   = rep(diffzeta[2], num_obs)
#   zeta1   = sample(diffzeta, num_obs, replace = TRUE)
#   zeta2   = sample(diffzeta, num_obs, replace = TRUE)
  theta_n = sample(theta, num_obs, replace=TRUE)
  theta_t = sample(theta, num_obs, replace=TRUE)
  N_i1    = sample(diffN_i, num_obs, replace=TRUE)
  N_i2    = sample(diffN_i, num_obs, replace=TRUE)
  mixDat1 = rbinom(num_obs, N_i1, zeta1*theta_n + (1-zeta1)*theta_t)
  mixDat2 = rbinom(num_obs, N_i2, zeta2*theta_n + (1-zeta2)*theta_t)
  rho1    = zeta1*theta_n + (1-zeta1)*theta_t
  rho2    = zeta2*theta_n + (1-zeta2)*theta_t
  data1   = data.frame(trials  = N_i1,
                       mixDat  = mixDat1,
                       zeta    = zeta1,
                       theta_n = theta_n,
                       theta_t = theta_t,
                       rho     = rho1)
  data2   = data.frame(trials  = N_i2, 
                       mixDat  = mixDat2,
                       zeta    = zeta2,
                       theta_n = theta_n,
                       theta_t = theta_t,
                       rho     = rho2)
  return(list(oridata1 = data1, oridata2 = data2))
}


