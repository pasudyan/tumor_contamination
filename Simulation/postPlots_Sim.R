## Script for Visualization of Posterior Density ## 
###################################################

## Label data for unidentifiability ##
label_id = which(datta$oridata1$theta_n != datta$oridata1$theta_t)

## Plot posterior density of Zeta ##
test_obs1    = sample(label_id, 10, replace = FALSE)
zeta_test_1  = store_zeta[[1]][test_obs1,(burnIn+1):numIter]
zeta_test_2  = store_zeta[[2]][test_obs1,(burnIn+1):numIter]

## Plot means of all zeta's 
mean_all_zeta_1 = melt(apply(store_zeta[[1]][label_id,(burnIn+1):numIter],
                             1,mean))
mean_all_zeta_2 = melt(apply(store_zeta[[2]][label_id,(burnIn+1):numIter],
                             1,mean))

b1 = qplot(value, data = mean_all_zeta_1, geom="histogram", 
           binwidth = 0.005) + 
  ggtitle("Histogram of Mean Zetas for Data Set 1")

jpeg(file = "plotmeanZeta1_01_09.jpeg")
print(b1)
dev.off()

b2 = qplot(value, data = mean_all_zeta_2, geom="histogram", binwidth = 0.005) +
  ggtitle("Histogram of Mean Zetas for Data Set 2")

jpeg(file = "plotmeanZeta2_01_09.jpeg")
print(b2)
dev.off()

## Plot means of all rho's
mean_all_rho_1 = melt(apply(rho[[1]][label_id,(burnIn+1):numIter],
                            1,mean))
mean_all_rho_2 = melt(apply(rho[[2]][label_id,(burnIn+1):numIter],
                            1,mean))

b = qplot(value, data = mean_all_rho_1, geom="histogram", binwidth = 0.01)+
  ggtitle("Histogram of Mean Rhos for Data Set 1")

#jpeg(file = "plotmeanRho1_01_09.jpeg")
print(b)
#dev.off()

b = qplot(value, data = mean_all_rho_2, geom="histogram", binwidth = 0.01)+
  ggtitle("Histogram of Mean Rhos for Data Set 2")

#jpeg(file = "plotmeanRho2_01_09.jpeg")
print(b)
#dev.off()

## Plot posterior density of theta ##
unlist(lapply(1:length(test_obs1), function(x) {
  store_theta_n[test_obs1[x], (burnIn+1):numIter]
})) %>% 
  matrix(nrow=(numIter-burnIn), ncol=10) ->
  dens_theta_n

unlist(lapply(1:length(test_obs1), function(x) {
  store_theta_t[test_obs1[x], (burnIn+1):numIter]
})) %>% 
  matrix(nrow=(numIter-burnIn), ncol=10) ->
  dens_theta_t

name = sapply(1:length(test_obs1), function(x){paste("Obs ",x)}) 
colnames(dens_theta_n) = name
colnames(dens_theta_t) = name

lapply(1:10, function(i) {
  data.frame(
    obs = i,
    theta_n = store_theta_n[test_obs1[i], (burnIn+1):numIter],
    theta_t = store_theta_t[test_obs1[i], (burnIn+1):numIter]
  ) %>%
    group_by(obs, theta_n, theta_t) %>%
    tally()
}) %>% 
  do.call(rbind, .) ->
  categories

categories$theta_n_ <- paste("Normal = ",theta[categories$theta_n], 
                             sep = "")
categories$theta_t_ <- paste("Tumor = ", theta[categories$theta_t], 
                             sep = "")

b = ggplot(data = categories) + 
  theme_bw() + 
  geom_bar(aes(x = as.factor(obs), y = n, fill = as.factor(obs)), stat = "identity") + 
  facet_grid(theta_t_ ~ theta_n_) + 
  scale_fill_brewer(palette = "Set3") + 
  labs(x = "Observation", y = "Count", fill = "Observation") +
  ggtitle("Joint Cluster Assignment Frequencies for Theta Normal and Tumor")

jpeg(file = "plottheta_01_09.jpeg")
print(b)
dev.off()

### Create violin plot of zeta distribution ###
zeta_test_1df = as.data.frame(t(zeta_test_1))
zeta_test_2df = as.data.frame(t(zeta_test_2))

names(zeta_test_1df) = test_obs1
names(zeta_test_2df) = test_obs1
zeta_test_1_melt = melt(zeta_test_1df)
zeta_test_2_melt = melt(zeta_test_2df)

b1 = qplot(factor(variable), value, data=zeta_test_1_melt, 
      geom="violin") +
  stat_summary(fun.y = mean, geom="point", shape = 23, size = 2) +
  ggtitle("Posterior Distribution of Zeta for Data 1") +
  ylab("Value") +
  xlab("Observations")

jpeg(file = "plotzeta1_01_09.jpeg")
print(b1)
dev.off()

b2 = qplot(factor(variable), value, data=zeta_test_2_melt, 
      geom="violin") + 
  stat_summary(fun.y = mean, geom="point", shape = 23, size = 1) +
  ggtitle("Posterior Distribution of Zeta for Data 2") +
  ylab("Value") +
  xlab("Observations")

jpeg(file = "plotzeta2_01_09.jpeg")
print(b2)
dev.off()

### Plot for Rho ###
rho_1_df = as.data.frame(t(rho[[1]][1:length(test_obs1), 
                                    (burnIn+1):numIter]))
rho_2_df = as.data.frame(t(rho[[2]][1:length(test_obs1), 
                                    (burnIn+1):numIter]))

names(rho_1_df) = names(rho_2_df) = test_obs1

rho_1_melt = melt(rho_1_df)
rho_2_melt = melt(rho_2_df)

b = ggplot(rho_1_melt, aes(factor(variable), value)) +
  geom_violin(scale = "width") +
  stat_summary(fun.y = mean, geom="point", shape = 23, size = 1) +
  ggtitle("Posterior of Rho for Data 1") + 
  ylab("Value") + 
  xlab("Observations")

jpeg(file = "plotrho1_01_09.jpeg")
print(b)
dev.off()

b = ggplot(rho_2_melt, aes(factor(variable), value)) +
  geom_violin(scale = "width") +
  stat_summary(fun.y = mean, geom="point", shape = 23, size = 1) + 
  ggtitle("Posterior of Rho for Data 2") + 
  ylab("Value") + 
  xlab("Observations")

jpeg(file = "plotrho2_01_09.jpeg")
print(b)
dev.off()

############ Tabulating the estimates ##############
theta_n_hat = theta[sapply(1:10, function(x) 
  which.max(tabulate(dens_theta_n[label_id,x])))]
theta_t_hat = theta[sapply(1:10, function(x) 
  which.max(tabulate(dens_theta_t[label_id,x])))]

mean_zeta_1 = round(apply(zeta_test_1df,2,mean),4)
mean_zeta_2 = round(apply(zeta_test_2df,2,mean),4)

mean_rho_1 = round(apply(rho_1_df,2,mean),4)
mean_rho_2 = round(apply(rho_2_df,2,mean),4)

sd_rho_1 = round(apply(rho_1_df,2,sd),4)
sd_rho_2 = round(apply(rho_2_df,2,sd),4)

est_data1 = cbind(theta_n_hat, theta_t_hat, 
                  mean_zeta_1, mean_rho_1, sd_rho_1)
est_data2 = cbind(theta_n_hat, theta_t_hat, 
                  mean_zeta_2, mean_rho_2, sd_rho_2)

print(est_data1)
print(est_data2)

est_data = list(est_data1, est_data2)

################## Number of clusters ####################
unique_zeta <- lapply(1:2, function(x){
  sapply(1:numIter, function(s){
    length(unique(store_zeta[[x]][,s]))})
})

b = qplot(unique_zeta[[1]], geom="histogram", binwidth = 1) + 
  xlab("Number of unique clusters") +
  ggtitle("Histogram of Number of Cluster for Data Set 1")

#jpeg(file = "chrm1_0.01_plotnumclust1.jpeg")
print(b)
#dev.off()
  
qplot(unique_zeta[[2]], geom="histogram", binwidth = 1) + 
  xlab("Number of unique clusters") +
  ggtitle("Histogram of Number of Cluster for Data Set 2")

#jpeg(file = "chrm1_0.01_plotnumclust2.jpeg")
print(b)
#dev.off()

####### For Simulated Data #######
cbind(datta$oridata1[test_obs1,], est_data1) %>% print()
cbind(datta$oridata2[test_obs1,], est_data2) %>% print()

###### Calculate error rates #####
error_zeta_1 = sqrt((mean_all_zeta_1 - (1-diffzeta[1]))^2)
error_zeta_2 = sqrt((mean_all_zeta_2 - (1-diffzeta[2]))^2)

error_rho_1 = sqrt((mean_all_rho_1 - datta$oridata1$rho[label_id])^2)
error_rho_2 = sqrt((mean_all_rho_2 - datta$oridata2$rho[label_id])^2)

## Plotting the errors
jpeg(file = "plotErrorZeta1_01_09.jpeg")
q = plot(sqrt(error_zeta_1), type = "l", ylab = "Error", main = "Error for Zeta Data Set 1")
print(q)
dev.off()

jpeg(file = "plotErrorZeta2_01_09.jpeg")
q = plot(sqrt(error_zeta_2), type = "l", ylab = "Error", main = "Error for Zeta Data Set 2")
print(q)
dev.off()

jpeg(file = "plotErrorRho1_01_09.jpeg")
q = plot(sqrt(error_rho_1), type = "l", ylab = "Error", main = "Error for Rho Data Set 1")
print(q)
dev.off()

jpeg(file = "plotErrorRho2_01_09.jpeg")
q = plot(error_rho_2, type = "l", ylab = "Error", main = "Error for Rho Data Set 2")
print(q)
dev.off()

## Quantifying error rate
mse_zeta1 = sum((mean_all_zeta_1 - (1-diffzeta[1]))^2)/(num_obs - 1)
mse_zeta2 = sum((mean_all_zeta_2 - (1-diffzeta[2]))^2)/(num_obs - 1)

mse_rho1  = sum((mean_all_rho_1 - datta$oridata1$rho[label_id])^2)/
  (num_obs - 1)
mse_rho2  = sum((mean_all_rho_1 -  datta$oridata1$rho[label_id])^2)/
  (num_obs - 1)

unlist(lapply(1:num_obs, function(x) {
  store_theta_n[test_obs1[x], (burnIn+1):numIter]
})) %>% 
  matrix(nrow=(numIter-burnIn), ncol=10) ->
  dens_theta_n

mode_theta_n_hat = theta[unlist(lapply(1:num_obs, function(x){
  a <- count(store_theta_n[x,(burnIn +1):numIter])
  a$x[which.max(a$freq)]
}))]

mode_theta_t_hat = theta[unlist(lapply(1:num_obs, function(x){
  a <- count(store_theta_t[x,(burnIn +1):numIter])
  a$x[which.max(a$freq)]
}))]

error_dif_theta_n <- mode_theta_n_hat - datta$oridata1$theta_t
error_dif_theta_t <- mode_theta_t_hat - datta$oridata1$theta_n

jpeg(file = "plotErrorThetaN_01_09.jpeg")
q = plot(error_dif_theta_n, type = "l", ylab = "Error", main = "Error for Theta_n")
print(q)
dev.off()

jpeg(file = "plotErrorThetaT_01_09.jpeg")
q = plot(error_dif_theta_t, type = "l", ylab = "Error", main = "Error for Theta_t")
print(q)
dev.off()

error_rate_theta_n <- length(which(
  (mode_theta_n_hat - datta$oridata1$theta_t) != 0))
error_rate_theta_t <- length(which(
  (mode_theta_t_hat - datta$oridata1$theta_n) != 0))

print(error_rate_theta_n)
print(error_rate_theta_t)
print(cbind(mse_zeta1, mse_zeta2))
print(cbind(mse_rho1, mse_rho2))

mcmc(t(store_theta_n), start=burnIn+1, end=numIter) %>% 
  effectiveSize() -> eff_size_theta_n
print(eff_size_theta_n)
mcmc(t(store_theta_t), start=burnIn+1, end=numIter) %>% 
  effectiveSize() -> eff_size_theta_t
print(eff_size_theta_t)

