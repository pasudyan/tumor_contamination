########################################################
##  Script to read and merge data ready for analysis  ##
##  Strand Data 11_10_010                             ##
########################################################

num_chrm <- 2
col_names <- c("Chr", "Loc", "Fwd", "Backwd", "Depth", "Pval")

temp_dat <- lapply(1:num_chrm, function(x){
  file1 <- paste("C:/Users/pasudyan/OneDrive/PhD_School/Research/Cancer_Project/Data/Strand/11_10_mix_0.05chr",x,".out", sep="")
  file2 <- paste("C:/Users/pasudyan/OneDrive/PhD_School/Research/Cancer_Project/Data/Strand/11_10_mix_0.95chr",x,".out", sep="")
  
  raw_dat1 <- read.table(file1)
  raw_dat2 <- read.table(file2)
  
  names(raw_dat1) <- col_names
  names(raw_dat2) <- col_names
  
  raw_dat1$Var  <- raw_dat1$Fwd + raw_dat1$Backwd
  raw_dat2$Var  <- raw_dat2$Fwd + raw_dat2$Backwd
  
  raw_dat1$Freq <- raw_dat1$Var/raw_dat1$Depth
  raw_dat2$Freq <- raw_dat2$Var/raw_dat2$Depth
  
  #Subseting to data for frequencies between 0.1 and 0.9
  raw_sub1 <- subset(raw_dat1, raw_dat1$Freq < 0.9 & raw_dat1$Freq > 0.1
                     & raw_dat1$Pval > 0.001)
  raw_sub2 <- subset(raw_dat2, raw_dat2$Freq < 0.9 & raw_dat2$Freq > 0.1
                     & raw_dat2$Pval > 0.001)
  
  merge_dat <- merge(raw_sub1, raw_sub2, by = "Loc")
  
  return(merge_dat)
})

l_ply(1:num_chrm, function(x){
  data1 <- temp_dat[[x]][, c("Loc", "Depth.x", "Var.x", "Freq.x", "Pval.x")]
  data2 <- temp_dat[[x]][, c("Loc", "Depth.y", "Var.y", "Freq.y", "Pval.y")]
  file.name <- paste("11_10_0.05chrm",x,".RData",sep = "")
  merge_loc <- list(data1 = data1, data2 = data2) 
  save(merge_loc, file= file.name)
})

######## Merging all 22 Chromosomes #########
load("11_10_0.1chrm1.Rdata")
merge_chrm_dta1 <- merge_loc$data1
merge_chrm_dta2 <- merge_loc$data2

for (s in 2:22){
  file_name <- paste("11_10_0.1chrm",s,".RData", sep="")
  load(file_name)
  merge_chrm_dta1 <- rbind(merge_chrm_dta1, merge_loc$data1)
  merge_chrm_dta2 <- rbind(merge_chrm_dta2, merge_loc$data2)
}

merge_chrm <- list(merge_chrm_dta1 = merge_chrm_dta1,
                   merge_chrm_dta2 = merge_chrm_dta2)

save(merge_chrm, file = "11_10_1merge_chrm.RData")

theta = c(0, 0.5, 1)
theta_t = rep(theta, 3)
theta_n = rep(theta, c(3,3,3))
zeta_est = 0.1516
rho = zeta_est*theta_n + (1-zeta_est)*theta_t
