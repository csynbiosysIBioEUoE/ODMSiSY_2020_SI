############################################ GAUSSIAN MIXTURE M1 ############################################

# Function to compute the Gaussian Mixture fit to the parameter draws from stan to be used for the computation of the evidence
# using the Laplace approximation. As input takes the name of the tag that the inference result from stan has. 
# This script only works for Model 1 (Lugagne et.al.)


GMM1 <- function(experiment){
  
  ############## Extraction of all posterior samples from MCMC (STAN) results
  
  fitName <- paste("fit_", experiment, "_Model1.stan.rds", sep="") 
  x <- as.array(readRDS(fitName))
  
  k_in_IPTG <- c(x[,1,17], x[,2,17], x[,3,17], x[,4,17])
  k_out_IPTG <- c(x[,1,18], x[,2,18], x[,3,18], x[,4,18])
  k_in_aTc <- c(x[,1,19], x[,2,19], x[,3,19], x[,4,19])
  k_out_aTc <- c(x[,1,20], x[,2,20], x[,3,20], x[,4,20])
  k_L_pm0 <- c(x[,1,21], x[,2,21], x[,3,21], x[,4,21])
  k_L_pm <- c(x[,1,22], x[,2,22], x[,3,22], x[,4,22])
  theta_T <- c(x[,1,23], x[,2,23], x[,3,23], x[,4,23])
  theta_aTc <- c(x[,1,24], x[,2,24], x[,3,24], x[,4,24])
  n_aTc <- c(x[,1,25], x[,2,25], x[,3,25], x[,4,25])
  n_T <- c(x[,1,26], x[,2,26], x[,3,26], x[,4,26])
  k_T_pm0 <- c(x[,1,27], x[,2,27], x[,3,27], x[,4,27])
  k_T_pm <- c(x[,1,28], x[,2,28], x[,3,28], x[,4,28])
  theta_L <- c(x[,1,29], x[,2,29], x[,3,29], x[,4,29])
  theta_IPTG <- c(x[,1,30], x[,2,30], x[,3,30], x[,4,30])
  n_IPTG <- c(x[,1,31], x[,2,31], x[,3,31], x[,4,31])
  n_L <- c(x[,1,32], x[,2,32], x[,3,32], x[,4,32])
  
  ############## Introduce MCMC results into a dataframe to work with 
  
  s <- matrix(data = 0, nrow = length(k_IPTG), ncol = 16)
  
  s[,1] <- k_in_IPTG
  s[,2] <- k_out_IPTG
  s[,3] <- k_in_aTc
  s[,4] <- k_out_aTc
  s[,5] <- k_L_pm0
  s[,6] <- k_L_pm
  s[,7] <- theta_T
  s[,8] <- theta_aTc
  s[,9] <- n_aTc
  s[,10] <-n_T
  s[,11] <-k_T_pm0
  s[,12] <-k_T_pm
  s[,13] <-theta_L
  s[,14] <-theta_IPTG
  s[,15] <-n_IPTG
  s[,16] <-n_L
  
  y <- data.frame(s)
  
  ############### Gaussian Mixtures with MClust
  library(mclust)
  
  mvpdfPost <- densityMclust(y, G=1:30)
  
  # Save MClust results
  mcres <- paste("mvpdfPost_", experiment, "_Model1.rds", sep = "")
  saveRDS(mvpdfPost, mcres)
  
}



