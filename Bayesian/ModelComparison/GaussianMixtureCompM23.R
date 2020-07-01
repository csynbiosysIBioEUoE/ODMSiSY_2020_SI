############################################ GAUSSIAN MIXTURE M1 ############################################

# Function to compute the Gaussian Mixture fit to the parameter draws from stan to be used for the computation of the evidence
# using the Laplace approximation. As input takes the name of the tag that the inference result from stan has. 
# This script only works for Model 2 (intermediate) and 3 (new)


GMM23 <- function(experiment){
  
  ############## Extraction of all posterior samples from MCMC (STAN) results
  
  fitName <- paste("fit_", experiment, ".rds", sep="") 
  x <- as.array(readRDS(fitName))
  
  k_IPTG <- c(x[,1,15], x[,2,15], x[,3,15], x[,4,15])
  k_aTc <- c(x[,1,16], x[,2,16], x[,3,16], x[,4,16])
  k_L_pm0 <- c(x[,1,17], x[,2,17], x[,3,17], x[,4,17])
  k_L_pm <- c(x[,1,18], x[,2,18], x[,3,18], x[,4,18])
  theta_T <- c(x[,1,19], x[,2,19], x[,3,19], x[,4,19])
  theta_aTc <- c(x[,1,20], x[,2,20], x[,3,20], x[,4,20])
  n_aTc <- c(x[,1,21], x[,2,21], x[,3,21], x[,4,21])
  n_T <- c(x[,1,22], x[,2,22], x[,3,22], x[,4,22])
  k_T_pm0 <- c(x[,1,23], x[,2,23], x[,3,23], x[,4,23])
  k_T_pm <- c(x[,1,24], x[,2,24], x[,3,24], x[,4,24])
  theta_L <- c(x[,1,25], x[,2,25], x[,3,25], x[,4,25])
  theta_IPTG <- c(x[,1,26], x[,2,26], x[,3,26], x[,4,26])
  n_IPTG <- c(x[,1,27], x[,2,27], x[,3,27], x[,4,27])
  n_L <- c(x[,1,28], x[,2,28], x[,3,28], x[,4,28])
  
  ############## Introduce MCMC results into a dataframe to work with 
  
  s <- matrix(data = 0, nrow = length(k_IPTG), ncol = 14)
  
  s[,1] <- k_IPTG
  s[,2] <- k_aTc
  s[,3] <- k_L_pm0
  s[,4] <- k_L_pm
  s[,5] <- theta_T
  s[,6] <- theta_aTc
  s[,7] <-n_aTc
  s[,8] <-n_T
  s[,9] <-k_T_pm0
  s[,10] <-k_T_pm
  s[,11] <-theta_L
  s[,12] <-theta_IPTG
  s[,13] <-n_IPTG
  s[,14] <-n_L
  
  y <- data.frame(s)
  
  ############### Gaussian Mixtures with MClust
  library(mclust)
  
  mvpdfPost <- densityMclust(y, G=1:30)
  
  # Save MClust results
  mcres <- paste("mvpdfPost_", experiment, "_Model1.rds", sep = "")
  saveRDS(mvpdfPost, mcres)
  
}
