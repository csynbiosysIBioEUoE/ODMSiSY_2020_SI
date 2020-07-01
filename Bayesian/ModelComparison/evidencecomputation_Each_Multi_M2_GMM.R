# This script computes the Evidence that the data in experiment j was generated using model i 
# Considering the Laplace approximation, the evidence is given by:
# Evidence ~ BestFitlikelihood * Occam factor
# where, 
# BestFitLikelihood = P(Dj|theta_MAP,Mi)
# Occam factor = P(theta_MAP|Dj,Mi)*det(A/2pi)^(-0.5) and A is the Hessian of the posterior evaluated at theta_MAP

# We therefore need to identify theta_MAP and the Hessian 
# theta_MAP is the maximum a posteriori estimate, and can be obtained as the solution of an optimisation problem
# searching for the maximum of the posterior distribution (which we have approximated with gaussian mixtures)

FileNames = c("Calibration_4","Calibration_5","Calibration_6","DynStim_1","DynStim_2", "DynStim_3","DynStim_8", "DynStim_9","DynStim_11", "DynStim_14")
l2_single <- c()
hdet <- c()

GM_data <- readRDS(file=paste("GaussianMixtures/mvpdfPost_ALL","_Model2.rds",sep=""))
MU <- GM_data$parameters$mean
E <- GM_data$parameters$variance$sigma
w <- GM_data$parameters$pro


LogGaussMix <- function(x, MU, E, w){
  comp <- length(w)
  FX <- c()
  
  for(i in 1:comp){
    
    fx <- w[i]*(1/sqrt(det(E[,,i])*(2*pi)^2))*exp((-1/2)*t(x-MU[,i])%*%solve(E[,,i])%*%(x-MU[,i]))
    FX <- c(FX, fx)
  }
  
  return(log(sum(FX))) # Added log, because we would need to compute the Hessian 
}

hessdet <- c()
Post_MAP_vect <-c()
theta_MAP_matrix <- matrix(data = 0, nrow = length(w),ncol = 16)
for(i in 1:length(w)){
  
  fit3 <- optim(MU[,i], LogGaussMix, control = list(fnscale = -1), method = "BFGS", hessian = TRUE, MU = MU, E =E, w = w)
  theta_MAP_matrix[i,] <- fit3$par
  hessdet <- c(hessdet,det(fit3$hessian))
  Post_MAP_vect <- c(Post_MAP_vect,fit3$value)
}

theta_MAP <- theta_MAP_matrix[which.max(Post_MAP_vect),]
Post_MAP <- Post_MAP_vect[which.max(Post_MAP_vect)]


for(e in 1:length(FileNames)){ # iterate over the experiments
   # Load the results of the gaussian mixture approximation
    FL <- FileNames[e]
    
    # Now we need to compute the best fit likelihood
    
    # Function to simulate the system for a given experimental profile (expert) and vector of parameters (param)
    Simul <- function(exper, param){
      library(rstan)
      expose_stan_functions("Simulations/Model2_Function_FakeData.stan")
      
      # Load all required data from the csv files for each experimental profile
      fileName <- exper
      fileInputs <- paste("/Users/lucia/Documents/Projects/BayGame/ExperimentalDatacsv/",FileNames[e], "_Events_Inputs.csv", sep="")
      
      # Extract input data and stored into global variables
      inputs <- read.csv(file=fileInputs, header=TRUE, sep=",")
      evnT <- c(round(inputs[,1]),round(inputs[1,2]))
      u_IPTG<- inputs[,5]
      u_aTc<-inputs[,6]
      preI<-inputs[1,3]
      preA<-inputs[1,4]
      inp<-c()
      # Inputs
      for (x in 1:length(u_IPTG)){
        inp <- c(inp, u_IPTG[x], u_aTc[x])
      }
      inp <- inp+1e-7
      
      et <- round(inputs[1,2])
      
      # Time series for ON incubation 
      toni <- seq(1e-9, 24*60)
      
      # Time series
      ts <- seq(1e-9, et, length=(et+1))
      
      fileObservables <- paste(fileName, "_Observables.csv", sep="")
      # Extract observables data and stored into global variables
      observables <- read.csv(file=fileObservables, header=TRUE, sep=",")
      samplingT <- round(observables[,1])
      GFPmean <- observables[,2]
      GFPstd <- observables[,3]
      RFPmean <- observables[,5]
      RFPstd <- observables[,6]
      
      ivss <- c(preI, preA, RFPmean[1], GFPmean[2])
      pre <- c(preI, preA)
      
      # Solve ODEs system
      # 1 = IPTG, 2 = aTc, 3 = RFP, 4 = GFP
      fd2 <- solve_coupled_ode(ts, param, 0, 0, evnT, inp, toni,ivss, pre)
      
      res <- matrix(data=NaN, nrow = length(GFPmean), ncol = 2)
      w = 1
      for(y in seq(1, et+1, 5)){
        res[w,1] = fd2[y,3] 
        res[w,2] = fd2[y,4]
        w = w+1
      }
      
      return(res)
      
    }
    
    simRes <- Simul(FL,type.convert(theta_MAP))
    
    BestFitLikelihood <- function(simRes,exper){
      fileName <- exper
      fileInputs <- paste("/Users/lucia/Documents/Projects/BayGame/ExperimentalDatacsv/",fileName, "_Events_Inputs.csv", sep="")
      fileObservables <- paste(fileName, "_Observables.csv", sep="")
      # Extract observables data and stored into global variables
      observables <- read.csv(file=fileObservables, header=TRUE, sep=",")
      samplingT <- round(observables[,1])
      GFPmean <- observables[,2]
      GFPstd <- observables[,3]
      RFPmean <- observables[,5]
      RFPstd <- observables[,6]
      
    
      LlkRFP <- -sum(log(RFPstd))-length(RFPstd)/2*log(2*pi)-0.5*sum((simRes[,1]-RFPmean)^2/RFPstd^2) # Loglikelihood RFP
      LlkGFP <- -sum(log(GFPstd))-length(GFPstd)/2*log(2*pi)-0.5*sum((simRes[,2]-GFPmean)^2/GFPstd^2) # Loglikelihood GFP
      Llk <- LlkRFP+LlkGFP
      
      return(Llk)
      
    }
    
    Llk <-BestFitLikelihood(simRes, FL)
    
    # Compute the logarithm of the prior term 
    pr<-numeric(length(theta_MAP))
    LogPriorTerm <-function(theta_MAP){
      # from each marginal prior, compute the p at the MAP
      pr[1] = dlnorm(theta_MAP[1], meanlog = -3.2188758248682, sdlog = 1.15129254649702) 
      pr[2] = dlnorm(theta_MAP[2], meanlog = -3.2188758248682, sdlog = 1.15129254649702)
      pr[3] = dlnorm(theta_MAP[3], meanlog = -2.30258509299405, sdlog = 1.15129254649702)
      pr[4] = dlnorm(theta_MAP[4], meanlog = -2.30258509299405, sdlog = 1.15129254649702)
      pr[5] = dlnorm(theta_MAP[5], meanlog = -3.50655789731998, sdlog = 1.15129254649702)
      pr[6] = dlnorm(theta_MAP[6],meanlog = 2.30258509299405,sdlog = 1.15129254649702)
      pr[7] = dlnorm(theta_MAP[7],meanlog = 3.40119738166216, sdlog = 1.15129254649702)
      pr[8] = dlnorm(theta_MAP[8],meanlog=2.30258509299405,sdlog = 1.15129254649702)
      pr[9] = dnorm(theta_MAP[9],mean = 2.5,sd = 1.25);
      pr[10] = dnorm(theta_MAP[10],mean = 2.5,sd = 1.25);
      pr[11] = dlnorm(theta_MAP[11],meanlog=2.30258509299405,sdlog = 1.15129254649702)
      pr[12] = dlnorm(theta_MAP[12],meanlog=2.22044604925031e-16,sdlog = 1.15129254649702)
      pr[13] = dlnorm(theta_MAP[13],meanlog=3.40119738166216,sdlog = 1.15129254649702)
      pr[14] = dlnorm(theta_MAP[14],meanlog=-2.30258509299405,sdlog = 1.15129254649702)
      pr[15] = dnorm(theta_MAP[15],mean = 2.5,sd = 1.25);
      pr[16] = dnorm(theta_MAP[16],mean = 2.5,sd = 1.25);
      
      LPT <- sum(log(pr))
      return(LPT)
    }
    
    Prior_MAP <- LogPriorTerm(theta_MAP)

# The Hessian is the inverse of the covariance matrix: 
# H = C^-1
# det(H) = det(C^-1) = 1/det(C)
# The covariance matrix can be computed from the samples of the mcmc (as we did in the past). 
# Would it be better to compute the Hessian using the approximation of the posterior?


    ### Try to compute Hessian considering analytical solution for Gaussian Mixtures
    ## Using reference "Mode-finding for mixtures of Gaussian distributions"
    GaussMix <- function(x, MU, E, w){
      comp <- length(w)
      FX <- c()
      
      for(i in 1:comp){
        
        fx <- w[i]*(1/sqrt(det(E[,,i])*(2*pi)^2))*exp((-1/2)*t(x-MU[,i])%*%solve(E[,,i])%*%(x-MU[,i]))
        FX <- c(FX, fx)
      }
      
      return(sum(FX)) 
    }
    
    p<-GaussMix(theta_MAP,MU,E,w)
    
    hx <- array(data = NA, dim = c(length(theta_MAP),length(theta_MAP),length(w)))
    HessianComp <- function(x, MU, E, w){
      comp <- length(w)
      HX <- matrix(data = 0,nrow = length(theta_MAP),ncol = length(theta_MAP))
      for(i in 1:comp){
        
        
        NormComponent <-as.vector(w[i]*(1/sqrt(det(E[,,i])*(2*pi)^2))*exp((-1/2)*t(x-MU[,i])%*%solve(E[,,i])%*%(x-MU[,i])))
        m <- (solve(E[,,i])%*%((MU[,i]-x)%*%t(MU[,i]-x)-E[,,i])%*%solve(E[,,i]))
        hx[,,i] <- NormComponent*m
        
        HX <- HX+ hx[,,i]
        print(HX)
      }
      return(HX)
    }
    
    H<- HessianComp(theta_MAP,MU,E,w)
    
    GradientComp <- function(x, MU, E, w){
      comp <- length(w)
      GX <- c(0)
  
      for(i in 1:comp){
        
        
        NormComponent <- as.vector(w[i]*(1/sqrt(det(E[,,i])*(2*pi)^2))*exp((-1/2)*t(x-MU[,i])%*%solve(E[,,i])%*%(x-MU[,i])))
        gx <- NormComponent*solve(E[,,i])%*%(MU[,i]-x)
        GX <- GX+ gx
      }
      
      return(GX) 
    }
    
    g<- GradientComp(theta_MAP,MU,E,w)
    
    MinusHessianLogP <- function(p,g,H){
      mhlogp <- 1/p*g%*%t(g)-1/p*H
      return(mhlogp)
      
    }
    
    Hlooked <-MinusHessianLogP(p,g,H)
    hdet <- c(hdet,det(Hlooked))
    
    LogEvidence <- function(Llk,hess,Prior_MAP,theta_MAP){
      logEv <- Llk+Prior_MAP+length(theta_MAP)/2*log(2*pi)-0.5*log(det(hess)) 
      # log(|A/(2*pi)|^-0.5) = log(1/sqrt(|A/(2*pi)|)) = N/2*log(2*pi)-1/2*log(|A|), where N is the dimension of the parameter space. We considered that |kA| = k^N*|A|
      return(logEv)
    }
    
    l<-LogEvidence(Llk,Hlooked,Prior_MAP,theta_MAP)
    
    l2_single<- c(l2_single,l)
}

 # Save results as a CSV file
 tit <- paste("lnEvidenceModel2_GMM_Each_Multi.csv", sep = "")
 
 write.table(t(l2_single), file = tit, col.names =c("Calibration_4","Calibration_5","Calibration_6","DynStim_1","DynStim_2", "DynStim_3","DynStim_8", "DynStim_9","DynStim_11", "DynStim_14"), 
             row.names = "lnEvidenceModel2_GMM", sep=",")
 
 # tit <- paste("DetHessianModel2_GMM.csv", sep = "")
 # 
 # write.table(t(hdet), file = tit, col.names =FileNames, 
 #             row.names = "DetHessianModel2_GMM", sep=",")


