
# --------------------- POSTERIOR PREDICTIVE CHECKS FUNCTION --------------------- #

# Function to simulate the ODE system for a determined experimental profile (input named experiment) using all the samples
# obtained by a determined inference result in the form of a stanfit object (input named experiment2). If both inputs are the
# same, predictions for all experimental results using all the posteriors obtained will be performed. This script only works for
# Model 1. 

expert <- c("Calibration_4","Calibration_5","Calibration_6",
            "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_8","DynStim_9", 
            "DynStim_11", "DynStim_14", "ALL_Long_Model3.stan")
expert2 <- c("Calibration_4","Calibration_5","Calibration_6",
             "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_8","DynStim_9", 
             "DynStim_11", "DynStim_14", "ALL_Long_Model3.stan")

PPCSimul1 <- function(experiment, experiment2){
  
  expose_stan_functions("OED_Model1_Function.stan")
  
  # Iterations over stanfit objects
  for(i in 1:length(experiment2)){
    
    stanRes <- paste("fit_", experiment2[i], ".rds", sep = "")
    past <- as.array(readRDS(stanRes))
    
    # Iterations over experimental profiles
    for(j in 1:length(experiment)){
      
      # Load all required data from the csv files for each experimental profile
      
      fileName <- experiment[j]
      
      fileInputs <- paste(fileName, "_Events_Inputs.csv", sep="")
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
      
      # Empty matrices to include the simulation results
      chr <- matrix(data = 0, ncol = 8000, nrow = et+1)
      chg <- matrix(data = 0, ncol = 8000, nrow = et+1)
      
      # Simulation of the ODE system for each sample from the stanfit object
      
      for(x in 1:2000) {
        
        q <- c(past[x,1,17],past[x,1,18],past[x,1,19],past[x,1,20],past[x,1,21],
               past[x,1,22],past[x,1,23],past[x,1,24], past[x,1,25], past[x,1,26],
               past[x,1,27],past[x,1,28],past[x,1,29],past[x,1,30],past[x,1,31],past[x,1,32])
        fd2 <- solve_coupled_ode(ts, q, 0, 0, evnT, inp, toni,ivss, pre)
        chr[,x] <- fd2[,3]
        chg[,x] <- fd2[,4]
      }  
      
      
      for(x in 1:2000){
        q <- c(past[x,1,17],past[x,1,18],past[x,1,19],past[x,1,20],past[x,1,21],
               past[x,1,22],past[x,1,23],past[x,1,24], past[x,1,25], past[x,1,26],
               past[x,1,27],past[x,1,28],past[x,1,29],past[x,1,30],past[x,1,31],past[x,1,32])
        fd2 <- solve_coupled_ode(ts, q, 0, 0, evnT, inp, toni,ivss, pre)
        chr[,(x+2000)] <- fd2[,3]
        chg[,(x+2000)] <- fd2[,4]
        
      }
      
      for(x in 1:2000){
        q <- c(past[x,1,17],past[x,1,18],past[x,1,19],past[x,1,20],past[x,1,21],
               past[x,1,22],past[x,1,23],past[x,1,24], past[x,1,25], past[x,1,26],
               past[x,1,27],past[x,1,28],past[x,1,29],past[x,1,30],past[x,1,31],past[x,1,32])
        fd2 <- solve_coupled_ode(ts, q, 0, 0, evnT, inp, toni,ivss, pre)
        chr[,(x+4000)] <- fd2[,3]
        chg[,(x+4000)] <- fd2[,4]
        
      }
      
      for(x in 1:2000){
        q <- c(past[x,1,17],past[x,1,18],past[x,1,19],past[x,1,20],past[x,1,21],
               past[x,1,22],past[x,1,23],past[x,1,24], past[x,1,25], past[x,1,26],
               past[x,1,27],past[x,1,28],past[x,1,29],past[x,1,30],past[x,1,31],past[x,1,32])
        fd2 <- solve_coupled_ode(ts, q, 0, 0, evnT, inp, toni,ivss, pre)
        chr[,(x+6000)] <- fd2[,3]
        chg[,(x+6000)] <- fd2[,4]
        
      }
      
      
      # Modify output to save only the sampling times for posterior comparisons with real data
      chrST <- matrix(data = 0, ncol = 8000, nrow = (et/5)+1)
      chgST <- matrix(data = 0, ncol = 8000, nrow = (et/5)+1)
      
      for(q in 1:8000){
        w = 1
        for(y in seq(1, et+1, 5)){
          chgST[w,q] = chg[y,q]
          chrST[w,q] = chr[y,q]
          w = w+1
        }
        
      }
      
      # Save results indicating in the name of the file from which stanfit object the samples come from and which experimental 
      # profile has been simulated.
      
      resnam1 <- paste("Model1_InferenceResults_", experiment2[i], "_Simulation_", experiment[j], "_RFP.csv", sep = "")
      resnam2 <- paste("Model1_InferenceResults_", experiment2[i], "_Simulation_", experiment[j], "_GFP.csv", sep = "")
      
      write.table(chrST, file = resnam1, sep=",")
      write.table(chgST, file = resnam2, sep=",")
      
    }
    
  }
  
}

