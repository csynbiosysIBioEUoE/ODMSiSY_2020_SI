################################# Generate Data Files for Optimum Experiments #################################

# Function used to generate the necessary CSV Data Files (Events_Inputs, Inputs and Observables) in order to generate plots
# and simulate the three different models to a specific set of inputs using the functions enclosed in Predictions&Analysis.
# The CSV files generated have the same structure as the ones extracted from Lugagne et.al. 
# As arguments it takes fileName (name to be taged to the different CSV files as a marker) and iputs (a vector with the input
# values ordered as IPTG1, aTc1, IPTG2, aTc2, ..., IPTGn, aTcn) and this can have any desired length (from 1 step experiment
# to n steps experiment) extracted from the optimisation results.



OptimFileGener <- function(fileName = "Optim", inputs = c(0.8100001, 6.0000001)){

  
  ########################### Generate Observables file
  baseFile <- read.csv("Optim_Observables.csv")
  
  write.csv(baseFile, file = paste(fileName, "_Observables.csv", sep = ""), row.names = FALSE)
  
  
  ########################### Generate Events_Inputs file
  
  FinalTime <- rep(24*60, length(inputs)/2) # Time vetor of the experiment (24h)
  Switchingtimes <- seq(0, 24*60, length = ((length(inputs)/2)+1)) # Switching times where each step has the same duration
  Switchingtimes <- Switchingtimes[1:length(FinalTime)] 
  IPTGpre <- rep(1, length(inputs)/2) # Input values for the 24h incubation
  aTcpre <- rep(0, length(inputs)/2)
  IPTG <- c() # Input values for the experiment
  aTc <- c()
  for(j in seq(1, length(inputs), 2)){
    IPTG <- c(IPTG, inputs[j])
    aTc <- c(aTc, inputs[j+1])
  }
  
  # Write results in a CSV file
  allTog <- cbind(Switchingtimes, FinalTime, IPTGpre, aTcpre, IPTG, aTc)
  
  write.csv(allTog, file = paste(fileName, "_Events_Inputs.csv", sep = ""), row.names = FALSE)
  
  
  ########################### Generate Inputs file
  # Same as the previous file, but instead of ordered by events ordered by time points
  
  time <- seq(0, 24*60, length = (24*60)+1)
  IPTGpre <- rep(1, (24*60)+1)
  aTcpre <- rep(0, (24*60)+1)
  IPTG <- c()
  aTc <- c()
  j <- 1
  Switchingtimes <- seq(0, 24*60, length = ((length(inputs)/2)+1))
  for(i in seq(1, length(inputs), 2)){
    if(i == 1){
      ip <- rep(inputs[i], Switchingtimes[j+1]-Switchingtimes[j]+1)
      at <- rep(inputs[i+1], Switchingtimes[j+1]-Switchingtimes[j]+1)
    } else {
      ip <- rep(inputs[i], Switchingtimes[j+1]-Switchingtimes[j])
      at <- rep(inputs[i+1], Switchingtimes[j+1]-Switchingtimes[j])
    }
    
    IPTG <- c(IPTG, ip)
    aTc <- c(aTc, at)
    
    j <- j+1
  }
  
  allTog2 <- cbind(time, IPTGpre, aTcpre, IPTG, aTc)
  
  write.csv(allTog2, file = paste(fileName, "_Inputs.csv", sep = ""), row.names = FALSE)

}


######### Generate Posterior Predictive Distributions files if required

#### MODEL 1 ####

# Path needs to be modified according to the user
source('D:/PhD/GitHub/FOSBE2019_Paper/Predictions&Analysis/PostPredCheckSimulM1.R')
# Inputs to the function need to be modified according to the user results
PPCSimul1(fileName, "ALL_Model1.stan")

#### MODEL 2 ####

# Path needs to be modified according to the user
source('D:/PhD/GitHub/FOSBE2019_Paper/Predictions&Analysis/PostPredCheckSimulM2.R')
# Inputs to the function need to be modified according to the user results
PPCSimul2(fileName, "ALL_Model2.stan")

#### MODEL 3 ####

# Path needs to be modified according to the user
source('D:/PhD/GitHub/FOSBE2019_Paper/Predictions&Analysis/PostPredCheckSimulM3.R')
# Inputs to the function need to be modified according to the user results
PPCSimul3(fileName, "ALL_Model3.stan")





