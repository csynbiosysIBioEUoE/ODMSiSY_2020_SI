
# ----------------------- DATA EXTRACTION/PROCESSING FUNCTION ----------------------- #

# Function to process the experimental data from the paper "Balancing a genetic toggle switch by real-time feedback 
# control and periodic forcing" into a list of objects to be passed to the STAN model for inference. Processing of the
# data can be done for a single experimental scheme (vectors) or for more than one experimental scheme at once (matrices)
# for a multiexperimental inference


# List with the names for the experiments presented in the paper "Balancing a genetic toggle switch by real-time feedback 
# control and periodic forcing". The function can take either a single element or else a vector from the following 
# experimental names

files<-c("Calibration_1", "Calibration_2","Calibration_3",
  "Calibration_4","Calibration_5","Calibration_6", "BangBang_1", "BangBang_2",
  "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_4", "DynStim_5", "DynStim_6", "DynStim_7", "DynStim_8", "DynStim_9", "DynStim_10", 
  "DynStim_11", "DynStim_12", "DynStim_13", "DynStim_14", "DynStimTooFast", "PI_1", "PI_2", "PI_3")



data_extraction_multiexperiment <- function (fileNamesVector){
  
  ######################## OBSERVABLES ########################
  
  mdp <- c() # Maximum data point for each file
  for(i in fileNamesVector){ # Loop that opens a file for each iteration to extract the file with maximum data points
    fileObs <- paste(i, "_Observables.csv", sep="")
    obs <- read.csv(file=fileObs, header=TRUE, sep=",")
    mdp = c(mdp, length(obs[,1]))
  }
  mrow <- max(mdp) # Maximum number of rows
  mcol <- length(fileNamesVector) # Maximum number of columns
  
  # Definition of the matrices with 0s on them
  samplingT <- matrix(data=0, nrow=mrow, ncol=mcol)
  GFPmean <- matrix(data=0, nrow=mrow, ncol=mcol)
  GFPstd <- matrix(data=0, nrow=mrow, ncol=mcol)
  RFPmean <- matrix(data=0, nrow=mrow, ncol=mcol)
  RFPstd <- matrix(data=0, nrow=mrow, ncol=mcol)
  
  for(i in 1:length(fileNamesVector)){ # Loop that opens a file for each iteration to extract the data
    fileObs <- paste(fileNamesVector[i], "_Observables.csv", sep="")
    obs <- read.csv(file=fileObs, header=TRUE, sep=",")
    for(j in 1:length(obs[,1])){
      samplingT[j,i] <- round(obs[j,1])
      GFPmean[j,i] <- obs[j,2]
      GFPstd[j,i] <- obs[j,3]
      RFPmean[j,i] <- obs[j,5]
      RFPstd[j,i] <- obs[j,6]
    }
  }
  mdp2 <- matrix(data=mdp, nrow=1, ncol = mcol)
  
  ######################## INPUTS ########################
  
  mdpI <- c() # Maximum data point for each file
  mtimes <- c() # Maximum time point for each file
  
  for(i in fileNamesVector){ # Loop that opens a file for each iteration to extract the file with maximum data points
    fileInp <- paste(i, "_Events_Inputs.csv", sep="")
    inp <- read.csv(file=fileInp, header=TRUE, sep=",")
    mdpI = c(mdpI, length(inp[,1]))
    mtimes = c(mtimes, round(inp[1,2]))
    
  }
  mrowI <- max(mdpI) # Maximum number of rows
  mcolI <- length(fileNamesVector) # Maximum number of columns
  mt <- max(mtimes) # Maximum number of time points
  
  mdpI2 <- matrix(data=mdpI, nrow=1, ncol = mcol)
  
  evnT <- matrix(data=0, nrow=mrowI+1, ncol=mcolI) # Event chance time points
  u_IPTG <- matrix(data=0, nrow=mrowI, ncol=mcolI) # External inducer concentrations
  u_aTc <- matrix(data=0, nrow=mrowI, ncol=mcolI)
  preI <- matrix(data=0, nrow=1, ncol=mcolI) # Inducer initial values
  preA <- matrix(data=0, nrow=1, ncol=mcolI)
  time <- matrix(data=0, nrow=mt+1, ncol=mcolI) # Time points vector
  inps <- matrix(data=0, nrow=mrowI*2, ncol=mcolI) # Inputs vector
  
  ltimes <- c() # Length of time series
  
  for(i in 1:length(fileNamesVector)){ # Loop that opens a file for each iteration to extract the data
    fileInp <- paste(fileNamesVector[i], "_Events_Inputs.csv", sep="")
    inp <- read.csv(file=fileInp, header=TRUE, sep=",")
    tempT <- seq(1e-9, round(inp[1,2]), length=round(inp[1,2])+1)
    ltimes = c(ltimes, length(tempT))
    ind = 1
    for(j in 1:length(inp[,1])){
      evnT[j,i] <- round(inp[j,1])
      u_IPTG[j,i] <- inp[j,5]
      u_aTc[j,i] <- inp[j,6]
      inps[ind,i] <- inp[j,5]
      inps[ind+1,i] <- inp[j,6]
      ind = ind+2
    }
    
    for(l in 1:length(tempT)){
      
      time[l,i] = tempT[l]
    }
    evnT[(length(inp[,1])+1),i] = round(inp[1,2])
    preI[,i] <- inp[1,3]
    preA[,i] <- inp[1,4]
    
  }
  ltimes2 <- matrix(data=ltimes, nrow=1, ncol = length(ltimes))
  
  toni <- seq(1e-9, 24*60) # 24h incuvation time vector
  
  data_multi <<- list (
    
    elm = mrowI, # Maximum length of the rows of the matrices except for time and evnT and pres
    tml = trunc(mt+1), # Maximum length of the rows for the time matrix
    ts = time+(1e-9), # Experimental time vector
    tsl = ltimes2, # length of time series per event
    tsmax = round(mtimes), # maximum time per event
    preIPTG = preI, # Inducer initial values
    preaTc = preA,
    IPTG = u_IPTG, # Inducer values
    aTc = u_aTc,
    Nsp = (mdpI2+1), # length(evnT) or number of event switching times
    inputs = inps+1e-7, # Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
    evnT = round(evnT), # event switching times
    m = mcol, # Number of time series
    stsl = mdp2, # Number of elements at each time series
    stslm = mrow, # Maximum number of time points
    sts = round(samplingT), # sampling time points
    GFPmean = GFPmean, # Means and standard deviations for the state variables LacI-RFP and TetR-RFP
    RFPmean = RFPmean,
    GFPstd = GFPstd,
    RFPstd = RFPstd,
    toni = toni, # Time points and time vector for the 24h steady state incubation 
    tonil = length(toni)
  )
}

