
################################# CONFIDENCE INTERVAL PLOTS #################################

# Function to plot the 95% confidence intervals plots for the inference results on the Toggle sqwitch Model. First argument 
# is the name of the experimental design to plot, and second argument aditional information that will be included 
# in the title of the plot and the file name. The third argument 
# is the name of the rstan object (result of inference) to be used for the analysis. All files need to be in the same directory.
# ODE_Model2_Function.stan is a stan function file necessary to simulate the model for each experimental profile
# This script only works for Model 2. 

ConfidenceInterval2 <- function(experiment, aditInfo = "", stanRes){
  
  expose_stan_functions("ODE_Model2_Function.stan")
  
  past <- as.array(readRDS(paste("fit_", stanRes, ".rds", sep="")))
  
  fileName <- experiment
  
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
    
    q <- c(past[x,1,15],past[x,1,16],past[x,1,17],past[x,1,18],past[x,1,19],past[x,1,20],past[x,1,21],
           past[x,1,22],past[x,1,23],past[x,1,24], past[x,1,25], past[x,1,26], past[x,1,27], 
           past[x,1,28])
    fd2 <- solve_coupled_ode(ts, q, 0, 0, evnT, inp, toni,ivss, pre)
    chr[,x] <- fd2[,3]
    chg[,x] <- fd2[,4]
  }  
  
  
  for(x in 1:2000){
    q <- c(past[x,2,15],past[x,2,16],past[x,2,17],past[x,2,18],past[x,2,19],past[x,2,20],past[x,2,21],
           past[x,2,22],past[x,2,23],past[x,2,24], past[x,2,25], past[x,2,26], past[x,2,27],
           past[x,2,28])
    fd2 <- solve_coupled_ode(ts, q, 0, 0, evnT, inp, toni,ivss, pre)
    chr[,(x+2000)] <- fd2[,3]
    chg[,(x+2000)] <- fd2[,4]
    
  }
  
  for(x in 1:2000){
    q <- c(past[x,3,15],past[x,3,16],past[x,3,17],past[x,3,18],past[x,3,19],past[x,3,20],past[x,3,21],
           past[x,3,22],past[x,3,23],past[x,3,24], past[x,3,25], past[x,3,26], past[x,3,27],
           past[x,3,28])
    fd2 <- solve_coupled_ode(ts, q, 0, 0, evnT, inp, toni,ivss, pre)
    chr[,(x+4000)] <- fd2[,3]
    chg[,(x+4000)] <- fd2[,4]
    
  }
  
  for(x in 1:2000){
    q <- c(past[x,4,15],past[x,4,16],past[x,4,17],past[x,4,18],past[x,4,19],past[x,4,20],past[x,4,21],
           past[x,4,22],past[x,4,23],past[x,4,24], past[x,4,25], past[x,4,26], past[x,4,27],
           past[x,4,28])
    fd2 <- solve_coupled_ode(ts, q, 0, 0, evnT, inp, toni,ivss, pre)
    chr[,(x+6000)] <- fd2[,3]
    chg[,(x+6000)] <- fd2[,4]
    
  }
  
  # Extraction of the median
  
  chrMED <- c()
  chgMED <- c()
  
  for(q in 1:length(chr[,1])){
    M1 <- quantile(chr[q,], probs = 0.5)
    M2 <- quantile(chg[q,], probs = 0.5)
    chrMED <- c(chrMED, M1[[1]])
    chgMED <- c(chgMED, M2[[1]])
  }
  
  ## 95% confidence regions
  chr25 <- c()
  chg25 <- c()
  
  for(q in 1:length(chr[,1])){
    M1 <- quantile(chr[q,], probs = 0.025)
    M2 <- quantile(chg[q,], probs = 0.025)
    chr25 <- c(chr25, M1[[1]])
    chg25 <- c(chg25, M2[[1]])
  }
  
  chr975 <- c()
  chg975 <- c()
  
  for(q in 1:length(chr[,1])){
    M1 <- quantile(chr[q,], probs = 0.975)
    M2 <- quantile(chg[q,], probs = 0.975)
    chr975 <- c(chr975, M1[[1]])
    chg975 <- c(chg975, M2[[1]])
  }
  
  # Save results into CSV files
  
  resPlotR <- matrix(data = 0, nrow = length(chrMED), ncol = 3)
  resPlotG <- matrix(data = 0, nrow = length(chgMED), ncol = 3)
  
  resPlotR[,1] <- chr25
  resPlotR[,2] <- chrMED
  resPlotR[,3] <- chr975
  colnames(resPlotR) <- c("2.5CI", "Median", "97.5CI")
  colnames(resPlotG) <- c("2.5CI", "Median", "97.5CI")
  
  resPlotG[,1] <- chg25
  resPlotG[,2] <- chgMED
  resPlotG[,3] <- chg975
  
  write.table(resPlotR, file = paste("Model2_ODEs_ConfInter_Params", stanRes, "_Exper",experiment,"_RFP.csv", sep = ""), col.names = TRUE, sep = ",")
  write.table(resPlotG, file = paste("Model2_ODEs_ConfInter_Params", stanRes, "_Exper",experiment,"_GFP.csv", sep = ""), col.names = TRUE, sep = ",")
  
  #### Plot GFP
  
  a = data.frame(time = ts, x = chgMED)
  b = data.frame(time = ts, x = chg975)
  c = data.frame(time = ts, x = chg25)
  
  min_a <- pmin(a, b, c)
  max_a <- pmax(a, b, c)
  
  
  png(paste("CI_Mm_GFP_Model2_", fileName, "_Fit", aditInfo,".png", sep=""), width = 6, height = 4, units = 'in', res = 300)
  plot(samplingT, GFPmean, col = "black", pch=20, cex = 0.1, xlab="time (min)", ylab="GFP", 
       main=paste("Results Inference for GFP ", fileName, ", ODE_Model, ",aditInfo, sep=""), 
       ylim = c(min(GFPmean)-GFPstd[which.min(GFPmean)]-10, max(GFPmean)+GFPstd[which.max(GFPmean)]+10))
  arrows(samplingT, GFPmean-GFPstd, samplingT, GFPmean+GFPstd, length=0.03, angle=90, code=3)
  
  # Uncoment next two lines if legend needs to be inlcuded. X and Y values for its location in the plot need to be specified
  # legend(500, 2000, legend=c("Confidence Interval", "HDI"),
  # col=c("green", "blue"), lty=1, cex=0.5) # for GFP
  
  polygon(c(c$time, rev(c$time)), c(max_a$x ,rev(min_a$x)), col = rgb(0.1, 1, 0.1,0.6), border = "green")
  lines(a$time, a$x, type = 'l',col="darkgreen")
  
  dev.off()
  
  
  #### Plot RFP
  
  a = data.frame(time = ts, x = chrMED)
  b = data.frame(time = ts, x = chr975)
  c = data.frame(time = ts, x = chr25)
  
  min_a <- pmin(a, b, c)
  max_a <- pmax(a, b, c)
  
  png(paste("CI_Mm_RFP_Model2_", fileName, "_Fit", aditInfo,".png", sep=""), width = 6, height = 4, units = 'in', res = 300)
  plot(samplingT, RFPmean, col = "black", pch=20, cex = 0.1, xlab="time (min)", ylab="RFP", 
       main=paste("Results Inference for RFP ", fileName, ", ODE_Model, ",aditInfo , sep=""), 
       ylim = c(min(RFPmean)-RFPstd[which.min(RFPmean)]-10, max(RFPmean)+RFPstd[which.max(RFPmean)]+10))
  arrows(samplingT, RFPmean-RFPstd, samplingT, RFPmean+RFPstd, length=0.03, angle=90, code=3)
  
  # Uncoment next two lines if legend needs to be inlcuded. X and Y values for its location in the plot need to be specified
  # legend(900, 1500, legend=c("Confidence Interval", "HDI"),
  #        col=c("red", "blue"), lty=1, cex=0.5) # for RFP
  
  polygon(c(c$time, rev(c$time)), c(max_a$x ,rev(min_a$x)), col = rgb(1, 0, 0,0.6), border = "red")
  lines(a$time, a$x, type = 'l',col="brown4")
  
  dev.off()
  
}
