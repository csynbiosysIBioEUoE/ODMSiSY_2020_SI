
# -------------------------- INFERENCE RUN SCRIPT -------------------------- #

# Script to perform multi-experimental or single experiment inference without the optimisation step contained in masterRunOptim.R. 
# As inputs takes the name of the experimental profile/s (experiments), the stan model name (ended with .stan) to use (modelName) and
# a character string (infP) that can be either "single" for a single experimental fit (which can be in series) or else "multi" for a 
# multiexperimental fit (all data profiles at once)

runAllExp2 <- function (files, modelName, infP = "single"){
    
  if(infP=="multi"){
    # Data Extraction (script MultiExtractExp with the function needs to be run first)
    data_extraction_multiexperiment(files)
    
    # Inference
    
    inter <- paste("fit_ALL_", modelName, ".rds", sep="")
    fit <- stan(file=modelName, data = data_multi, iter = 3000, warmup = 1000, chains = 4, init = p, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
    # Save fit result as rds file
    fit.dso <- new("cxxdso")
    saveRDS(fit, file = inter)
    
  } else if (infP=="single"){
    
      for(x in files){
        
        data_extraction_multiexperiment(x)
        inter <- paste("fit_", x,"_",modelName, ".rds", sep="")
        fit <- stan(file=modelName, data = data_multi, iter = 3000, warmup = 1000, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
        # Save fit result as rds file
        fit.dso <- new("cxxdso")
        saveRDS(fit, file = inter)
        
        mes <- paste("Fit for ", x , " has finished!", sep = "")
        print(mes)
      }
    
  } else {
    print("Please, select either single or else multi")
  }
  
}