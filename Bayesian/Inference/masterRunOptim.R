
# -------------------------- INFERENCE RUN SCRIPT OPTIMISATION -------------------------- #

# Script to perform single or multi experiment inference with the calculation of the initial guesses for the four chains using the 
# stan function optimizing for the computation of a point estimate. As inputs takes the name of the experimental profile (experiments), 
# the stan model to use (modelName) and the threshold value for the optimizing function (OptimVal) to finalise the while loop. The number
# of parameters of the model needs to be passed into parNu (with default value of 14 for our case). The last input of the duntion is 
# character string (infP) that can be either "single" for a single experimental fit (which can be in series) or else "multi" for a 
# multiexperimental fit (all data profiles at once). Function can take as experiments input a string vector to perform inference
# in series for the experimental results indicated. 

files<-c(
  "Calibration_1", "Calibration_2","Calibration_3",
  "Calibration_4","Calibration_5","Calibration_6", "BangBang_1", "BangBang_2",
  "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_4", "DynStim_5", "DynStim_6", "DynStim_7", "DynStim_8", "DynStim_9", "DynStim_10", 
  "DynStim_11", "DynStim_12", "DynStim_13", "DynStim_14", "DynStimTooFast", "PI_1", "PI_2", "PI_3")

runExp2 <- function (experiments, modelName, OptimVal, infP = "single", parNu = 14){
  if(infP=="multi"){
    
    # Data Extraction (script MultiExtractExp with the function needs to be run first)
    data_extraction_multiexperiment(experiments)
    
    # Initial value optimisation
    m <- stan_model(modelName)
    s1 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
    s2 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
    s3 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
    s4 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
    
    # Refinment of the initial value guess
    while (s1$value < OptimVal || s2$value < OptimVal || s3$value < OptimVal || s4$value < OptimVal){
      if (s1$value < OptimVal){
        print("s1")
        s1 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
      }
      else if (s2$value < OptimVal){
        print("s2")
        s2 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
      }
      else if (s3$value < OptimVal){
        print("s3")
        s3 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
      }
      else if (s4$value < OptimVal){
        print("s4")
        s4 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
      }
      
    }
    
    # Formating initial values to be passed to the stan fit function as for different lists (one for each chain)
    p1 <- list()
    p2 <- list()
    p3 <- list()
    p4 <- list()
    
    q <- c(1:parNu)
    
    for(x in q){
      val1 <- s1$par[[x]]
      nam1 <- names(s1$par[x])
      val2 <- s2$par[[x]]
      nam2 <- names(s2$par[x])
      val3 <- s3$par[[x]]
      nam3 <- names(s3$par[x])
      val4 <- s4$par[[x]]
      nam4 <- names(s4$par[x])
      
      p1[[nam1]] <- val1
      p2[[nam2]] <- val2
      p3[[nam3]] <- val3
      p4[[nam4]] <- val4
    }
    
    p <- list( p1=p1, p2=p2, p3=p3, p4=p4)
    
    # Inference
    
    inter <- paste("fit_ALL_", modelName, ".rds", sep="")
    fit <- stan(file=modelName, data = data_multi, iter = 3000, warmup = 1000, chains = 4, init = p, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
    # Save fit result as rds file
    fit.dso <- new("cxxdso")
    saveRDS(fit, file = inter)
    
    
  } else if (infP=="single"){
    
      for(x in experiments){
        # Data Extraction (script MultiExtractExp with the function needs to be run first)
        data_extraction_multiexperiment(x)
        
        # Initial value optimisation
        m <- stan_model(modelName)
        s1 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
        s2 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
        s3 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
        s4 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
        
        # Refinment of the initial value guess
        while (s1$value < OptimVal || s2$value < OptimVal || s3$value < OptimVal || s4$value < OptimVal){
          if (s1$value < OptimVal){
            print("s1")
            s1 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
          }
          else if (s2$value < OptimVal){
            print("s2")
            s2 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
          }
          else if (s3$value < OptimVal){
            print("s3")
            s3 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
          }
          else if (s4$value < OptimVal){
            print("s4")
            s4 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
          }
          
        }
        
        q <- c(1:parNu)
        
        # Formating initial values to be passed to the stan fit function as for different lists (one for each chain)
        p1 <- list()
        p2 <- list()
        p3 <- list()
        p4 <- list()
        
        for(i in q){
          val1 <- s1$par[[i]]
          nam1 <- names(s1$par[i])
          val2 <- s2$par[[i]]
          nam2 <- names(s2$par[i])
          val3 <- s3$par[[i]]
          nam3 <- names(s3$par[i])
          val4 <- s4$par[[i]]
          nam4 <- names(s4$par[i])
          
          p1[[nam1]] <- val1
          p2[[nam2]] <- val2
          p3[[nam3]] <- val3
          p4[[nam4]] <- val4
        }
        
        p <- list( p1=p1, p2=p2, p3=p3, p4=p4)
        
        # Inference
        
        inter <- paste("fit_", x,"_", modelName, ".rds", sep="")
        fit <- stan(file=modelName, data = data_multi, iter = 2300, warmup = 300, chains = 4, init = p, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
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
