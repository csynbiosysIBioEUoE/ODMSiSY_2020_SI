# Optimally Designed Model Selection in Synthetic Biology

## Here the scripts used to generate, analyse and visualise the data presented in the paper are made available. 

## The use of the scripts requires the RStan Package, available at the link: 
https://cran.r-project.org/web/packages/rstan/index.html 

## And the use of the Bayesian Optimisation package available at the link: 
https://github.com/fmfn/BayesianOptimization

## For the Matlab toolbox AMIGO2, refer to: 
https://sites.google.com/site/amigo2toolbox/

## The data for the *Bayesian* case is organised in the following subfolders (Bayesian Directory):

-	**Inference:**
	- *ODE_Model1.stan*, stan statistical model script with Model 1 (Lugagne et.al.) from the paper used to perform Bayesian Inference of the experimental data from [1].
	- *ODE_Model2.stan*, stan statistical model script with Model 2 (Intermediate) from the paper used to perform Bayesian Inference of the experimental data from [1].
	- *ODE_Model3.stan*, stan statistical model script with Model 3 (New) from the paper used to perform Bayesian Inference of the experimental data from [1].
	- *MultiExtractExp.R*, script designed to access the experimental data and experimental schemes from [1] to generate an appropriate list of objects to be passed to the stan model to perform the inference. The csv files are generated using the script DataExtraction.m.
	- *DataExtraction.m*, script to extract the desired experimental data and experimental profiles from [1]. 
	- *masterRun.R*, script to perform inference through RStan using the designed model ODE_Model.stan and the list of data extracted from MultiExtractExp.R. The script allows performing inference on single datasets in series or on the combined set. 
	- *masterRunOptim.R*, script designed as the masterRun.R script but including an initial optimisation process for the initialisation of the 4 MCMC chains used in the inference. 
  
-	**PriorDefinition:**
	- *ExtractingInitialPriorsLugagneLog.mat*, Matlab script to compute the mean and standard deviation of our priors (10 lognormal and 4 normal distributions) based on the results of the fit obtained in [1].

-   **Predictions&Analysis**
      - *ODE_Model1_Function.stan*, stan script containing the proposed ODE system and the implementation of the event-based representation of the inputs to simulate the response to a selected input (processed with the function MultiExtractExp.R). This Script is for Model 1 (Lugagne et.al.).
      - *ODE_Model2_Function.stan*, stan script containing the proposed ODE system and the implementation of the event-based representation of the inputs to simulate the response to a selected input (processed with the function MultiExtractExp.R). This Script is for Model 2 (Intermediate).
      - *ODE_Model3_Function.stan*, stan script containing the proposed ODE system and the implementation of the event-based representation of the inputs to simulate the response to a selected input (processed with the function MultiExtractExp.R). This Script is for Model 3 (New).
      -	*PostPredCheckSimulM1.R*, function to simulate the ODEs for a determined experimental profile selected using all the MCMC samples from a stanfit object result selected and save the results in CSV format. This Script is for Model 1 (Lugagne et.al.).
      -	*PostPredCheckSimulM2.R*, function to simulate the ODEs for a determined experimental profile selected using all the MCMC samples from a stanfit object result selected and save the results in CSV format. This Script is for Model 2 (Intermediate).
      -	*PostPredCheckSimulM3.R*, function to simulate the ODEs for a determined experimental profile selected using all the MCMC samples from a stanfit object result selected and save the results in CSV format. This Script is for Model 3 (New).
      -	*ConfidenceIntervalPlotsFunctionM1.R*, function to extract all the MCMC samples from a stanfit object and simulate a determined experimental profile, obtaining the 95% confidence intervals and saving the plot of the simulation for checks. This Script is for Model 1 (Lugagne et.al.).
      -	*ConfidenceIntervalPlotsFunctionM2.R*, function to extract all the MCMC samples from a stanfit object and simulate a determined experimental profile, obtaining the 95% confidence intervals and saving the plot of the simulation for checks. This Script is for Model 2 (Intermediate).
      -	*ConfidenceIntervalPlotsFunctionM3.R*, function to extract all the MCMC samples from a stanfit object and simulate a determined experimental profile, obtaining the 95% confidence intervals and saving the plot of the simulation for checks. This Script is for Model 3 (New).
      

-	**ModelComparison:**
      -	*GaussianMixtureCompM1.R*, R script used to produce a Gaussian Mixture fit on the posteriors obtained by the Rstan inference in model 1. 
      -	*GaussianMixtureCompM23.R*, R script used to produce a Gaussian Mixture fit on the posteriors obtained by the Rstan inference in model 2 and 3. 
      - *BayesPredictions.ipynb*, Python script used to assess the goodness of fit computing the nRMSE
	  - *BayesPredictions-Training.ipynb*, Python script used to assess the goodness of fit on unseen data computing the nRMSE
	  - *evidencecomputation_Each_Multi_M2_GMM*, R script to compute the approximation of the model evidence for Model 1 ussing the Laplace approximation. 
	  - *evidencecomputation_Each_Multi_M3_GMM*, R script to compute the approximation of the model evidence for Model 3 ussing the Laplace approximation. 
	  - *evidencecomputation_Each_Multi_M4_GMM*, R script to compute the approximation of the model evidence for Model 2 ussing the Laplace approximation. 
  
-	**BEDms:**
      -	*ExtractReducedThetaDraws.R*, R script to extract all the draws for the stanfit results on the three models and also generate separate files with only 400 randomly picked draws used for the optimisation. 
      - *OptimDataFilesGenerator.R*, R script used to generate the necessary CSV files to simulate the three different models according to any input setting obtained during the optimisation. 
      - *OptimisationScriptsSingle*, directory containing all the Jupyter Notebooks used to perform Bayesian Experimental Design for model selection with the default settings from the BayesianOptimisation package (kernel = Matern2.5, aqf = ucb) for experiments of 2,4,6 and 8 steps considering as a utility function the Bhattacharyya distance from RFP, GFP, the average between the two or the product between the two. 
      - *OptimisationScriptsVariants*, directory containing the Jupyter Notebooks used o perform Bayesian Experimental Design for model selection with varying kernels for the Gaussian Processes and/or different acquisition functions with different parameters for 4 step experiments. It also contains the script to perform Dynamic Bayesian Experimental Design for model selection on an 8 step experiment where the optimisation is performed every 2 steps. 
      - *Packages*, PDF file containing the packages used for Gaussian Process Regression and Bayesian Optimisation and where to access them. 
	  - *ComputationalImprovement*, directory containing an example script from OptimisationsScriptsSingle but with modifications (c compilations and parallelisation) to improve computational speed (as an example, done after all computations)
	  - *CompareMultiplicativeResults.ipynb*, Python script used to compare the Bayesian and Frequentist results in a Bayesian scheme (use to assess differences between the two procedures). 
      
	  
## The data for the *Frequentist* case is organised in the following subfolders (Frequentist Directory):  

-	**Inference:**
	- *Models*, directory containing all the Matlab/AMIGO functions containing the mathematical models of the Toggle Switch.
	- *AnalysisPECrossValidation.m*, Matlab script to extract the best run from the 100 trials using cross-validation.
	- *AnalysisPEResultsV2.m*, Matlab script to extract the best run from the 100 trials using the values of the cost function. 
	- *ExtractGaussParamsTheta.m*, Matlab script to extract the mean and covariance matrix from the Bayesian Posterior distributions. 
	- *Masterfit_to_AllModels.m*, main Matlab script to run Parameter Estimation (NO cross-validation considered).
	- *Masterfit_to_AllModels_CV.m*, main Matlab script to run Parameter Estimation (cross-validation considered).
	- *MasterRun_fit_to_AllModels.m*, Matlab script to run Masterfit_to_AllModels.m with desired inputs.
	- *MasterRun_fit_to_AllModels_CV.m*,  Matlab script to run Masterfit_to_AllModels_CV.m with desired inputs.
 
-	**PriorDefinition:**
	- *getThetaBounds.m*, Matlab script to get the parameter bounds for all parameters (hardcoded according to the Bayesian case). 

-	**ModelComparison:**
      -	  *AMIGO_PEPostAnalysisCustom.m*, Matlab code from AMIGO modified to compute all model comparison desired metrics and compute AIC for the individual experiments. 
	  -   *ModelComparisonStatistics_V2.m*, Matlab script used to simulate the best PE result (NO Cross-Validation) for the 3 models on the training data and compute a series of metrics to assess the goodness of fit. 
	  -   *ModelComparisonStatistics_V2_CV.m*, Matlab script used to simulate the best PE result (Cross-Validation) for the 3 models on the validation data and compute a series of metrics to assess the goodness of fit. 
	  -   *FrequentistPredictions.ipynb*, Python script used to simulate the validation experiments for the three models and compute the nRMSE. Results are used in Supplementary Figure 1. 
	  
-   **BEDms:**
	-   *AMIGOChanged*, directory containing modified AMIGO functions. The either add elements into the AMIGO scripts or delete console printing and pauses to increase speed on computations. 
	-   *CostFunctions*, directory containing the cost-function files for all the cases considered. A series of empty constraint files are also included since AMIGO requires them. 
	-   *CostFunctions_CompareBayesFreq*, directory containing modifications on the cost-function files so to be used outside the optimisation scheme, only introducing a matrix with the simulation results. 
	-   *AnalyseResultsOED.m*, Matlab script used to analyse the results from OED for model selection. 
	-   *ComparisonBayesianVsFrequentistAverage.m*, Matlab script used to compare the Bayesian and Frequentist results in a Frequentist scheme (use to assess differences between the two procedures). 
	-   *OEDModelSelection_Master.m*, main Matlab script to run the OED for model selection. 
	-   *OEDModelSelection_Master_ShortIter.m*, main Matlab script to run the OED for model selection where maximum number of iterations has been set to 800 to be comparable with the Bayesian counterpart. 
	-   *Run_OEDModelSelection_Master.m*, Matlab script used to run OEDModelSelection_Master.m or OEDModelSelection_Master_ShortIter.m with desired inputs. 


## The data for the *Stability Analysis* is organised in the following subfolders (Stability_Analysis Directory):

-   *AnalyticalSteadyState*, directory containing the Matlab function to compute the analytical steady-state for the three models. 
-   *CostUtilityComputation*, directory containing Python and Matlab files used for the analysis of the results of the section and computation of the Utility/Cost function for a series of given inputs for the competing models. 
-   *BistabilityAnalysisBayesian.m*, Matlab script used to compute the stability matrix in the Bayesian case for the three models. 
-   *BiStabilityAnalysisFrequentist.m*, Matlab script used to compute the stability matrix in the frequentist case for the three models. 
-   *InterX.m*, Link to the function/package used to compute the line intersection of the nullclines. 

## The data for the *Paper Figures* is organised in the following subfolders (Images Directory):

-   *Figure2d0_Size.ipynb*, python script used to generate figure 2 from the paper. 
-   *Figure3Sketch.ipynb*, python script used to generate figure 3 from the paper. 
-   *Figure5Optimisation2d0.ipynb*, python script used to generate figure 5 from the paper. 
-   *Figure6SketchVersion2.ipynb*, python script used to generate figure 6 from the paper. 
-   *FigureS1Sketch.ipynb*, python script used to generate supplementary figure 1 from the paper. 


To run the scripts plase place all the scripts and data required in the R, Matlab or Jupyter Notebook working directory. 

The data associated with these scripts can be found at:
TO MODIFY!
https://datasync.ed.ac.uk/index.php/s/oyrtYQSCoARaxhR (pwd: IAC_CIBCB2019_SI_Data)

References:

[1] Jean-Baptiste Lugagne, Sebastián Sosa Carrillo, Melanie Kirch, Agnes Köhler, Gregory Batt & Pascal Hersen, 2017. Balancing a genetic toggle switch by real-time feedback control and periodic forcing. Nature Communications, 8 (1671), pp. 1-7.
