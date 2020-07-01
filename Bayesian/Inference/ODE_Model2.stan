// ------------------------- TOGGLE SWITCH STAN MODEL 3 ------------------------- //


// Stan model script containing the inducer exchange model for the intermediate Toggle Switch where there is a single difussion rate for the inducers but NO degradation due to cell growth. 
// The script can be used for inference on a single experimental result or for a multiexperimental inference

functions{
  // Function containing the ODE to be used for the inference
  real[] Toogle_one(real t, real[] y, real[] p, real[] x_r, int[] x_i){
    // Inputs
    real u_IPTG = x_r[1];
    real u_aTc = x_r[2];
    // Parameters
    real k_IPTG = p[1];
    real k_aTc = p[2];
    real k_L_pm0 = p[3];
    real k_L_pm = p[4];
    real theta_T = p[5];
    real theta_aTc = p[6];
    real n_aTc = p[7];
    real n_T = p[8];
    real k_T_pm0 = p[9];
    real k_T_pm = p[10];
    real theta_L = p[11];
    real theta_IPTG = p[12];
    real n_IPTG = p[13];
    real n_L = p[14];
    
    //Equations
    real dInd_dt[4];
    
    dInd_dt[1] = k_IPTG*(x_r[1]-y[1]);
    dInd_dt[2] = k_aTc*(x_r[2]-y[2]);

    dInd_dt[3] = ((1/0.1386)*(k_L_pm0+(k_L_pm/(1+(y[4]/theta_T*1/(1+(y[2]/theta_aTc)^n_aTc))^n_T))))-0.0165*y[3];
    dInd_dt[4] = ((1/0.1386)*(k_T_pm0+(k_T_pm/(1+(y[3]/theta_L*1/(1+(y[1]/theta_IPTG)^n_IPTG))^n_L))))-0.0165*y[4];
    
    // RESULTS
    return dInd_dt;
  }
  
  // Function type vector containing the equations where the root needs to be calculated for the steady states
  vector SteadyState(vector init, vector p, real[] x_r, int[] x_i){
    vector[2] alpha;
    // Parameters
    real k_IPTG = p[1];
    real k_aTc = p[2];
    real k_L_pm0 = p[3];
    real k_L_pm = p[4];
    real theta_T = p[5];
    real theta_aTc = p[6];
    real n_aTc = p[7];
    real n_T = p[8];
    real k_T_pm0 = p[9];
    real k_T_pm = p[10];
    real theta_L = p[11];
    real theta_IPTG = p[12];
    real n_IPTG = p[13];
    real n_L = p[14];
    // Equations
    alpha[1] = ((1/0.1386)*(k_L_pm0+(k_L_pm/(1+(init[2]/theta_T*1/(1+(x_r[2]/theta_aTc)^n_aTc))^n_T))))/0.0165;
    
    alpha[2] = ((1/0.1386)*(k_T_pm0+(k_T_pm/(1+(init[1]/theta_L*1/(1+(x_r[1]/theta_IPTG)^n_IPTG))^n_L))))/0.0165;
    // Results
    return alpha;
  }
}

data {
    // Observables
    int m; // Total number of data series
    int stslm; // Maximum number of rows for all the observable matrices
    int stsl[1,m]; // Number of elements at each time series for each series m
    int sts[stslm,m]; // Sampling times for each series m
    real GFPmean[stslm,m]; // estimated observables for tetR+GFP at each sampling time
    real RFPmean[stslm,m]; // estimated observables for LacI+RFP at each sampling time
    real GFPstd[stslm,m]; // standard error for tetR+GFP at each sampling time
    real RFPstd[stslm,m]; // standard error for LacI+RFP at each sampling time
    
    // Inputs
    int elm; // Number of rows in the matrices for IPTG and aTc, half for the inputs and -1 for the total number of events
    int tml; // Maximum length of the rows for the sampling times matrix
    int Nsp[1,m]; // Number of event switching points (including final time) times for each series m
    real ts[tml, m]; // Time series for each serie m
    int tsl[1,m]; // Length of sampling time series per serie m
    real preIPTG[1,m]; // Values of inputs for each serie m for the ON incubation 
    real preaTc[1,m];
    real IPTG[elm,m]; // Values of inputs at each event for each serie m
    real aTc[elm,m];
    real inputs[(elm*2),m]; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
    int evnT[(elm+1),m]; // Event change time points for each serie m
    
    // Over night incubation times
    int tonil;
    real toni[tonil];
}

transformed data {
  int nParms = 14; // Number of parameters of the model
  int Neq = 4; // Total number of equations of the model
  int x_i[0]; // Empty x_i object (needs to be deffined)
  real x_r[(elm*2),m]=inputs; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
  real ivss[Neq-2,m]; // Initial experimental values for the calculation of the steady state ordered as LacI+RFP, TetR+GFP
  real pre[2,m]; // Input values during the ON incubation ordered as IPTG, aTc
  
  for(i in 1:m){
    ivss[1,i] = RFPmean[1,i];
    ivss[2,i] = GFPmean[1,i];
    pre[1,i] = preIPTG[1,i];
    pre[2,i] = preaTc[1,i];
  };
}

parameters {
    // Parameters to be infered in the model (unit scale)
    real<lower=-2.1,upper=2.1> k_IPTG_raw;
    real<lower=-2.1,upper=2.1> k_aTc_raw;
    real<lower=-2.1,upper=2.1> k_L_pm0_raw;
    real<lower=-2.1,upper=2.1> k_L_pm_raw;
    real<lower=-2.1,upper=2.1> theta_T_raw;
    real<lower=-2.1,upper=2.1> theta_aTc_raw;
    real<lower=-2.1,upper=2.1> n_aTc_raw;
    real<lower=-2.1,upper=2.1> n_T_raw;
    real<lower=-2.1,upper=2.1> k_T_pm0_raw;
    real<lower=-2.1,upper=2.1> k_T_pm_raw;
    real<lower=-2.1,upper=2.1> theta_L_raw;
    real<lower=-2.1,upper=2.1> theta_IPTG_raw;
    real<lower=-2.1,upper=2.1> n_IPTG_raw;
    real<lower=-2.1,upper=2.1> n_L_raw;
}

transformed parameters {
  // Introduction of the paramemters in an indexed object and transformation to real paramater values that will be introduced into the ODEs
  real theta[nParms];

  theta[1] = exp(((k_IPTG_raw)*(1.15129254649702))+(-3.2188758248682));
  theta[2] = exp(((k_aTc_raw)*(1.15129254649702))+(-2.30258509299405));
  theta[3] = exp(((k_L_pm0_raw)*(1.15129254649702))+(-3.50655789731998));
  theta[4] = exp(((k_L_pm_raw)*(1.15129254649702))+(2.30258509299405));
  theta[5] = exp(((theta_T_raw)*(1.15129254649702))+(3.40119738166216));
  theta[6] = exp(((theta_aTc_raw)*(1.15129254649702))+(2.30258509299405));
  theta[7] = (((n_aTc_raw)*(1.25))+(2.5));
  theta[8] = (((n_T_raw)*(1.25))+(2.5));
  theta[9] = exp((((k_T_pm0_raw))*(1.15129254649702))-2.30258509299405);
  theta[10] = exp(((k_T_pm_raw)*(1.15129254649702))+(2.22044604925031e-16));
  theta[11] = exp(((theta_L_raw)*(1.15129254649702))+3.40119738166216);
  theta[12] = exp(((theta_IPTG_raw)*(1.15129254649702))-2.30258509299405);
  theta[13] = ((n_IPTG_raw)*(1.25))+2.5;
  theta[14] = ((n_L_raw)*(1.25))+2.5;
}

model {
  
  // Intermediate parameters
  int i; // Increasing index for the inputs
  vector[2] ing; // Vector that will include the solutio of the algebraic solution for the steady state of the model
  real ssv[tonil,Neq]; // Real that will include the solution of the ODE for the ON incubation (24h)
  real Y0[Neq,m]; // Initial values for the ODEs variables at the first event
  real yhat3[stslm,m]; // Reals that will include the values of RFP and GFP over time separately
  real yhat4[stslm,m];
  
  // Priors definition (Unit Normal)
  k_IPTG_raw ~ normal(0,1);
  k_aTc_raw ~ normal(0,1);
  k_L_pm0_raw ~ normal(0,1);
  k_L_pm_raw ~ normal(0,1);
  theta_T_raw ~ normal(0,1);
  theta_aTc_raw ~ normal(0,1);
  n_aTc_raw ~ normal(0,1);
  n_T_raw ~ normal(0,1);
  k_T_pm0_raw ~ normal(0,1);
  k_T_pm_raw ~ normal(0,1);
  theta_L_raw ~ normal(0,1);
  theta_IPTG_raw ~ normal(0,1);
  n_IPTG_raw ~ normal(0,1);
  n_L_raw ~ normal(0,1);
  
  // Likelihood
  for (j in 1:m){
    real ivst[Neq]; // Initial value of the states 
    real y_hat[(tsl[1,j]),Neq]; // Were results will be included
    
    // Calculation of initial guesses
    ing = SteadyState(to_vector(ivss[1:2,j]), to_vector(theta), pre[1:2,j], x_i); // Calculation of initial guesses for steady state
    Y0[1,j] = preIPTG[1,j];
    Y0[2,j] = preaTc[1,j];
    Y0[3,j] = ing[1];
    Y0[4,j] = ing[2];
    ssv = integrate_ode_bdf(Toogle_one, Y0[,j],0,toni,theta,pre[1:2,j], x_i, 1e-9, 1e-9, 1e7); // ON incubation calculation for the steady state
    
    Y0[,j] = ssv[tonil];
    i = 1;
    
    // Loop (over the number of events) to solve the ODE for each event stopping the solver and add them to the final object y_hat
    for (q in 1:Nsp[1,j]-1){
      
      int itp = evnT[q,j];  // Initial time points of each event
      int lts = num_elements(ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j]);  // Length of the time series for each event
      real part1[lts,Neq]; // Temporary object that will include the solution of the ODE for each event at each loop
      // Calculation of the solution for the ODEs where for events that are not the firt one the time series starts one minute before the original point of the time serie overlaping with the last point of the previous event with same state values at the time
      if (q == 1){
        ivst = Y0[,j];
        part1 = integrate_ode_bdf(Toogle_one,ivst,itp,ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i:(i+1),j]), x_i, 1e-9, 1e-9, 1e7);
      }
      else{
        part1 = integrate_ode_bdf(Toogle_one, ivst,(itp-1e-7),ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i:(i+1),j]), x_i, 1e-9, 1e-9, 1e7);
      }

      // Modification of the initial state values for the next event
      ivst = part1[lts];
      // Increase index for inputs
      i=i+2;
      
      // Introduction of the result of part1 into the object y_hat
      for (y in (itp+1):(itp+lts)){
        y_hat[(y),]=(part1)[(y-itp),];
      };
    };
    
      // Likelihood at each sampling time
    for (t in 1:stsl[1,j]){
      yhat3[t,j] = y_hat[(sts[t,j]+1),3];
      yhat4[t,j] = y_hat[(sts[t,j]+1),4];
      RFPmean[t,j] ~ normal(yhat3[t,j],RFPstd[t,j]);
      GFPmean[t,j] ~ normal(yhat4[t,j],GFPstd[t,j]);
    }
  };
}

