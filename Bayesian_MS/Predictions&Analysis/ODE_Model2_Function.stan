// ------------------------- TOGGLE SWITCH MODEL ------------------------- //


// Stan model script containing the inducer exchange model for the intermediate Toggle Switch where there is a single difussion rate for the inducers but NO degradation due to cell growth.
// The script is used to simulate the response of the model (ODEs) for a determined event based experimental profile


functions{
  //ODE
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
  
  matrix solve_coupled_ode(real[] ts, real[] p, real[] x_r, int[] x_i, int[] sp, real[] inputs, real[] toni, real[] ivss, real[] pre){
    int maxtime = num_elements(ts);
    int Nsp = num_elements(sp);
    int Nevents = num_elements(sp)-1;
    int tonil = num_elements(toni);
    int Neq = 4;
    
    // matrix[maxtime,Neq] total;
    real final[maxtime,Neq];
    real initialV[Neq];
    int i;
    vector[2] y_al; // Vector that will include the solutio of the algebraic solution for the steady state of the model
    real ssv[tonil,Neq];
    real y0[Neq]; // Initial values for the ODEs variables at the first event
    
    // Calculation of initial guesses
  y_al = SteadyState(to_vector(ivss), to_vector(p), pre, x_i); // Calculation of initial guesses for steady state
  y0[1] = pre[1];
  y0[2] = pre[2];
  y0[3] = y_al[1];
  y0[4] = y_al[2];
  
  ssv = integrate_ode_rk45(Toogle_one, y0,0,toni,p,pre, x_i, 1e-9, 1e-9, 1e7); // ON incubation calculation for the steady state
  y0 = ssv[tonil];
  
    
    initialV = y0;
    i = 1;
    
    for (q in 1:Nevents){
      int itp = sp[q];  // General way to extract the initial time points of each event
      int lts = num_elements(ts[(sp[q]+1):sp[q+1]+1]);  // General way to define the number of elements in each event series
      real input = inputs[q]; // General way to extract the input values
      real Tevent[lts] = ts[(sp[q]+1):sp[q+1]+1];  // General way to extract the times of each event
      
      real part1[lts,Neq];
      real temp[lts,Neq];
      
      
      if (q == 1){part1 = integrate_ode_rk45(Toogle_one, initialV,itp,ts[(sp[q]+1):sp[q+1]+1],p,to_array_1d(inputs[i:(i+1)]), x_i, 1e-9, 1e-9, 1e7);}
      else{part1 = integrate_ode_rk45(Toogle_one, initialV,(itp-1e-7),ts[(sp[q]+1):sp[q+1]+1],p,to_array_1d(inputs[i:(i+1)]), x_i, 1e-9, 1e-9, 1e7);}  
      
      initialV = part1[lts];
      i=i+2;
      
      for (y in (itp+1):(itp+lts)){
        
        final[(y),]=(part1)[(y-itp),];
      };
      
    };
    
    return(to_matrix(final));
  }
}
