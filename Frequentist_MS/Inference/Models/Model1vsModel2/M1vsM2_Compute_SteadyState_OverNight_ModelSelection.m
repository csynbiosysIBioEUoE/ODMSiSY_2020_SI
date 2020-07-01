function [InitialConditions]= M1vsM2_Compute_SteadyState_OverNight_ModelSelection(inputs,params,InitialExpData,Initial_u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate the initial state of the ToggleSwitch model after the
% ON inclubation. We assume that the system is at steady state
% concentration specified in Run_MPLacr_in_silico_experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      

% Inputs:
%       - Inputs: inputs matlab structure required by AMIGO (paths, plot
%       and model need to be defined)
%       - params: theta vector to be used for the simulation
%       - InitialExpData: Experimental values for RFP, GFP, RFP and GFP at the first
%       time point of the experiment
%       - Initial_u: IPTG and aTc values used during the ON incuvation.


% Fixed parts of the experiment
duration = 24*60;     % Duration of the experiment in min

clear newExps;
newExps.n_exp = 1;
newExps.n_obs{1}=4;                                  % Number of observables per experiment                         
newExps.obs_names{1} = char('LacI_M1','TetR_M1', 'LacI_M3','TetR_M3');
newExps.obs{1} = char('LacI_M1 = L_RFP','TetR_M1 = T_GFP', 'LacI_M3 = L_RFP2','TetR_M3 = T_GFP2');% Name of the observables 
newExps.exp_y0{1}=M1vsM2_compute_steady_state_Analytical_ModelSelection(params,InitialExpData,Initial_u);     

newExps.t_f{1}=duration;               % Experiment duration
    
newExps.u_interp{1}='sustained';
newExps.u{1}=[Initial_u(1); Initial_u(2)];
newExps.t_con{1}=[0,duration];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mock an experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


inputs.exps = newExps;
inputs.plotd.plotlevel='noplot';

% SIMULATION
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='cvodes';
inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;
   
% AMIGO_Prep(inputs);
sim = AMIGO_SModel_NoVer(inputs);

InitialConditions = sim.sim.states{1,1}(end,:);

end