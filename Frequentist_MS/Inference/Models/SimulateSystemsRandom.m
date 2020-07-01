%% Simulate random step experiment
% This script simulate a step experiment with random input profile for a
% selected model (Model that is a combination of Model1 and 3 is also
% included in the scheme)

% Inputs:
%       - epccOutputResultFileNameBase: String as an identifier for the
%       results
%       - mdl: Integer from 1 to 4 to indicate the model that is desired to
%       be simulated. 4 indicated the combination of model 1 and 3. 


function [simMo] = SimulateSystemsRandom(epccOutputResultFileNameBase, mdl)

epcc_exps = 1;

resultFileName = [strcat(epccOutputResultFileNameBase),'.dat'];
rng shuffle;
rngToGetSeed = rng;

% Write the header information of the .dat file . 
fid = fopen(resultFileName,'w');
fprintf(fid,'HEADER DATE %s\n',datestr(datetime()));
fprintf(fid,'HEADER RANDSEED %d\n',rngToGetSeed.Seed);
fclose(fid);

startTime = datenum(now);

clear model;
clear exps;
clear best_global_theta;
clear pe_results;
clear pe_inputs;

% Specify folder name and short_name
results_folder = strcat('MTS',datestr(now,'yyyy-mm-dd-HHMMSS'));
short_name     = strcat('MTS_CIC_Random',int2str(epcc_exps));

% Read the model into the model variable
switch mdl
    case 1
        model = ToggleSwitch_load_model_M1;
    case 2
        model = ToggleSwitch_load_model_M2;
    case 3
        model = ToggleSwitch_load_model_M3;
    case 4
        model = ToggleSwitch_load_model_M1vsM3_ModelSelection;
end
global_theta_guess = model.par;

% Start with no experiments
exps.n_exp=0;

global_theta_guess = global_theta_guess';

% Compile the model
clear inputs;
inputs.model = model;
inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = 'initial_setup';

switch mdl
    case 1
        y0 = M1_Compute_SteadyState_OverNight(inputs,global_theta_guess,[23, 1400],[1 0]+1e-7);
    case 2
        y0 = M2_Compute_SteadyState_OverNight(inputs,global_theta_guess,[23, 1400],[1 0]+1e-7);
    case 3
        y0 = M3_Compute_SteadyState_OverNight(inputs,global_theta_guess,[23, 1400],[1 0]+1e-7);
    case 4
        y0 = M1vsM3_Compute_SteadyState_OverNight_ModelSelection(inputs,global_theta_guess,[23, 1400, 23, 1400],[1 0]+1e-7);
end

% Fixed parts of the experiment
duration = 24*60;               % Duration in of the experiment (minutes)
Stepd = 180;                                         % Duration of each step (minutes). Note that this value equals the response time, quantified in 80 mins for MPLac,r
inducer = [rand(1,round(duration/Stepd)); randi([0 100],1,round(duration/Stepd))]+1e-7;  % Extract 2 random integer values from the Uniform distribution in (0,1) ng/mL for IPTG and (0,100) for aTc to be used as low and high values in each cycle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a new experiment to simulate with the step input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newExps;
newExps.n_exp = 1;                                         % Number of experiments 
newExps.n_obs{1}=2;                                        % Number of observables per experiment        
switch mdl
    case 1
        newExps.obs_names{1} = char('LacI_M1','TetR_M1');
        newExps.obs{1} = char('LacI_M1 = L_RFP','TetR_M1 = T_GFP');% Name of the observables 
    case 2
        newExps.obs_names{1} = char('LacI_M2','TetR_M2');
        newExps.obs{1} = char('LacI_M2 = L_RFP','TetR_M2 = T_GFP');% Name of the observables 
    case 3
        newExps.obs_names{1} = char('LacI_M3','TetR_M3');
        newExps.obs{1} = char('LacI_M3 = L_RFP','TetR_M3 = T_GFP');% Name of the observables 
    case 4
        newExps.obs_names{1} = char('LacI_M1','TetR_M1', 'LacI_M3','TetR_M3');
        newExps.obs{1} = char('LacI_M1 = L_RFP','TetR_M1 = T_GFP', 'LacI_M3 = L_RFP2','TetR_M3 = T_GFP2');% Name of the observables 
end

newExps.exp_y0{1}=y0;                                      % Initial condition for the experiment    

newExps.t_f{1}=duration;                                   % Experiment duration
newExps.n_s{1}=duration/5 + 1;                             % Number of sampling times
newExps.t_s{1}=0:5:duration ;                              % Times of samples

newExps.u_interp{1}='step';                                % Interpolating function for the input
newExps.n_steps{1}=round(duration/Stepd);                  % Number of steps in the input
newExps.u{1}= inducer;                                     % IPTG and aTc values for the input
newExps.t_con{1}=[0:180:1440];                     % Switching times

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mock the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear inputs;
inputs.model = model;
inputs.exps = newExps;

inputs.exps.data_type='pseudo';
inputs.exps.noise_type='hetero_proportional';
inputs.exps.std_dev{1}=[0.0 0.0];

% SIMULATION
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='cvodes';
inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;

inputs.plotd.plotlevel='medium';

inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = strcat('sim-',int2str(epcc_exps));

AMIGO_Prep(inputs);

% simDa = AMIGO_SData(inputs);
simMo = AMIGO_SModel(inputs);



