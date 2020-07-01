
%% Posteior distribution
mdl=1;
% Read Posterior
d2 = 'F:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\';
posterior = readtable([d2, 'PosteriorDistributionsBayes\draws_ALL_Model',num2str(mdl),'.csv']);
% Extract data
thetaNames = posterior.Properties.VariableNames;
thetaMatrix = posterior.Variables;

if ~isfolder('ReultsSimulations')
   mkdir('ReultsSimulations') 
end
%% Simulation Model 1 Training

% Training
Dat = load('F:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\LugagneData\BayesSelectDataLugagne.mat');

% Validation
% Dat = load('F:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\LugagneData\ValidationSetDataLugagne.mat');

Mod1 = ToggleSwitch_load_model_M1;

% Model 1
clear model;
clear exps;
clear best_global_theta;
clear pe_results;
clear pe_inputs;

% Specify folder name and short_name
results_folder = strcat('MTS',datestr(now,'yyyy-mm-dd-HHMMSS'));
short_name     = strcat('MTS_Statistics');

clear inputs;
inputs.model = Mod1;
inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = 'initial_setup';


AMIGO_Prep(inputs)
clear newExps;
newExps.n_exp = length(Dat.Data.exp_data);                                         % Number of experiments 
for i=1:length(Dat.Data.exp_data)
    
    y0 = M1_Compute_SteadyState_OverNight(inputs,mean(thetaMatrix),Dat.Data.exp_data{i}(end,:),[Dat.Data.Initial_IPTG{i} Dat.Data.Initial_aTc{i}]+1e-7);
    
    newExps.n_obs{i}=2;                                        % Number of observables per experiment        
    newExps.obs_names{i} = char('LacI_M1','TetR_M1');
    newExps.obs{i} = char('LacI_M1 = L_RFP','TetR_M1 = T_GFP');% Name of the observables 
            
    newExps.exp_y0{i}=y0;                                      % Initial condition for the experiment    

    newExps.t_f{i}=round(Dat.Data.t_samples{i}(end));                                   % Experiment duration
    newExps.n_s{i}=length(Dat.Data.t_samples{i});                             % Number of sampling times
    newExps.t_s{i}=round(Dat.Data.t_samples{i}(:,1))';                              % Times of samples

    newExps.u_interp{i}='step';                                % Interpolating function for the input
    % newExps.u_interp{2}='step';                                % Interpolating function for the input
    newExps.n_steps{i}=length(Dat.Data.input{i}(:,1));                  % Number of steps in the input
    newExps.u{i}= Dat.Data.input{i}'+1e-7;                                     % IPTG and aTc values for the input
    newExps.t_con{i}=round(Dat.Data.t_con{i}(:,1)');                     % Switching times
    %     newExps.t_con{1}=[0 3000];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mock the experiment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
    
inputs.exps = newExps;

inputs.exps.data_type='pseudo';
inputs.exps.noise_type='hetero_proportional';
inputs.exps.std_dev{1}=[0.0 0.0];

% SIMULATION
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='cvodes';
inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;

inputs.plotd.plotlevel='noplot';

inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = strcat('sim-Comp');

% simDa = AMIGO_SData(inputs);
AMIGO_Prep(inputs);

Simulations=cell(1,800);
for i=1:800
    inputs.model.par = thetaMatrix(i,:);
    simMod1 = AMIGO_SModel_NoVer(inputs);
    Simulations{1,i} = simMod1.sim.states;
end
save('ReultsSimulations\SimulationsBayesTrainingSetModel1.mat', 'Simulations')

mkdir('ReultsSimulations\PlotsModel1')
for j=1:10
    h = figure;
    subplot(6,1,1:2)
    hold on
    errorbar(simMod1.sim.tsim{j}, Dat.Data.exp_data{j}(:,1), Dat.Data.standard_dev{j}(:,1), 'r')
    for i=1:800
        plot(simMod1.sim.tsim{j}, Simulations{1,i}{j}(:,3), 'b')
    end
    ylabel('RFP')
    xlabel('time(min)')
    title(num2str(j))
    
    subplot(6,1,3:4)
    hold on
    errorbar(simMod1.sim.tsim{j}, Dat.Data.exp_data{j}(:,2), Dat.Data.standard_dev{j}(:,2), 'g')
    for i=1:800
        plot(simMod1.sim.tsim{j}, Simulations{1,i}{j}(:,4), 'b')
    end
    ylabel('GFP')
    xlabel('time(min)')
    
    subplot(6,1,5)
    hold on
    stairs(Dat.Data.t_con{j}(:,1), [Dat.Data.input{j}(:,1); Dat.Data.input{j}(end,1)])
    ylabel('IPTG')
    xlabel('time(min)')
    
    subplot(6,1,6)
    hold on
    stairs(Dat.Data.t_con{j}(:,1), [Dat.Data.input{j}(:,2); Dat.Data.input{j}(end,2)])
    ylabel('aTc')
    xlabel('time(min)')
    
    saveas(h, ['ReultsSimulations\PlotsModel1\BayesSimulationModel1_TrainingSet_',num2str(j),'.png'])
end



%% Simulation Model 1 Validation

% Training
% Dat = load('F:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\LugagneData\BayesSelectDataLugagne.mat');

% Validation
Dat = load('F:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\LugagneData\ValidationSetDataLugagne.mat');

Mod1 = ToggleSwitch_load_model_M1;

% Model 1
clear model;
clear exps;
clear best_global_theta;
clear pe_results;
clear pe_inputs;

% Specify folder name and short_name
results_folder = strcat('MTS',datestr(now,'yyyy-mm-dd-HHMMSS'));
short_name     = strcat('MTS_Statistics');

clear inputs;
inputs.model = Mod1;
inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = 'initial_setup';


AMIGO_Prep(inputs)
clear newExps;
newExps.n_exp = length(Dat.Data.exp_data);                                         % Number of experiments 
for i=1:length(Dat.Data.exp_data)
    
    y0 = M1_Compute_SteadyState_OverNight(inputs,mean(thetaMatrix),Dat.Data.exp_data{i}(end,:),[Dat.Data.Initial_IPTG{i} Dat.Data.Initial_aTc{i}]+1e-7);
    
    newExps.n_obs{i}=2;                                        % Number of observables per experiment        
    newExps.obs_names{i} = char('LacI_M1','TetR_M1');
    newExps.obs{i} = char('LacI_M1 = L_RFP','TetR_M1 = T_GFP');% Name of the observables 
            
    newExps.exp_y0{i}=y0;                                      % Initial condition for the experiment    

    newExps.t_f{i}=round(Dat.Data.t_samples{i}(end));                                   % Experiment duration
    newExps.n_s{i}=length(Dat.Data.t_samples{i});                             % Number of sampling times
    newExps.t_s{i}=round(Dat.Data.t_samples{i}(:,1))';                              % Times of samples

    newExps.u_interp{i}='step';                                % Interpolating function for the input
    % newExps.u_interp{2}='step';                                % Interpolating function for the input
    newExps.n_steps{i}=length(Dat.Data.input{i}(:,1));                  % Number of steps in the input
    newExps.u{i}= Dat.Data.input{i}'+1e-7;                                     % IPTG and aTc values for the input
    newExps.t_con{i}=round(Dat.Data.t_con{i}(:,1)');                     % Switching times
    %     newExps.t_con{1}=[0 3000];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mock the experiment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
    
inputs.exps = newExps;

inputs.exps.data_type='pseudo';
inputs.exps.noise_type='hetero_proportional';
inputs.exps.std_dev{1}=[0.0 0.0];

% SIMULATION
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='cvodes';
inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;

inputs.plotd.plotlevel='noplot';

inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = strcat('sim-Comp');

% simDa = AMIGO_SData(inputs);
AMIGO_Prep(inputs);

Simulations=cell(1,800);
for i=1:800
    inputs.model.par = thetaMatrix(i,:);
    simMod1 = AMIGO_SModel_NoVer(inputs);
    Simulations{1,i} = simMod1.sim.states;
end
save('ReultsSimulations\SimulationsBayesValidationSetModel1.mat', 'Simulations')

mkdir('ReultsSimulations\PlotsModel1')
for j=1:10
    h = figure;
    subplot(6,1,1:2)
    hold on
    errorbar(simMod1.sim.tsim{j}, Dat.Data.exp_data{j}(:,1), Dat.Data.standard_dev{j}(:,1))
    for i=1:800
        plot(simMod1.sim.tsim{j}, Simulations{1,i}{j}(:,3), 'r')
    end
    ylabel('RFP')
    xlabel('time(min)')
    title(num2str(j))
    
    subplot(6,1,3:4)
    hold on
    errorbar(simMod1.sim.tsim{j}, Dat.Data.exp_data{j}(:,2), Dat.Data.standard_dev{j}(:,2))
    for i=1:800
        plot(simMod1.sim.tsim{j}, Simulations{1,i}{j}(:,4), 'g')
    end
    ylabel('GFP')
    xlabel('time(min)')
    
    subplot(6,1,5)
    hold on
    stairs(Dat.Data.t_con{j}(:,1), [Dat.Data.input{j}(:,1); Dat.Data.input{j}(end,1)])
    ylabel('IPTG')
    xlabel('time(min)')
    
    subplot(6,1,6)
    hold on
    stairs(Dat.Data.t_con{j}(:,1), [Dat.Data.input{j}(:,2); Dat.Data.input{j}(end,2)])
    ylabel('aTc')
    xlabel('time(min)')
    
    saveas(h, ['ReultsSimulations\PlotsModel1\BayesSimulationModel1_ValidationSet_',num2str(j),'.png'])
end





















