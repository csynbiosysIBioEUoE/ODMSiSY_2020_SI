%% Script used to assess the best bayesian OED in the frequentist scheme
% It compares it with the best frequentist case in order to assess the
% differences between the two approaches (focus on posterior uncertainty in
% bayesian case, lost in here). 


%% Define path of the results
foldw = 'E:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\OED_ModelSelection\Results\'; % Main directory
fold = [foldw, '\Final\F_M1vsM3\']; % Folder containing the results
tag = 'OED1_M1vsM3_ThetaAMIGOlsq_Rep1-Steps_'; % Common file name for all the runs
steps = [2,4,6,8]; % Number of steps for each run (needed)
trial = 'M1vsM3_Bistability_'; % Tag that will be added in the name of the saved resulting plots

%% Frequentist Input
i=1;
x = load([fold, tag, num2str(steps(i)), '-MultOED_ModelSelection.mat']); % Here used best 2 step in frequentist for comparison, however it can be changed 
disp('IPTG')
disp(x.oed_results.do.u(1,:))
disp(' ')
disp('aTc')
disp(x.oed_results.do.u(2,:))

%% Bayesian Input
bu = [0.4300001, 0.4500001; 20.0000001, 22.0000001]; % Hard coded from Bayesian results. See python scripts for this.

%% Simulations
% Simulations are needed
clear model;
clear exps;
clear best_global_theta;
clear pe_results;
clear pe_inputs;
clear inputs;
% Read the model into the model variable
model = ToggleSwitch_load_model_M1vsM3_ModelSelection;
global_theta_guess = model.par;
% Start with no experiments
exps.n_exp=2;
global_theta_guess = global_theta_guess';
% Compile the model
inputs.model = model;
inputs.pathd.results_folder = 'TestSimulation';                        
inputs.pathd.short_name     = 'TS1_Random';
inputs.pathd.runident       = 'initial_setup';
AMIGO_Prep(inputs)

y0 = M1vsM3_Compute_SteadyState_OverNight_ModelSelection(inputs,global_theta_guess,[23, 1400, 23, 1400],[1 0]+1e-7);

% Fixed parts of the experiment
duration = 24*60;               % Duration in of the experiment (minutes)
clear newExps;
newExps.n_exp = 2;                                         % Number of experiments 

for i=1:2
newExps.n_obs{i}=2;                                        % Number of observables per experiment        
newExps.obs_names{i} = char('LacI_M1','TetR_M1', 'LacI_M3','TetR_M3');
newExps.obs{i} = char('LacI_M1 = L_RFP','TetR_M1 = T_GFP', 'LacI_M3 = L_RFP2','TetR_M3 = T_GFP2');% Name of the observables 
newExps.exp_y0{i}=y0;                                      % Initial condition for the experiment    
newExps.t_f{i}=duration;                                   % Experiment duration
newExps.n_s{i}=duration/5 + 1;                             % Number of sampling times
newExps.t_s{i}=0:5:duration ;                              % Times of samples
newExps.u_interp{i}='step';                                % Interpolating function for the input

end




inputs.model = model;
inputs.exps = newExps;
inputs.exps.data_type='pseudo';
inputs.exps.noise_type='hetero_proportional';
inputs.exps.std_dev{1}=[0.0 0.0];
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='fdsens5';
inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;
inputs.plotd.plotlevel='noplot';

inputs.exps.n_steps{1} = 2;
inputs.exps.u{1} = x.oed_results.do.u;
inputs.exps.t_con{1} = x.oed_results.do.t_con;

inputs.exps.n_steps{2} = 2;
inputs.exps.u{2} = bu;
inputs.exps.t_con{2} = 0:(duration/2):duration;


simMo = AMIGO_SModel(inputs);

%% Cost Function
% Compute both cost function values
[ff,~,~] = OEDMSCostMultiplic_ComparisonBF(simMo.sim.states{1});
[fb,~,~] = OEDMSCostMultiplic_ComparisonBF(simMo.sim.states{2});

disp(' ')
disp(['The frequentist result has a cost function value of: ', num2str(ff)])
disp(' ')
disp(['The bayesian result has a cost function value of: ', num2str(fb)])
disp(' ')

if ff<fb
    disp(['The frequentist result is ',num2str(ff/fb),' times better'])
elseif ff>fb
    disp(['The bayesian result is ',num2str(fb/ff),' times better'])
end
disp(' ')


%% Plot Results
% Frequentist
h = figure('Renderer', 'painters', 'Position', [50 50 900 600]);
subplot(3,1,1)
title(['Multiplicative optimisation, ', num2str(steps(i)), ' Steps, Best Frequentist'])
hold on
plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,3), 'r')
plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,7), 'b')
ylim([0,2500])
ylabel('RFP')
legend('M1', 'M3')

subplot(3,1,2)
hold on
plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,4), 'g')
plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,8), 'b')
ylim([0,1500])
ylabel('GFP')
legend('M1', 'M3')

subplot(3,1,3)
yyaxis left
stairs(inputs.exps.t_con{1}, [inputs.exps.u{1}(1,:), inputs.exps.u{1}(1,end)])
ylabel('IPTG')
yyaxis right
stairs(inputs.exps.t_con{1}, [inputs.exps.u{1}(2,:), inputs.exps.u{1}(2,end)])
ylabel('aTc')
xlabel('time (min)')

% Bayesian
h = figure('Renderer', 'painters', 'Position', [50 50 900 600]);
subplot(3,1,1)
title(['Multiplicative optimisation, ', num2str(steps(1)), ' Steps, Best Bayesian'])
hold on
plot(simMo.sim.tsim{2}, simMo.sim.states{2}(:,3), 'r')
plot(simMo.sim.tsim{2}, simMo.sim.states{2}(:,7), 'b')
ylim([0,2500])
ylabel('RFP')
legend('M1', 'M3')

subplot(3,1,2)
hold on
plot(simMo.sim.tsim{2}, simMo.sim.states{2}(:,4), 'g')
plot(simMo.sim.tsim{2}, simMo.sim.states{2}(:,8), 'b')
ylim([0,1500])
ylabel('GFP')
legend('M1', 'M3')

subplot(3,1,3)
yyaxis left
stairs(inputs.exps.t_con{2}, [inputs.exps.u{2}(1,:), inputs.exps.u{2}(1,end)])
ylabel('IPTG')
yyaxis right
stairs(inputs.exps.t_con{2}, [inputs.exps.u{2}(2,:), inputs.exps.u{2}(2,end)])
ylabel('aTc')
xlabel('time (min)')











