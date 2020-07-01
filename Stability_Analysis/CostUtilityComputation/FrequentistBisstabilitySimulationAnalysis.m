%% Script to compute the OED for model selection cost function for a given set of inputs
% Script used to generate data for Figure 6

%% Inducers to be used in simulations
inducersrandom = [0.83, 10; 0.1, 47; 0.76, 86; 0.1, 7; 0.25 23; 0.08 14; 0.6, 18; 0.12 6;  0.74 16; 0.5 39; 0.94 51; 0.78 43]; % Randomly sampled points in stability matrix from Figure 6
indbest2sR = [0.2402, 26.1299; 0.4320, 22.7485]; %Best 2 Step Frequentist OED RFP
indbest2sG = [0.4415, 100.0000; 0.4873, 26.5901]; %Best 2 Step Frequentist OED GFP
indbest2sA = [0.2532, 26.6272; 0.4313, 22.7305]; %Best 2 Step Frequentist OED Average
indbest2sM = [0.2623, 27.3593; 0.4298, 22.7023]; %Best 2 Step Frequentist OED Multiplicative

% BEst frequentist OED for RFP, GFP, Average and Multiplicative
indbestR = [0.2626    0.5140    0.3634    0.3878    0.3764    0.3729    0.3744    0.4108; 37.0855   22.1488   22.4032   22.7415   22.7849   22.7849   22.7849   22.7849];
indbestG = [0.3992    0.5396    0.4261    0.3693    0.3483    0.3479    0.3484    0.3490; 99.9997   23.5026   20.3263   19.1514   19.4201   19.5584   19.8114   19.8374];
indbestA = [0.2736    0.5143    0.3761    0.3877    0.3813    0.3958    0.3962    0.4003; 38.2720   21.9904   22.4132   22.7342   22.7852   22.7852   22.7852   22.7853];
indbestM = [0.2844    0.5145    0.3679    0.3687    0.3669    0.3661    0.3542    0.3892; 40.1364   21.7241   22.4238   22.7250   22.7861   22.7862   22.7862   22.7862];


%% Simulation definition for random picks
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
newExps.n_exp = 12;                                         % Number of experiments 

for i=1:12
newExps.n_obs{i}=2;                                        % Number of observables per experiment        
newExps.obs_names{i} = char('LacI_M1','TetR_M1', 'LacI_M3','TetR_M3');
newExps.obs{i} = char('LacI_M1 = L_RFP','TetR_M1 = T_GFP', 'LacI_M3 = L_RFP2','TetR_M3 = T_GFP2');% Name of the observables 
newExps.exp_y0{i}=y0;                                      % Initial condition for the experiment    
newExps.t_f{i}=duration;                                   % Experiment duration
newExps.n_s{i}=duration/5 + 1;                             % Number of sampling times
newExps.t_s{i}=0:5:duration ;                              % Times of samples
newExps.u_interp{i}='step';                                % Interpolating function for the input
newExps.n_steps{i} = 1;
newExps.t_con{i} = [0, duration];
newExps.u{i} = inducersrandom(i,:)';
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

simMoRan = AMIGO_SModel_NoVer(inputs);
save('SimulationsBistability\RandomPicks.mat','simMoRan')

%% 2 Step simulations (Best Bayesians)

indbest2s = [indbest2sR;indbest2sG;indbest2sA;indbest2sM];

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
newExps.n_exp = 8;                                         % Number of experiments 

for i=1:8
newExps.n_obs{i}=2;                                        % Number of observables per experiment        
newExps.obs_names{i} = char('LacI_M1','TetR_M1', 'LacI_M3','TetR_M3');
newExps.obs{i} = char('LacI_M1 = L_RFP','TetR_M1 = T_GFP', 'LacI_M3 = L_RFP2','TetR_M3 = T_GFP2');% Name of the observables 
newExps.exp_y0{i}=y0;                                      % Initial condition for the experiment    
newExps.t_f{i}=duration;                                   % Experiment duration
newExps.n_s{i}=duration/5 + 1;                             % Number of sampling times
newExps.t_s{i}=0:5:duration ;                              % Times of samples
newExps.u_interp{i}='step';                                % Interpolating function for the input
newExps.n_steps{i} = 1;
newExps.t_con{i} = [0, duration];
newExps.u{i} = indbest2s(i,:)';
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

simMo2S = AMIGO_SModel_NoVer(inputs);
save('SimulationsBistability\2SPicks.mat','simMo2S')


%% Best OED simulations (Best Frequentist)

indbest = [indbestR;indbestG;indbestA;indbestM];

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
newExps.n_exp = 4;                                         % Number of experiments 

r=1:2:8;
for i=1:4
newExps.n_obs{i}=2;                                        % Number of observables per experiment        
newExps.obs_names{i} = char('LacI_M1','TetR_M1', 'LacI_M3','TetR_M3');
newExps.obs{i} = char('LacI_M1 = L_RFP','TetR_M1 = T_GFP', 'LacI_M3 = L_RFP2','TetR_M3 = T_GFP2');% Name of the observables 
newExps.exp_y0{i}=y0;                                      % Initial condition for the experiment    
newExps.t_f{i}=duration;                                   % Experiment duration
newExps.n_s{i}=duration/5 + 1;                             % Number of sampling times
newExps.t_s{i}=0:5:duration ;                              % Times of samples
newExps.u_interp{i}='step';                                % Interpolating function for the input
newExps.n_steps{i} = 8;
newExps.t_con{i} = 0:180:duration;
newExps.u{i} = indbest(r(i):r(i)+1,:);
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

simMoBest = AMIGO_SModel_NoVer(inputs);
save('SimulationsBistability\BestPicks.mat','simMoBest')


%% Plots

% Random picks
for i=1:12
    h=figure;
    subplot(6,1,1:2)
    hold on
    plot(simMoRan.sim.tsim{1}, simMoRan.sim.states{i}(:,3), 'r')
    plot(simMoRan.sim.tsim{1}, simMoRan.sim.states{i}(:,7), 'b')
    legend('Model 1', 'Model 2')
    ylabel('RFP(A.U.)')
    title(['Random ', num2str(i)])
    ylim([0,3500])
    
    subplot(6,1,3:4)
    hold on
    plot(simMoRan.sim.tsim{1}, simMoRan.sim.states{i}(:,4), 'g')
    plot(simMoRan.sim.tsim{1}, simMoRan.sim.states{i}(:,8), 'b')
    legend('Model 1', 'Model 2')
    ylabel('GFP(A.U.)')
    ylim([0,2500])
    
    subplot(6,1,5)
    stairs([0,duration],[inducersrandom(i,1), inducersrandom(i,1)])
    ylabel('IPTG')
    
    subplot(6,1,6)
    stairs([0,duration],[inducersrandom(i,2), inducersrandom(i,2)])
    ylabel('aTc')
    xlabel('time (min)')
    
    saveas(h, ['SimulationsBistability\Plots\','RandomPicksSimulations_',num2str(i),'.png'])
end



% Best 2S picks
for i=1:8
    h=figure;
    subplot(6,1,1:2)
    hold on
    plot(simMo2S.sim.tsim{1}, simMo2S.sim.states{i}(:,3), 'r')
    plot(simMo2S.sim.tsim{1}, simMo2S.sim.states{i}(:,7), 'b')
    legend('Model 1', 'Model 2')
    ylabel('RFP(A.U.)')
    title(['Random ', num2str(i)])
    ylim([0,3500])
    
    subplot(6,1,3:4)
    hold on
    plot(simMo2S.sim.tsim{1}, simMo2S.sim.states{i}(:,4), 'g')
    plot(simMo2S.sim.tsim{1}, simMo2S.sim.states{i}(:,8), 'b')
    legend('Model 1', 'Model 2')
    ylabel('GFP(A.U.)')
    ylim([0,2500])
    
    subplot(6,1,5)
    stairs([0,duration],[indbest2s(i,1), indbest2s(i,1)])
    ylabel('IPTG')
    
    subplot(6,1,6)
    stairs([0,duration],[indbest2s(i,2), indbest2s(i,2)])
    ylabel('aTc')
    xlabel('time (min)')
    
    saveas(h, ['SimulationsBistability\Plots\','Best2SPicksSimulations_',num2str(i),'.png'])
end




% Best OED picks
for i=1:4
    h=figure;
    subplot(6,1,1:2)
    hold on
    plot(simMoBest.sim.tsim{1}, simMoBest.sim.states{i}(:,3), 'r')
    plot(simMoBest.sim.tsim{1}, simMoBest.sim.states{i}(:,7), 'b')
    legend('Model 1', 'Model 2')
    ylabel('RFP(A.U.)')
    title(['Random ', num2str(i)])
    ylim([0,3500])
    
    subplot(6,1,3:4)
    hold on
    plot(simMoBest.sim.tsim{1}, simMoBest.sim.states{i}(:,4), 'g')
    plot(simMoBest.sim.tsim{1}, simMoBest.sim.states{i}(:,8), 'b')
    legend('Model 1', 'Model 2')
    ylabel('GFP(A.U.)')
    ylim([0,2500])
    
    subplot(6,1,5)
    stairs([0:180:duration],[indbest(r(i),1), indbest(r(i),:)])
    ylabel('IPTG')
    
    subplot(6,1,6)
    stairs([0:180:duration],[indbest(r(i),1), indbest(r(i),:)])
    ylabel('aTc')
    xlabel('time (min)')
    
    saveas(h, ['SimulationsBistability\Plots\','Best2SPicksSimulations_',num2str(i),'.png'])
end


%% Compute Distance Metrics

% Average
RRan = zeros(1,12);
for i=1:12
    [ff,~,~] = OEDMSCostAver_ComparisonBF(simMoRan.sim.states{i});
    RRan(1,i) = ff;
end

R2S = zeros(1,8);
for i=1:8
    [ff,~,~] = OEDMSCostAver_ComparisonBF(simMo2S.sim.states{i});
    R2S(1,i) = ff;
end

RBt = zeros(1,4);
for i=1:4
    [ff,~,~] = OEDMSCostAver_ComparisonBF(simMoBest.sim.states{i});
    RBt(1,i) = ff;
end

all = abs([RRan, R2S, RBt]);

figure
X = categorical({'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
X = reordercats(X,{'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
bar(X, all)




% Multiplicative
RRan = zeros(1,12);
for i=1:12
    [ff,~,~] = OEDMSCostMultiplic_ComparisonBF(simMoRan.sim.states{i});
    RRan(1,i) = ff;
end

R2S = zeros(1,8);
for i=1:8
    [ff,~,~] = OEDMSCostMultiplic_ComparisonBF(simMo2S.sim.states{i});
    R2S(1,i) = ff;
end

RBt = zeros(1,4);
for i=1:4
    [ff,~,~] = OEDMSCostMultiplic_ComparisonBF(simMoBest.sim.states{i});
    RBt(1,i) = ff;
end

all = abs([RRan, R2S, RBt]);

figure
X = categorical({'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
X = reordercats(X,{'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
bar(X, all)


% RFP
RRan = zeros(1,12);
for i=1:12
    [ff,~,~] = OEDMSCostRFP_ComparisonBF(simMoRan.sim.states{i});
    RRan(1,i) = ff;
end

R2S = zeros(1,8);
for i=1:8
    [ff,~,~] = OEDMSCostRFP_ComparisonBF(simMo2S.sim.states{i});
    R2S(1,i) = ff;
end

RBt = zeros(1,4);
for i=1:4
    [ff,~,~] = OEDMSCostRFP_ComparisonBF(simMoBest.sim.states{i});
    RBt(1,i) = ff;
end

all = abs([RRan, R2S, RBt]);

figure
X = categorical({'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
X = reordercats(X,{'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
bar(X, all)



% GFP
RRan = zeros(1,12);
for i=1:12
    [ff,~,~] = OEDMSCostGFP_ComparisonBF(simMoRan.sim.states{i});
    RRan(1,i) = ff;
end

R2S = zeros(1,8);
for i=1:8
    [ff,~,~] = OEDMSCostGFP_ComparisonBF(simMo2S.sim.states{i});
    R2S(1,i) = ff;
end

RBt = zeros(1,4);
for i=1:4
    [ff,~,~] = OEDMSCostGFP_ComparisonBF(simMoBest.sim.states{i});
    RBt(1,i) = ff;
end

all = abs([RRan, R2S, RBt]);

figure
X = categorical({'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
X = reordercats(X,{'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
bar(X, all)


%% Independent Plots

%% RFP
% RFP
RRan = zeros(1,12);
for i=1:12
    [ff,~,~] = OEDMSCostRFP_ComparisonBF(simMoRan.sim.states{i});
    RRan(1,i) = ff;
end

R2S = zeros(1,2);
for i=1:8
    [ff,~,~] = OEDMSCostRFP_ComparisonBF(simMo2S.sim.states{i});
    R2S(1,i) = ff;
end

RBt = zeros(1,4);
for i=1:4
    [ff,~,~] = OEDMSCostRFP_ComparisonBF(simMoBest.sim.states{i});
    RBt(1,i) = ff;
end

all = abs([RRan, R2S, RBt]);

figure
X = categorical({'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
X = reordercats(X,{'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
bar(X, all)

allRFP = abs([RRan, R2S(1:2), RBt(1)]);
h = figure;
hold on
X = categorical({'GFP1','RFP1','RFP2','Both Bistable','M1Bi1','M1Bi2','M3Bi1','M3Bi2','Disc1','Disc2','Disc3','Disc4',...
    'OED2S1','OED2S2','BestOED'});
X = reordercats(X,{'GFP1','RFP1','RFP2','Both Bistable','M1Bi1','M1Bi2','M3Bi1','M3Bi2','Disc1','Disc2','Disc3','Disc4',...
    'OED2S1','OED2S2','BestOED'});
bar(X,allRFP, 'c')
indm = find(max(allRFP)==allRFP);
bar(X(indm), allRFP(indm), 'r')
ylabel('abs(Cost Function)')
saveas(h, ['SimulationsBistability\Plots\','CostFunctionValuesFrequentist_','RFP.png'])




%% GFP
% RFP
RRan = zeros(1,12);
for i=1:12
    [ff,~,~] = OEDMSCostGFP_ComparisonBF(simMoRan.sim.states{i});
    RRan(1,i) = ff;
end

R2S = zeros(1,2);
for i=1:8
    [ff,~,~] = OEDMSCostGFP_ComparisonBF(simMo2S.sim.states{i});
    R2S(1,i) = ff;
end

RBt = zeros(1,4);
for i=1:4
    [ff,~,~] = OEDMSCostGFP_ComparisonBF(simMoBest.sim.states{i});
    RBt(1,i) = ff;
end

all = abs([RRan, R2S, RBt]);

figure
X = categorical({'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
X = reordercats(X,{'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
bar(X, all)

allGFP = abs([RRan, R2S(3:4), RBt(2)]);
h = figure;
hold on
X = categorical({'GFP1','RFP1','RFP2','Both Bistable','M1Bi1','M1Bi2','M3Bi1','M3Bi2','Disc1','Disc2','Disc3','Disc4',...
    'OED2S1','OED2S2','BestOED'});
X = reordercats(X,{'GFP1','RFP1','RFP2','Both Bistable','M1Bi1','M1Bi2','M3Bi1','M3Bi2','Disc1','Disc2','Disc3','Disc4',...
    'OED2S1','OED2S2','BestOED'});
bar(X,allGFP, 'c')
indm = find(max(allGFP)==allGFP);
bar(X(indm), allGFP(indm), 'g')
ylabel('abs(Cost Function)')
saveas(h, ['SimulationsBistability\Plots\','CostFunctionValuesFrequentist_','GFP.png'])



%% Average
% RFP
RRan = zeros(1,12);
for i=1:12
    [ff,~,~] = OEDMSCostAver_ComparisonBF(simMoRan.sim.states{i});
    RRan(1,i) = ff;
end

R2S = zeros(1,2);
for i=1:8
    [ff,~,~] = OEDMSCostAver_ComparisonBF(simMo2S.sim.states{i});
    R2S(1,i) = ff;
end

RBt = zeros(1,4);
for i=1:4
    [ff,~,~] = OEDMSCostAver_ComparisonBF(simMoBest.sim.states{i});
    RBt(1,i) = ff;
end

all = abs([RRan, R2S, RBt]);

figure
X = categorical({'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
X = reordercats(X,{'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
bar(X, all)

allAve = abs([RRan, R2S(5:6), RBt(3)]);
h = figure;
hold on
X = categorical({'GFP1','RFP1','RFP2','Both Bistable','M1Bi1','M1Bi2','M3Bi1','M3Bi2','Disc1','Disc2','Disc3','Disc4',...
    'OED2S1','OED2S2','BestOED'});
X = reordercats(X,{'GFP1','RFP1','RFP2','Both Bistable','M1Bi1','M1Bi2','M3Bi1','M3Bi2','Disc1','Disc2','Disc3','Disc4',...
    'OED2S1','OED2S2','BestOED'});
bar(X,allAve, 'c')
indm = find(max(allAve)==allAve);
bar(X(indm), allAve(indm), 'y')
ylabel('abs(Cost Function)')
saveas(h, ['SimulationsBistability\Plots\','CostFunctionValuesFrequentist_','Average.png'])




%% Multiplicative
% RFP
RRan = zeros(1,12);
for i=1:12
    [ff,~,~] = OEDMSCostMultiplic_ComparisonBF([simMod1.sim.states{7}, simMod3.sim.states{7}]);
    RRan(1,i) = ff;
end

R2S = zeros(1,2);
for i=1:8
    [ff,~,~] = OEDMSCostMultiplic_ComparisonBF(simMo2S.sim.states{i});
    R2S(1,i) = ff;
end

RBt = zeros(1,4);
for i=1:4
    [ff,~,~] = OEDMSCostMultiplic_ComparisonBF(simMoBest.sim.states{i});
    RBt(1,i) = ff;
end

all = abs([RRan, R2S, RBt]);

figure
X = categorical({'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
X = reordercats(X,{'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12',...
    '2S1','2S2','2S3','2S4','2S5','2S6','2S7','2S8','B1','B2','B3','B4'});
bar(X, all)

allMul = abs([RRan, R2S(7:8), RBt(4)]);
h = figure;
hold on
X = categorical({'GFP1','RFP1','RFP2','Both Bistable','M1Bi1','M1Bi2','M3Bi1','M3Bi2','Disc1','Disc2','Disc3','Disc4',...
    'OED2S1','OED2S2','BestOED'});
X = reordercats(X,{'GFP1','RFP1','RFP2','Both Bistable','M1Bi1','M1Bi2','M3Bi1','M3Bi2','Disc1','Disc2','Disc3','Disc4',...
    'OED2S1','OED2S2','BestOED'});
bar(X,allMul, 'c')
indm = find(max(allMul)==allMul);
bar(X(indm), allMul(indm), 'b')
ylabel('abs(Cost Function)')
saveas(h, ['SimulationsBistability\Plots\','CostFunctionValuesFrequentist_','Multiplicative.png'])









