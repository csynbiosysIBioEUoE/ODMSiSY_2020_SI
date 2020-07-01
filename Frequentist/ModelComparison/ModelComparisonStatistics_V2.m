%% Model Comparison scripts (frequentist)
% The script generates a series of statistics for for the three models in order
% to assess the model that better descrives the data fitted (model
% selection). It is a simplified version of the script
% ModelComparisonStatistics.m so no inputs are required and the 3 models
% are assessed at once. 

% Computes different model comparison statistics on Training dataset

% References -->
%     - Methods and Criteria for Model Selection - (http://www.stat.cmu.edu/tr/tr759/tr759.pdf)
%     - http://statweb.stanford.edu/~jtaylo/courses/stats203/notes/selection.pdf



function [] = ModelComparisonStatistics_V2()

%% Simulate Models

Dat = load('F:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\LugagneData\BayesSelectDataLugagne.mat');
Mod1 = ToggleSwitch_load_model_M1;
Mod2 = ToggleSwitch_load_model_M2;
Mod3 = ToggleSwitch_load_model_M3;

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



clear newExps;
newExps.n_exp = length(Dat.Data.exp_data);                                         % Number of experiments 
for i=1:length(Dat.Data.exp_data)
    
    y0 = M1_Compute_SteadyState_OverNight(inputs,inputs.model.par,Dat.Data.exp_data{i}(end,:),[Dat.Data.Initial_IPTG{i} Dat.Data.Initial_aTc{i}]+1e-7);
    
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
inputs.ivpsol.senssolver='fdsens5';
inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;

inputs.plotd.plotlevel='noplot';

inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = strcat('sim-Comp');

% simDa = AMIGO_SData(inputs);
AMIGO_Prep(inputs);

simMod1 = AMIGO_SModel(inputs);
save(['AMIGOSsimModel','1','.mat'], 'simMod1')

% AMIGO Stats
resu.model = Mod1;
resu.sim = simMod1.sim;
resu.sim.exp_data = Dat;
statsM1 = AMIGO_PEPostAnalysisCustom(inputs,resu);
save(['AMIGOStatsModel','1','.mat'], 'statsM1')



% Model 2
clear model;
clear exps;
clear best_global_theta;
clear pe_results;
clear pe_inputs;

% Specify folder name and short_name
results_folder = strcat('MTS',datestr(now,'yyyy-mm-dd-HHMMSS'));
short_name     = strcat('MTS_Statistics');

clear inputs;
inputs.model = Mod2;
inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = 'initial_setup';



clear newExps;
newExps.n_exp = length(Dat.Data.exp_data);                                         % Number of experiments 
for i=1:length(Dat.Data.exp_data)
    

    y0 = M2_Compute_SteadyState_OverNight(inputs,inputs.model.par,Dat.Data.exp_data{i}(end,:),[Dat.Data.Initial_IPTG{i} Dat.Data.Initial_aTc{i}]+1e-7);
    
    newExps.n_obs{i}=2;                                        % Number of observables per experiment        
    newExps.obs_names{i} = char('LacI_M2','TetR_M2');
    newExps.obs{i} = char('LacI_M2 = L_RFP','TetR_M2 = T_GFP');% Name of the observables 

    newExps.exp_y0{i}=y0;                                      % Initial condition for the experiment    

    newExps.t_f{i}=round(Dat.Data.t_samples{i}(end,1));                                   % Experiment duration
    newExps.n_s{i}=Dat.Data.n_samples{i}(1); 
    newExps.t_s{i} = round(Dat.Data.t_samples{i}(:,1))'; % Number of sampling times
%     newExps.t_s{i}=round(Dat.Data.t_samples{i}(:,1))';                              % Times of samples

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
inputs.ivpsol.senssolver='fdsens5';
inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;

inputs.plotd.plotlevel='noplot';

inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = strcat('sim-Comp');

% simDa = AMIGO_SData(inputs);
AMIGO_Prep(inputs);

simMod2 = AMIGO_SModel(inputs);
save(['AMIGOSsimModel','2','.mat'], 'simMod2')


% AMIGO Stats
resu.model = Mod2;
resu.sim = simMod2.sim;
resu.sim.exp_data = Dat;
statsM2 = AMIGO_PEPostAnalysisCustom(inputs,resu);
save(['AMIGOStatsModel','2','.mat'], 'statsM2')




% Model 3
clear model;
clear exps;
clear best_global_theta;
clear pe_results;
clear pe_inputs;

% Specify folder name and short_name
results_folder = strcat('MTS',datestr(now,'yyyy-mm-dd-HHMMSS'));
short_name     = strcat('MTS_Statistics');

clear inputs;
inputs.model = Mod3;
inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = 'initial_setup';



clear newExps;
newExps.n_exp = length(Dat.Data.exp_data);                                         % Number of experiments 
for i=1:length(Dat.Data.exp_data)
    
    y0 = M3_Compute_SteadyState_OverNight(inputs,inputs.model.par,Dat.Data.exp_data{i}(end,:),[Dat.Data.Initial_IPTG{i} Dat.Data.Initial_aTc{i}]+1e-7);
    
    newExps.n_obs{i}=2;                                        % Number of observables per experiment        
    newExps.obs_names{i} = char('LacI_M3','TetR_M3');
    newExps.obs{i} = char('LacI_M3 = L_RFP','TetR_M3 = T_GFP');% Name of the observables 

    newExps.exp_y0{i}=y0;                                      % Initial condition for the experiment    

    newExps.t_f{i}=round(Dat.Data.t_samples{i}(end,1));                                   % Experiment duration
    newExps.n_s{i}=Dat.Data.n_samples{i}(1); 
    newExps.t_s{i} = round(Dat.Data.t_samples{i}(:,1))'; % Number of sampling times
%     newExps.t_s{i}=round(Dat.Data.t_samples{i}(:,1))';                              % Times of samples

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
inputs.ivpsol.senssolver='fdsens5';
inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;

inputs.plotd.plotlevel='noplot';

inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = strcat('sim-Comp');

% simDa = AMIGO_SData(inputs);
AMIGO_Prep(inputs);

simMod3 = AMIGO_SModel(inputs);
save(['AMIGOSsimModel','3','.mat'], 'simMod3')


% AMIGO Stats
resu.model = Mod3;
resu.sim = simMod3.sim;
resu.sim.exp_data = Dat;
statsM3 = AMIGO_PEPostAnalysisCustom(inputs,resu);
save(['AMIGOStatsModel','3','.mat'], 'statsM3')


%% Display Statistics

% Extract each one of the statistics of interest and print the results,
% including value for each model and which model would be the best
% according to the metric

[a11,a12,a13,a14,a15,a16,a17,a18,a19,a110,AIC_M1]=statsM1.AkaikeInformationCrit;
[a21,a22,a23,a24,a25,a26,a27,a28,a29,a210,AIC_M2]=statsM2.AkaikeInformationCrit;
[a31,a32,a33,a34,a35,a36,a37,a38,a39,a310,AIC_M3]=statsM3.AkaikeInformationCrit;

sm1 = [a11.AIC,a12.AIC,a13.AIC,a14.AIC,a15.AIC,a16.AIC,a17.AIC,a18.AIC,a19.AIC,a110.AIC];
sm2 = [a21.AIC,a22.AIC,a23.AIC,a24.AIC,a25.AIC,a26.AIC,a27.AIC,a28.AIC,a29.AIC,a210.AIC];
sm3 = [a31.AIC,a32.AIC,a33.AIC,a34.AIC,a35.AIC,a36.AIC,a37.AIC,a38.AIC,a39.AIC,a310.AIC];

y = zeros(length(sm1), 3);
for i=1:length(sm1)
    y(i,1) = sm1(i);
    y(i,2) = sm2(i);
    y(i,3) = sm3(i);
end
h=figure;
bar(y)
title('Single Experiment AIC')
xlabel('Training Experiment')
ylabel('AIC')
legend('M1', 'M2', 'M3')
set(gca, 'YScale', 'log')
saveas(h, ['Plots\AIC_TrainingSet_','Final','.png'])

disp(' ')
disp('*************************  AIC *************************')
disp(' ')
disp(['AIC for Model ', '1', ' is ', num2str(AIC_M1.AIC)])
disp(['AIC for Model ', '2', ' is ', num2str(AIC_M2.AIC)])
disp(['AIC for Model ', '3', ' is ', num2str(AIC_M3.AIC)])
disp(' ')

if AIC_M1.AIC == min([AIC_M1.AIC,AIC_M2.AIC,AIC_M3.AIC])
    disp(['Model ', '1', ' is better'])
elseif AIC_M2.AIC == min([AIC_M1.AIC,AIC_M2.AIC,AIC_M3.AIC])
    disp(['Model ', '2', ' is better'])
elseif AIC_M3.AIC == min([AIC_M1.AIC,AIC_M2.AIC,AIC_M3.AIC])
    disp(['Model ', '3', ' is better'])    
end
disp(' ')
disp('********************************************************')
disp(' ')



[b11,b12,b13,b14,b15,b16,b17,b18,b19,b110,BIC_M1]=statsM1.BayesianInformationCrit;
[b21,b22,b23,b24,b25,b26,b27,b28,b29,b210,BIC_M2]=statsM2.BayesianInformationCrit;
[b31,b32,b33,b34,b35,b36,b37,b38,b39,b310,BIC_M3]=statsM3.BayesianInformationCrit;

sm1 = [b11,b12,b13,b14,b15,b16,b17,b18,b19,b110];
sm2 = [b21,b22,b23,b24,b25,b26,b27,b28,b29,b210];
sm3 = [b31,b32,b33,b34,b35,b36,b37,b38,b39,b310];

y = zeros(length(sm1), 3);
for i=1:length(sm1)
    y(i,1) = sm1(i);
    y(i,2) = sm2(i);
    y(i,3) = sm3(i);
end
h=figure;
bar(y)
title('Single Experiment BIC')
xlabel('Training Experiment')
ylabel('BIC')
legend('M1', 'M2', 'M3')
set(gca, 'YScale', 'log')
saveas(h, ['Plots\BIC_TrainingSet_','Final','.png'])

disp(' ')
disp('*************************  BIC *************************')
disp(' ')
disp(['BIC for Model ', '1', ' is ', num2str(BIC_M1)])
disp(['BIC for Model ', '2', ' is ', num2str(BIC_M2)])
disp(['BIC for Model ', '3', ' is ', num2str(BIC_M3)])
disp(' ')

if BIC_M1 == min([BIC_M1,BIC_M2,BIC_M3])
    disp(['Model ', '1', ' is better'])
elseif BIC_M2 == min([BIC_M1,BIC_M2,BIC_M3])
    disp(['Model ', '2', ' is better'])
elseif BIC_M3 == min([BIC_M1,BIC_M2,BIC_M3])
    disp(['Model ', '3', ' is better'])    
end
disp(' ')
disp('********************************************************')
disp(' ')




[b11,b12,b13,b14,b15,b16,b17,b18,b19,b110,R2_M1]=statsM1.R2;
[b21,b22,b23,b24,b25,b26,b27,b28,b29,b210,R2_M2]=statsM2.R2;
[b31,b32,b33,b34,b35,b36,b37,b38,b39,b310,R2_M3]=statsM3.R2;

sm1 = [b11,b12,b13,b14,b15,b16,b17,b18,b19,b110];
sm2 = [b21,b22,b23,b24,b25,b26,b27,b28,b29,b210];
sm3 = [b31,b32,b33,b34,b35,b36,b37,b38,b39,b310];

y = zeros(length(sm1), 3);
for i=1:length(sm1)
    y(i,1) = sm1(i);
    y(i,2) = sm2(i);
    y(i,3) = sm3(i);
end
h=figure;
bar(abs(y))
title('Single Experiment R^2')
xlabel('Training Experiment')
ylabel('Absolute R^2')
legend('M1', 'M2', 'M3')
set(gca, 'YScale', 'log')
saveas(h, ['Plots\R2_TrainingSet_','Final','.png'])

disp(' ')
disp('*************************  R^2 *************************')
disp(' ')
disp(['R^2 for Model ', '1', ' is ', num2str(R2_M1)])
disp(['R^2 for Model ', '2', ' is ', num2str(R2_M2)])
disp(['R^2 for Model ', '3', ' is ', num2str(R2_M3)])
disp(' ')

if R2_M1 == max([R2_M1,R2_M2,R2_M3])
    disp(['Model ', '1', ' is better'])
elseif R2_M2 == max([R2_M1,R2_M2,R2_M3])
    disp(['Model ', '2', ' is better'])
elseif R2_M3 == max([R2_M1,R2_M2,R2_M3])
    disp(['Model ', '3', ' is better'])    
end
disp(' ')
disp('********************************************************')
disp(' ')




[b11,b12,b13,b14,b15,b16,b17,b18,b19,b110,RMSE_M1]=statsM1.RMSE;
[b21,b22,b23,b24,b25,b26,b27,b28,b29,b210,RMSE_M2]=statsM2.RMSE;
[b31,b32,b33,b34,b35,b36,b37,b38,b39,b310,RMSE_M3]=statsM3.RMSE;


sm1 = [sum(b11),sum(b12),sum(b13),sum(b14),sum(b15),sum(b16),sum(b17),sum(b18),sum(b19),sum(b110)];
sm2 = [sum(b21),sum(b22),sum(b23),sum(b24),sum(b25),sum(b26),sum(b27),sum(b28),sum(b29),sum(b210)];
sm3 = [sum(b31),sum(b32),sum(b33),sum(b34),sum(b35),sum(b36),sum(b37),sum(b38),sum(b39),sum(b310)];

y = zeros(length(sm1), 3);
for i=1:length(sm1)
    y(i,1) = sm1(i);
    y(i,2) = sm2(i);
    y(i,3) = sm3(i);
end
h=figure;
bar(abs(y))
title('Single Experiment RMSE')
xlabel('Training Experiment')
ylabel('RMSE')
legend('M1', 'M2', 'M3')
saveas(h, ['Plots\RMSE_TrainingSet_','Final','.png'])

disp(' ')
disp('*************************  RMSE *************************')
disp(' ')
disp(['RMSE for Model ', '1', ' is ', num2str(RMSE_M1)])
disp(['RMSE for Model ', '2', ' is ', num2str(RMSE_M2)])
disp(['RMSE for Model ', '3', ' is ', num2str(RMSE_M3)])
disp(' ')

if RMSE_M1 == min([RMSE_M1,RMSE_M2,RMSE_M3])
    disp(['Model ', '1', ' is better'])
elseif RMSE_M2 == min([RMSE_M1,RMSE_M2,RMSE_M3])
    disp(['Model ', '2', ' is better'])
elseif RMSE_M3 == min([RMSE_M1,RMSE_M2,RMSE_M3])
    disp(['Model ', '3', ' is better'])    
end
disp(' ')
disp('********************************************************')
disp(' ')


%% Mallow's Cp

SSEr_M1 = cell(1,length(Dat.Data.exp_data));
SSEg_M1 = cell(1,length(Dat.Data.exp_data));
SSEr_M2 = cell(1,length(Dat.Data.exp_data));
SSEg_M2 = cell(1,length(Dat.Data.exp_data));
SSEr_M3 = cell(1,length(Dat.Data.exp_data));
SSEg_M3 = cell(1,length(Dat.Data.exp_data));

N = [];
for exp=1:length(Dat.Data.exp_data)
    
    for samp=1:length(Dat.Data.exp_data{exp}(:,1))
        SSEr_M1{1,exp}(samp) = ((simMod1.sim.states{exp}(samp,3) - Dat.Data.exp_data{exp}(samp,1))^2)/Dat.Data.standard_dev{exp}(samp,1)^2;
        SSEg_M1{1,exp}(samp) = ((simMod1.sim.states{exp}(samp,4) - Dat.Data.exp_data{exp}(samp,2))^2)/Dat.Data.standard_dev{exp}(samp,2)^2;
        
        SSEr_M2{1,exp}(samp) = ((simMod2.sim.states{exp}(samp,3) - Dat.Data.exp_data{exp}(samp,1))^2)/Dat.Data.standard_dev{exp}(samp,1)^2;
        SSEg_M2{1,exp}(samp) = ((simMod2.sim.states{exp}(samp,4) - Dat.Data.exp_data{exp}(samp,2))^2)/Dat.Data.standard_dev{exp}(samp,2)^2;
        
        SSEr_M3{1,exp}(samp) = ((simMod3.sim.states{exp}(samp,3) - Dat.Data.exp_data{exp}(samp,1))^2)/Dat.Data.standard_dev{exp}(samp,1)^2;
        SSEg_M3{1,exp}(samp) = ((simMod3.sim.states{exp}(samp,4) - Dat.Data.exp_data{exp}(samp,2))^2)/Dat.Data.standard_dev{exp}(samp,2)^2;
        
    end
    N = [N, Dat.Data.n_samples{exp}(1)];
end

for i=1:length(Dat.Data.exp_data)
    SSE_M1(i) = sum(SSEr_M1{i})+sum(SSEg_M1{i});
    SSE_M2(i) = sum(SSEr_M2{i})+sum(SSEg_M2{i});
    SSE_M3(i) = sum(SSEr_M3{i})+sum(SSEg_M3{i});
    
    sm1(i) = SSE_M1(i) - (N(i))*2 + 2*Mod1.n_par;
    sm2(i) = SSE_M2(i) - (N(i))*2 + 2*Mod1.n_par;
    sm3(i) = SSE_M3(i) - (N(i))*2 + 2*Mod1.n_par;
end

y = zeros(length(sm1), 3);
for i=1:length(sm1)
    y(i,1) = sm1(i);
    y(i,2) = sm2(i);
    y(i,3) = sm3(i);
end
h=figure;
bar(abs(y))
title('Single Experiment Mallow''s Cp')
xlabel('Training Experiment')
ylabel('Mallow''s Cp')
legend('M1', 'M2', 'M3')
set(gca, 'YScale', 'log')
saveas(h, ['Plots\MCP_TrainingSet_','Final','.png'])

SSE_M1 = sum(SSE_M1);
SSE_M2 = sum(SSE_M2);
SSE_M3 = sum(SSE_M3);

Cp_M1 = SSE_M1 - sum(N)*2 + 2*Mod1.n_par;
Cp_M2 = SSE_M2 - sum(N)*2 + 2*Mod2.n_par;
Cp_M3 = SSE_M3 - sum(N)*2 + 2*Mod3.n_par;


disp(' ')
disp('*************************  Mallow''s Cp *************************')
disp(' ')
disp(['Mallow''s Cp for Model ', '1', ' is ', num2str(Cp_M1)])
disp(['Mallow''s Cp for Model ', '2', ' is ', num2str(Cp_M2)])
disp(['Mallow''s Cp for Model ', '3', ' is ', num2str(Cp_M3)])
disp(' ')

if Cp_M1 == min([Cp_M1,Cp_M2,Cp_M3])
    disp(['Model ', '1', ' is better'])
elseif Cp_M2 == min([Cp_M1,Cp_M2,Cp_M3])
    disp(['Model ', '2', ' is better'])
elseif Cp_M3 == min([Cp_M1,Cp_M2,Cp_M3])
    disp(['Model ', '3', ' is better'])    
end
disp(' ')
disp('********************************************************')
disp(' ')


%% SSE

SSEr_M1 = cell(1,length(Dat.Data.exp_data));
SSEg_M1 = cell(1,length(Dat.Data.exp_data));
SSEr_M2 = cell(1,length(Dat.Data.exp_data));
SSEg_M2 = cell(1,length(Dat.Data.exp_data));
SSEr_M3 = cell(1,length(Dat.Data.exp_data));
SSEg_M3 = cell(1,length(Dat.Data.exp_data));

for exp=1:length(Dat.Data.exp_data)
    
    for samp=1:length(Dat.Data.exp_data{exp}(:,1))
        SSEr_M1{1,exp}(samp) = ((simMod1.sim.states{exp}(samp,3) - Dat.Data.exp_data{exp}(samp,1))^2);
        SSEg_M1{1,exp}(samp) = ((simMod1.sim.states{exp}(samp,4) - Dat.Data.exp_data{exp}(samp,2))^2);
        
        SSEr_M2{1,exp}(samp) = ((simMod2.sim.states{exp}(samp,3) - Dat.Data.exp_data{exp}(samp,1))^2);
        SSEg_M2{1,exp}(samp) = ((simMod2.sim.states{exp}(samp,4) - Dat.Data.exp_data{exp}(samp,2))^2);
        
        SSEr_M3{1,exp}(samp) = ((simMod3.sim.states{exp}(samp,3) - Dat.Data.exp_data{exp}(samp,1))^2);
        SSEg_M3{1,exp}(samp) = ((simMod3.sim.states{exp}(samp,4) - Dat.Data.exp_data{exp}(samp,2))^2);
        
    end
end


for i=1:length(Dat.Data.exp_data)
    SSE_M1(i) = sum(SSEr_M1{i})+sum(SSEg_M1{i});
    SSE_M2(i) = sum(SSEr_M2{i})+sum(SSEg_M2{i});
    SSE_M3(i) = sum(SSEr_M3{i})+sum(SSEg_M3{i});
    
    sm1(i) = SSE_M1(i);
    sm2(i) = SSE_M2(i);
    sm3(i) = SSE_M3(i);
end

y = zeros(length(sm1), 3);
for i=1:length(sm1)
    y(i,1) = sm1(i);
    y(i,2) = sm2(i);
    y(i,3) = sm3(i);
end
h=figure;
bar(abs(y))
title('Single Experiment SSE')
xlabel('Training Experiment')
ylabel('SSE')
legend('M1', 'M2', 'M3')
set(gca, 'YScale', 'log')
saveas(h, ['Plots\SSE_TrainingSet_','Final','.png'])

SSE_M1 = sum(SSE_M1);
SSE_M2 = sum(SSE_M2);
SSE_M3 = sum(SSE_M3);

disp(' ')
disp('*************************  SSE  *************************')
disp(' ')
disp(['SSE for Model ', '1', ' is ', num2str(SSE_M1)])
disp(['SSE for Model ', '2', ' is ', num2str(SSE_M2)])
disp(['SSE for Model ', '3', ' is ', num2str(SSE_M3)])
disp(' ')

if SSE_M1 == min([SSE_M1,SSE_M2,SSE_M3])
    disp(['Model ', '1', ' is better'])
elseif SSE_M2 == min([SSE_M1,SSE_M2,SSE_M3])
    disp(['Model ', '2', ' is better'])
elseif SSE_M3 == min([SSE_M1,SSE_M2,SSE_M3])
    disp(['Model ', '3', ' is better'])    
end
disp(' ')
disp('********************************************************')
disp(' ')


%% nRMSE

nRMSEr_M1 = cell(1,length(Dat.Data.exp_data));
nRMSEg_M1 = cell(1,length(Dat.Data.exp_data));
nRMSEr_M2 = cell(1,length(Dat.Data.exp_data));
nRMSEg_M2 = cell(1,length(Dat.Data.exp_data));
nRMSEr_M3 = cell(1,length(Dat.Data.exp_data));
nRMSEg_M3 = cell(1,length(Dat.Data.exp_data));

N = zeros(1,length(Dat.Data.n_samples));
for exp=1:length(Dat.Data.exp_data)
    
    for samp=1:length(Dat.Data.exp_data{exp}(:,1))
        nRMSEr_M1{1,exp}(samp) = ((simMod1.sim.states{exp}(samp,3) - Dat.Data.exp_data{exp}(samp,1))^2)/Dat.Data.standard_dev{exp}(samp,1)^2;
        nRMSEg_M1{1,exp}(samp) = ((simMod1.sim.states{exp}(samp,4) - Dat.Data.exp_data{exp}(samp,2))^2)/Dat.Data.standard_dev{exp}(samp,2)^2;
        
        nRMSEr_M2{1,exp}(samp) = ((simMod2.sim.states{exp}(samp,3) - Dat.Data.exp_data{exp}(samp,1))^2)/Dat.Data.standard_dev{exp}(samp,1)^2;
        nRMSEg_M2{1,exp}(samp) = ((simMod2.sim.states{exp}(samp,4) - Dat.Data.exp_data{exp}(samp,2))^2)/Dat.Data.standard_dev{exp}(samp,2)^2;
        
        nRMSEr_M3{1,exp}(samp) = ((simMod3.sim.states{exp}(samp,3) - Dat.Data.exp_data{exp}(samp,1))^2)/Dat.Data.standard_dev{exp}(samp,1)^2;
        nRMSEg_M3{1,exp}(samp) = ((simMod3.sim.states{exp}(samp,4) - Dat.Data.exp_data{exp}(samp,2))^2)/Dat.Data.standard_dev{exp}(samp,2)^2;
    end
    N(1, exp) = Dat.Data.n_samples{exp}(1);
    nRMSEr_M1b(1,exp) = sqrt(sum(nRMSEr_M1{1,exp})/N(1, exp));
    nRMSEg_M1b(1,exp) = sqrt(sum(nRMSEg_M1{1,exp})/N(1, exp));
    
    nRMSEr_M2b(1,exp) = sqrt(sum(nRMSEr_M2{1,exp})/N(1, exp));
    nRMSEg_M2b(1,exp) = sqrt(sum(nRMSEg_M2{1,exp})/N(1, exp));
    
    nRMSEr_M3b(1,exp) = sqrt(sum(nRMSEr_M3{1,exp})/N(1, exp));
    nRMSEg_M3b(1,exp) = sqrt(sum(nRMSEg_M3{1,exp})/N(1, exp));
end

y = zeros(length(sm1), 3);
for i=1:length(sm1)
    y(i,1) = nRMSEr_M1b(i)+nRMSEg_M1b(i);
    y(i,2) = nRMSEr_M2b(i)+nRMSEg_M2b(i);
    y(i,3) = nRMSEr_M3b(i)+nRMSEg_M3b(i);
end
h=figure;
bar(abs(y))
title('Single Experiment nRMSE')
xlabel('Training Experiment')
ylabel('nRMSE')
legend('M1', 'M2', 'M3')
set(gca, 'YScale', 'log')
saveas(h, ['Plots\nRMSE_TrainingSet_','Final','.png'])

nRMSE_M1 = sum(nRMSEr_M1b)+sum(nRMSEg_M1b);
nRMSE_M2 = sum(nRMSEr_M2b)+sum(nRMSEg_M2b);
nRMSE_M3 = sum(nRMSEr_M3b)+sum(nRMSEg_M3b);


disp(' ')
disp('*************************  nRMSE  *************************')
disp(' ')
disp(['nRMSE for Model ', '1', ' is ', num2str(nRMSE_M1)])
disp(['nRMSE for Model ', '2', ' is ', num2str(nRMSE_M2)])
disp(['nRMSE for Model ', '3', ' is ', num2str(nRMSE_M3)])
disp(' ')

if nRMSE_M1 == min([nRMSE_M1,nRMSE_M2,nRMSE_M3])
    disp(['Model ', '1', ' is better'])
elseif nRMSE_M2 == min([nRMSE_M1,nRMSE_M2,nRMSE_M3])
    disp(['Model ', '2', ' is better'])
elseif nRMSE_M3 == min([nRMSE_M1,nRMSE_M2,nRMSE_M3])
    disp(['Model ', '3', ' is better'])    
end
disp(' ')
disp('********************************************************')
disp(' ')
end












