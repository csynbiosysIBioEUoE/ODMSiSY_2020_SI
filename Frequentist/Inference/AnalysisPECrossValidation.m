%% Script to extract the best PE run
% Extract the iteration yielding the minimum SSE over the test set. 

f1 = 'E:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\ParameterEstimation\'; % Main directory
Folder = [f1, 'PE_Final\PE_M3_CV_Rep2']; % Directory where the results are at

filePattern = fullfile(Folder, strcat('PE_Model3_CrosValidation_Rep2-*PEMultiTest.mat')); % Common string for PE results
Files = dir(filePattern); 

filePattern2 = fullfile(Folder, strcat('PE_Model3_CrosValidation_Rep2-*sim.mat')); % Common string for simulation and SSE results from each run
Files2 = dir(filePattern2); 

tags = 'M3CVRep2'; % Tag added in the saved files
%% Plot convergence Curves

figure
hold on
for i=1:length(Files)
    x = load([Folder, '\', Files(i).name]);
    stairs(x.pe_results.nlpsol.conv_curve(:,1), x.pe_results.nlpsol.conv_curve(:,2))
end
xlabel('Iteration')
ylabel('CFV')
% set(gca, 'YScale', 'log')


%% Create a structure with the SSE values over the test set for each PE iteration
SSE_Mat = cell(length(Files2),1);
SSE_vect = zeros(length(Files2),1);
for i=1:length(Files2)
    FileName = Files2(i).name;
    x = load([Folder, '\', FileName]);
    SSE_Mat{i,:} = x.SSE;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp = [];
    for j=1:length(x.SSE)
        if j~=2 && j~=13 && j~=15 && j~= 16 && j~=4 && j~=5 && j~=8
            temp = [temp; x.SSE(j,:)];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SSE_vect(i,1) = sum(temp(:), 'omitnan');
%     SSE_vect(i,1) = sum(x.SSE(:), 'omitnan');
end

%% Plot the values of the SSE, summed over the test set, for each iteration
figure; 
plot([1:1:length(Files2)]',SSE_vect);hold on; 
plot(find(SSE_vect == min(SSE_vect)),SSE_vect(find(SSE_vect == min(SSE_vect))),'*r')
xlabel('iteration index')
ylabel('SSE on the test set')

%% Load the file corresponding to the minimum SSE

load([Folder, '\', Files2(find(SSE_vect == min(SSE_vect))).name]);
exps_indexTest = 1:16;

if ~isfolder([Folder,'\Plots'])
    mkdir([Folder,'\Plots']);
end

% Plot simulations over the test set for the best parameter estimate
for i=1:sim_exps.n_exp
    h=figure;
    subplot(2,1,1)
    errorbar(sim_exps.t_s{1,i}/60,sim_exps.exp_data{1,i}(:,1),sim_exps.error_data{1,i}(:,1),'ok'); hold on; 
    plot(sim_results.sim.tsim{1,i}/60,sim_results.sim.obs{1,i}(:,1),'r','LineWidth',2)
    title(int2str(exps_indexTest(i)))
    legend('experimental data', 'best fit')
    xlabel('Time (hours)')
    ylabel('RFP (AU)')
    
    subplot(2,1,2)
    errorbar(sim_exps.t_s{1,i}/60,sim_exps.exp_data{1,i}(:,2),sim_exps.error_data{1,i}(:,2),'ok'); hold on; 
    plot(sim_results.sim.tsim{1,i}/60,sim_results.sim.obs{1,i}(:,2),'g','LineWidth',2)
    title(int2str(exps_indexTest(i)))
    legend('experimental data', 'best fit')
    xlabel('Time (hours)')
    ylabel('RFP (AU)')
    
%     saveas(h, [Folder,'\Plots\BestFitExper_ValidationSet_',num2str(i),'_',tags,'.png'])
end
%%
disp('The best estimate is')
best_global_theta




%% Plot simulations agains training data 

bestx = load([Folder, '\', Files(find(SSE_vect == min(SSE_vect))).name]);
for i=1:bestx.pe_inputs.exps.n_exp
    h = figure;
    subplot(2,1,1)
    hold on 
    errorbar(bestx.pe_inputs.exps.t_s{i}, bestx.pe_inputs.exps.exp_data{i}(:,1), bestx.pe_inputs.exps.error_data{i}(:,1), 'r')
%     plot(bestx.pe_results.sim.tsim{i}, simMo.sim.states{i}(:,3), 'b')% For UpdatedY0
    plot(bestx.pe_results.sim.tsim{i}, bestx.pe_results.sim.states{i}(:,3), 'b')
    ylabel('RFP')
    xlabel('time(min)')

    subplot(2,1,2)
    hold on 
    errorbar(bestx.pe_inputs.exps.t_s{i}, bestx.pe_inputs.exps.exp_data{i}(:,2), bestx.pe_inputs.exps.error_data{i}(:,2), 'g')
%     plot(bestx.pe_results.sim.tsim{i}, simMo.sim.states{i}(:,4), 'b')% For UpdatedY0
    plot(bestx.pe_results.sim.tsim{i}, bestx.pe_results.sim.states{i}(:,4), 'b')
    ylabel('GFP')
    xlabel('time(min)')
    
%     saveas(h, [Folder,'\Plots\BestFitExper_TrainingSet_',num2str(i),'_',tags,'.png'])
end

%% Plot Distribution Of parameters compared to Bayesian Results
mdl=3;
% Read Posterior
d2 = 'F:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\';
posterior = readtable([d2, 'PosteriorDistributionsBayes\draws_ALL_Model',num2str(mdl),'.csv']);
% Extract data
thetaNames = posterior.Properties.VariableNames;
thetaMatrix = posterior.Variables;

% Frequentist
m2 = bestx.pe_results.fit.thetabest';
c2 = bestx.pe_results.fit.g_var_cov_mat;
sd1 = zeros(1,length(m2));
for i=1:length(m2)
    sd1(1,i) = sqrt(c2(i,i));
end

for th=1:length(m2)
    h = figure;
    hold on 
    ksdensity(thetaMatrix(:,th))
    top = (m2(th)+5*sd1(th));
    bot = (m2(th)-5*sd1(th));
    x = bot:((top-bot)/200):top;
    y = normpdf(x,m2(th),sd1(th));
    plot(x,y)
    ylabel('Density')
    thetanam = erase(thetaNames{th},"_raw");
    xlabel(strrep(thetanam, '_', '\_'))
    legend('Bayesian', 'Frequentist')
    hold off
%     saveas(h, [Folder,'\Plots\KernelDensity_',num2str(th),'_',thetanam,'_',tags,'.png'])
end



%% Plot Distribution Of parameters compared to Bayesian Results assuming normality in Bayesian results

mdl=3;
% Read Posterior
d2 = 'F:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\';
[m1, c1, thetaNames] = ExtractGaussParamsTheta(1);

% Frequentist
m2 = bestx.pe_results.fit.thetabest';
c2 = bestx.pe_results.fit.g_var_cov_mat;
sd1 = zeros(1,length(m1));
sd2 = zeros(1,length(m2));
for i=1:length(m1)
    sd1(1,i) = sqrt(c1(i,i));
    sd2(1,i) = sqrt(c2(i,i));
end

 [parNams, parMax, parMin] = getThetaBounds(1);

for th=1:length(m2)
    h = figure;
    hold on 
    % Bayesian
    top1 = (m1(th)+5*sd1(th));
    bot1 = (m1(th)-5*sd1(th));
    x1 = bot1:((top1-bot1)/200):top1;
    y1 = normpdf(x1,m1(th),sd1(th));
    % Frequentist
    top2 = (m2(th)+5*sd2(th));
    bot2 = (m2(th)-5*sd2(th));
    x2 = bot2:((top2-bot2)/200):top2;
    y2 = normpdf(x2,m2(th),sd2(th));
    
    plot(x1,y1)
    plot(x2,y2)
    
    % Bounds
    plot([parMax(th), parMax(th)],[0,max([y1, y2])])
    plot([parMin(th), parMin(th)],[0,max([y1, y2])])
    
    ylabel('Density')
    thetanam = erase(thetaNames{th},"_raw");
    xlabel(strrep(thetanam, '_', '\_'))
    legend('Bayesian', 'Frequentist')
    set(gca, 'XScale', 'log')
    hold off
%     saveas(h, [Folder,'\Plots\KernelDensity_',num2str(th),'_',thetanam,'_',tags,'.png'])
end
























