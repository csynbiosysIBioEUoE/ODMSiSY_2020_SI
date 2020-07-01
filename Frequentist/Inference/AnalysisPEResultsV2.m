%% Script to extract the best PE run
% Extract the iteration yielding the best cost function Value without
% considering the validation step

d1 = 'F:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\ParameterEstimation'; % Main Directory
d2 = 'F:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\'; % MAin directory to save
dirm = [d1, '\PE_Final\PE_M3_CV_Rep2\']; % Directory of the PE results
tags = 'PE_Model3_CrosValidation_Rep2-'; % Common string of the result files
mdl = 3; % Model being analysed 

%% Cost Function Analysis

% Plot cost function
figure
hold on
for i=1:100
    try
        x = load([dirm,'\',tags,num2str(i),'PEMultiTest']);
        stairs(x.pe_results.nlpsol.conv_curve(:,1), x.pe_results.nlpsol.conv_curve(:,2))
    catch
    end
end
% set(gca, 'YScale', 'log')

% Extract run with the best cost function value
cfv = nan(1,100);
for i=1:100
    try
        x = load([dirm,'\',tags,num2str(i),'PEMultiTest']);
        cfv(i) = x.pe_results.nlpsol.fbest;
    catch
    end
end

bcfv = min(cfv);
bestind = find(cfv==bcfv);

bestx = load([dirm,'\',tags,num2str(bestind),'PEMultiTest']);

 % Save the results
BESTPE.Model = mdl;
BESTPE.indx = bestind;
BESTPE.bcfv = bcfv;
BESTPE.theta =  bestx.pe_results.fit.thetabest;
BESTPE.cov = bestx.pe_results.fit.g_var_cov_mat;
BESTPE.fullresults = bestx;

save([dirm,'\BestResults_',tags,'.mat'],'BESTPE')


%% Plot simulations agains data 
if ~isfolder([dirm,'\Plots'])
    mkdir([dirm,'\Plots'])
end

for i=1:bestx.pe_inputs.exps.n_exp
    h = figure;
    subplot(2,1,1)
    hold on 
    errorbar(bestx.pe_inputs.exps.t_s{i}, bestx.pe_inputs.exps.exp_data{i}(:,1), bestx.pe_inputs.exps.error_data{i}(:,1), 'r')
    plot(bestx.pe_results.sim.tsim{i}, simMo.sim.states{i}(:,3), 'b')% For UpdatedY0
%     plot(bestx.pe_results.sim.tsim{i}, bestx.pe_results.sim.states{i}(:,3), 'b')
    ylabel('RFP')
    xlabel('time(min)')

    subplot(2,1,2)
    hold on 
    errorbar(bestx.pe_inputs.exps.t_s{i}, bestx.pe_inputs.exps.exp_data{i}(:,2), bestx.pe_inputs.exps.error_data{i}(:,2), 'g')
    plot(bestx.pe_results.sim.tsim{i}, simMo.sim.states{i}(:,4), 'b')% For UpdatedY0
%     plot(bestx.pe_results.sim.tsim{i}, bestx.pe_results.sim.states{i}(:,4), 'b')
    ylabel('GFP')
    xlabel('time(min)')
    
%     saveas(h, [dirm,'\Plots\BestFit_Exp',num2str(i),'_', tags,'.png'])
    
end


%% Plot Distribution Of parameters compared to Bayesian Results

% Read Posterior
posterior = readtable([d2, 'PosteriorDistributionsBayes\draws_ALL_Model',num2str(mdl),'.csv']);
% Extract data
thetaNames = posterior.Properties.VariableNames;
thetaMatrix = posterior.Variables;

% Frequentist
m2 = BESTPE.theta';
c2 = BESTPE.cov;
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
%     saveas(h, [dirm,'\Plots\KernelDensity_',num2str(th),'_',thetanam,'_',tags,'.png'])
end











