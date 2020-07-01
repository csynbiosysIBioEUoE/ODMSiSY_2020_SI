% ------------------- PRIORS DEFINITION ------------------------

% Definition of priors based on parameter estimates in Lugagne et al.
% Parameters obtained upon fitting on the calibration data in Nat.Com 2017
%% 
clc, clear all, close all;

%% % Mean estimates for the parameters (same order of magnitude of estimates in Lugagne )
k_Lm0 = 3e-2;
k_Lm = 10;
k_Tm0 = 1e-1;
k_Tm = 1;
k_Lp = 1;
k_Tp = 1;

g_Lm = 1.386e-1;
g_Tm = g_Lm;

g_Lp = 1.65e-2;
g_Tp = g_Lp;

Theta_L = 30;
Theta_T = 30;
Theta_IPTG = 0.5;
Theta_aTc = 50;
n_L = 2;
n_T = 2;
n_aTc = 2;
n_IPTG = 2;

k_IPTG = 4e-2;
k_in_IPTG = 4e-2;
k_out_IPTG = 4e-2;

k_aTc = 1e-1;
k_in_aTc = 1e-1;
k_out_aTc = 1e-1;

%% For the unknown parameters, build normal prior centered on the mean and covering 0.1-10 times their value
Params = [k_Lm0*k_Lp k_Lm*k_Lp Theta_T Theta_aTc n_aTc n_T k_Tm0*k_Tp k_Tm*k_Tp Theta_L Theta_IPTG n_IPTG n_L k_IPTG k_in_IPTG k_out_IPTG k_aTc k_in_aTc k_out_aTc];
Mean = zeros(1,length(Params)+4);
Std = Mean;
Type = {};
Parameter = {'k_L_pm0','k_L_pm','Theta_T','Theta_aTc','n_aTc','n_T','k_T_pm0','k_T_pm','Theta_L','Theta_IPTG','n_IPTG','n_L','k_IPTG','k_in_IPTG','k_out_IPTG','k_aTc','k_in_aTc','k_out_aTc','g_Lm','g_Tm','g_Lp','g_Tp'};
hill = [5,6,11,12];
for i=1:length(Params)
    if ismember(i,hill)
        Parameter{i}
        start = 0;
        stop = 5;        
    else
        start = log(0.1*Params(i));
        stop = log(10*Params(i));
    end

    Mean(i) = mean([start,stop]);
    Std(i) = (stop-start)/4;
    Type = [Type; 'unknown'];
    
end
%% Modify priors for theta_aTc and theta_IPTG
% Theta_aTc
Mean(4) = mean([log(1),log(100)]);
Std(4) = (log(100)-log(1))/4;

% Theta_IPTG
Mean(10) = mean([log(0.01),log(1)]); 
Std(10) = (log(1)-log(0.01))/4;

%% Adding info on the fixed parameters
Mean(i+1:end) = [g_Lm, g_Tm,g_Lp,g_Tp];
Type{19,1} = 'fixed';
Type{20,1} = 'fixed';
Type{21,1} = 'fixed';
Type{22,1} = 'fixed';

%% Creating the table
index = linspace(1,length(Parameter),length(Parameter));
rowsi = strread(num2str(index),'%s');
varNames = {'Parameter_Name','Type','Mean','Std'};

Data= table(Parameter',Type,Mean',Std','RowNames',rowsi,'VariableNames',varNames);
writetable(Data,'InitialPriors_AllModels_Relaxed_Log.csv');    



