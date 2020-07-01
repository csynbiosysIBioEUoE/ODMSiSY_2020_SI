%% Cost function for Average Case
% This script computs the cost function as the product of euclidean distance
% between the RFP and GFP simulations of the 2 models. 

% Inputs: 
%       - od: Input profile to be tested
%       - inputs,results,privstruct: AMIGO structures

function [f,g,h] = OEDModelSelectionCostMultiplic(od,inputs,results,privstruct)


% Set the experimental structure
   
    inputsSIM.model = inputs.model;
    inputsSIM.pathd.results_folder = inputs.pathd.results_folder;                        
    inputsSIM.pathd.short_name     = inputs.pathd.short_name;
    inputsSIM.pathd.runident       = inputs.pathd.runident;
    clear newExps;
    newExps.n_exp = 1;                                         % Number of experiments 
    newExps.n_obs{1}=2;                                        % Number of observables per experiment        
    newExps.obs_names{1} = char('LacI_M1','TetR_M1', 'LacI_M3','TetR_M3');
    newExps.obs{1} = char('LacI_M1 = L_RFP','TetR_M1 = T_GFP', 'LacI_M3 = L_RFP2','TetR_M3 = T_GFP2');% Name of the observables 
    newExps.exp_y0{1}=inputs.exps.exp_y0{1};                                      % Initial condition for the experiment    

    newExps.t_f{1}=inputs.exps.t_f{1};                                   % Experiment duration
    newExps.n_s{1}=inputs.exps.t_f{1}/5;                             % Number of sampling times
    newExps.t_s{1}=0:5:inputs.exps.t_f{1} ;                              % Times of samples

    newExps.u_interp{1}='step';                                % Interpolating function for the input
    % newExps.u_interp{2}='step';                                % Interpolating function for the input
    newExps.n_steps{1}=inputs.exps.n_steps{1};                  % Number of steps in the input
    newExps.u{1}= [od(1:inputs.exps.n_steps{1}); od(inputs.exps.n_steps{1}+1:end)]+1e-7;                                     % IPTG and aTc values for the input
    newExps.t_con{1}=inputs.exps.t_con{1};                     % Switching times
    %     newExps.t_con{1}=[0 3000];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mock the experiment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inputsSIM.exps = newExps;

    % SIMULATION
    inputsSIM.ivpsol.ivpsolver=inputs.ivpsol.ivpsolver;
    inputsSIM.ivpsol.senssolver=inputs.ivpsol.senssolver;
    inputsSIM.ivpsol.rtol=inputs.ivpsol.rtol;
    inputsSIM.ivpsol.atol=inputs.ivpsol.atol;

    inputsSIM.plotd.plotlevel='noplot';

    % Simulate the models 
    y = AMIGO_SModel_NoVer(inputsSIM);
    
    % Extract each simualtion vector
	IPTGi =y.sim.states{1}(:,1);
	aTci  =y.sim.states{1}(:,2);
	L_RFP =y.sim.states{1}(:,3);
	T_GFP =y.sim.states{1}(:,4);
	IPTGi2=y.sim.states{1}(:,5);
	aTci2 =y.sim.states{1}(:,6);
	L_RFP2=y.sim.states{1}(:,7);
	T_GFP2=y.sim.states{1}(:,8);
% 	u_IPTG=[	 od(1)	 od(2)	];
% 	u_aTc =[	 od(3)	 od(4)	];

% Definition of the cost function to be minimised (average euclidean
% distance)

    % Euclidean distance for RFP simulations
    subs = L_RFP-L_RFP2;
    sqr = zeros(length(subs),1);
    for i=1:length(sqr)
        sqr(i,1) = subs(i,1)^2;
    end

	f1 = sqrt(sum(sqr));
    
    % Euclidean distance for GFP simulations
    subs2 = T_GFP-T_GFP2;
    sqr2 = zeros(length(subs2),1);
    for i=1:length(sqr2)
        sqr2(i,1) = subs2(i,1)^2;
    end

	f2 = sqrt(sum(sqr2));
    
    % Cost function
    f = -(f1*f2);
	 h(1)=0;
	 g(1)=0;

return


















end


