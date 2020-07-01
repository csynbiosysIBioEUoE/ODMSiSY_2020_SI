%% OED for Model Selection main script
% This script specifies all the experimental conditions and inputs to be
% introduced to AMIGO2 in order to perform the OED for model selection.

% Inputs: 
%       - epccOutputResultFileNameBase: string as a tag for the
%       optimisation results file. 
%       - steps: Integer indicating the number of steps of the experiment
%       to be optimised
%       - initGuess: Initial guess for the input profile (random initial
%       point for the optimisation process)
%       - casee: string indicating the type of combination between the
%       system outputs to be followed.
%           * RFP: Only difference in RFP levels of the 2 models is
%           considered
%           * GFP: Only difference in GFP levels of the 2 models is
%           considered
%           * Aver: Average between GFP and RFP levels of the 2 models is
%           considered
%           * Mult: Product between GFP and RFP levels of the 2 models is
%           considered
%           * MultiObj: Multi-Objective optimisation is considered


function [out] = OEDModelSelection_Master(epccOutputResultFileNameBase,steps,initGuess, casee)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Specify folder name and short_name
        results_folder = strcat('OEDModelSelection_',casee,'_',datestr(now,'yyyy-mm-dd-HHMMSS'));
        short_name     = strcat('OEDMS',int2str(steps));

        % Load model. 
        
        model = ToggleSwitch_load_model_M1vsM3_ModelSelection;
        
        % Initial guesses for theta 
        
        global_theta_guess = model.par';

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Run OED
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Compile the model
        clear inputs;
        inputs.model = model;
        

        inputs.pathd.results_folder = results_folder;                        
        inputs.pathd.short_name     = short_name;
        inputs.pathd.runident       = 'initial_setup';
        

        % Compute the steady state considering the initial theta guess, u_IPTG and
        % u_aTc
        y0 = zeros(1,model.n_st);
        y0(1,:) = M1vsM3_Compute_SteadyState_OverNight_ModelSelection(inputs,model.par,...
                        [28.510, 1363.193, 28.510, 1363.193],...
                        [1, 0]+1e-7);  
        
        inputs.model.par = global_theta_guess;

        inputs.DOsol.y0=y0;                               %Initial conditions
        inputs.DOsol.tf_type='fixed';                          %Process duration type: fixed or free
        inputs.DOsol.tf_guess=24*60;                               %Process duration
        
        % Specify number of cost functions to be evaluated
        switch casee
            case {'RFP','GFP','Aver','Mult'}
                inputs.DOsol.N_DOcost = 1;
            case 'MultiObj'
                inputs.DOsol.N_DOcost = 2;
        end

        %COST FUNCTION
        inputs.DOsol.DOcost_type='min';                        %Type of problem: max/min
        inputs.DOsol.user_cost = 1;
        
        %CVP (Control Vector Parameterization) DETAILS
        inputs.DOsol.u_interp='stepf';                         %Control definition
                                                               %'sustained' |'stepf'|'step'|'linear'|
        inputs.DOsol.n_steps=steps;
        inputs.DOsol.u_guess=initGuess;                             % Initial guess for the input
        inputs.DOsol.u_min= [0*ones(1,inputs.DOsol.n_steps); 0*ones(1,inputs.DOsol.n_steps)]+1e-7;
        inputs.DOsol.u_max=[1*ones(1,inputs.DOsol.n_steps); 100*ones(1,inputs.DOsol.n_steps)];% Minimum and maximum value for the input
        inputs.DOsol.t_con=0:inputs.DOsol.tf_guess/inputs.DOsol.n_steps:24*60;         % Input swithching times, including intial and
                                                  
        %==================================
        % NUMERICAL METHDOS RELATED DATA
        %==================================

        % SIMULATION
        inputs.ivpsol.ivpsolver='cvodes';
        inputs.ivpsol.senssolver='cvodes';%'fdsens5'

        inputs.ivpsol.rtol=1.0D-13;
        inputs.ivpsol.atol=1.0D-13;

        %OPTIMIZATION

        switch casee
            case {'RFP','GFP','Aver','Mult'}
                inputs.nlpsol.reopt='off'; 
                inputs.nlpsol.nlpsolver='eSS';                    % [] NLP solver:

                inputs.nlpsol.eSS.maxeval = 60000;               % Maximum number of cost function evaluations
                inputs.nlpsol.eSS.maxtime = 50000;                  % Maximum computational time in seconds
                inputs.nlpsol.eSS.local.solver = 'fmincon';       % Local solver- SQP
                inputs.nlpsol.eSS.local.finish = 'fmincon';    % Local solver- Direct method

            case 'MultiObj'
                inputs.nlpsol.nlpsolver='nsga2';                      % Solves the problem using a multi-objective
                                                      % optimizer- nsga2
                inputs.nlpsol.nsga2.popsize = 10*inputs.DOsol.n_steps;% nsga2 population size
                inputs.nlpsol.nsga2.maxGen  = 10*inputs.DOsol.n_steps;% nsga2 max number of generations
                inputs.nlpsol.nsga2.plotInterval =50;                 % nsga2 interval between two calls of "plotnsga".
                %inputs.nlpsol.nsga2.initfun={@initpop, 'populations.txt'} %nsga2 use previous results to refine
                inputs.nlpsol.nsga2.mutation={'gaussian',0.1, 0.5};   % nsga2 type of mutation
        end
        
        
        %================================
        % CALL AMIGO2 from COMMAND LINE
        %================================
        % It is recommended to keep all inputs in a 'problem_file'.m.

        AMIGO_Prep(inputs);
        
        % Specify the user defined cost function for each one of the cases.
        % AMIGO_Prep does not recognise this entries and I do not know why
        % since they can be used by AMIGO_DO. At some point I will modify
        % the AMIGO_Prep scripts so it does allow them. 
        
        switch casee
            case 'RFP'
                inputs.pathd.DO_function = 'OEDModelSelectionCostRFP';
                inputs.pathd.DO_constraints = 'OEDModelSelectionConstRFP';
            case 'GFP'
                inputs.pathd.DO_function = 'OEDModelSelectionCostGFP';
                inputs.pathd.DO_constraints = 'OEDModelSelectionConstGFP';
            case 'Aver'
                inputs.pathd.DO_function = 'OEDModelSelectionCostAverage';
                inputs.pathd.DO_constraints = 'OEDModelSelectionConstAverage';
            case 'Mult'
                inputs.pathd.DO_function = 'OEDModelSelectionCostMultiplic';
                inputs.pathd.DO_constraints = 'OEDModelSelectionConstMultiplic';
            case 'MultiObj'
                inputs.pathd.DO_function = 'OEDModelSelectionCostMultiObject';
                inputs.pathd.DO_constraints = 'OEDModelSelectionConstMultiObject';
        end
        
        % Run optimisation
        oed_results = AMIGO_DO(inputs);
        
        % Save results
        save(strcat(epccOutputResultFileNameBase,'OED_ModelSelection.mat'),'oed_results','inputs','steps', 'casee');

out = 1;
end
