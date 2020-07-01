%% PE main script
% This script specifies all the experimental conditions and inputs to be
% introduced to AMIGO2 in order to perform Parameter Estimation with the
% Lugagne et al., 2017. data considered in the Bayesian inference

% This script does NOT consider crossvalidation

% Inputs: 
%       - epccOutputResultFileNameBase: string as a tag for the
%       optimisation results file. 
%       - epcc_exps: Integer indicating the number of the run
%       - global_theta_guess: Initial guess for theta
%       - mdl: Integer between 1 and 3 indicating the Toggle Switch model
%       for which the parameter bounds are desired. 

function [out] = Masterfit_to_AllModels(epccOutputResultFileNameBase,epcc_exps,global_theta_guess,mdl)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Specify folder name and short_name
        results_folder = strcat('MFit',datestr(now,'yyyy-mm-dd-HHMMSS'));
        short_name     = strcat('MF',int2str(epcc_exps));
        
        rng('default')
        resultFileName = [strcat(epccOutputResultFileNameBase),'.dat'];
        rng shuffle;
        rngToGetSeed = rng;

        % Write the header information of the .dat file 
        fid = fopen(resultFileName,'w');
        fprintf(fid,'HEADER DATE %s\n', datestr(datetime()));
        fprintf(fid,'HEADER RANDSEED %d\n',rngToGetSeed.Seed);
        fclose(fid);
        
        
        % Load experimental data. 
        LugDat = load('D:\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\LugagneData\BayesSelectDataLugagne.mat');

        % Read the model into the model variable. 
        switch mdl
            case 1
                model = ToggleSwitch_load_model_M1;
            case 2
                model = ToggleSwitch_load_model_M2;
            case 3
                model = ToggleSwitch_load_model_M3;
            otherwise
                disp('Model NOT selected correctly!!!')
                return
        end
        
        % Get theta bounds
        [~, global_theta_max, global_theta_min] = getThetaBounds(mdl);
        % Initial guesses for theta 
        
        global_theta_guess = global_theta_guess';

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Run  PE on the set
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Compile the model
        clear inputs;
        inputs.model = model;
        

        inputs.pathd.results_folder = results_folder;                        
        inputs.pathd.short_name     = short_name;
        inputs.pathd.runident       = 'initial_setup';
        

        % Compute the steady state considering the initial theta guess, u_IPTG and
        % u_aTc
        AMIGO_Prep(inputs);
        y0 = zeros(length(LugDat.Data.expName),model.n_st);
        switch mdl
            case 1
                for iexp=1:length(LugDat.Data.expName)
                    y0(iexp,:) = M1_Compute_SteadyState_OverNight(inputs,global_theta_guess,...
                        [LugDat.Data.exp_data{iexp}(1,1), LugDat.Data.exp_data{iexp}(1,2)],...
                        [LugDat.Data.Initial_IPTG{iexp} LugDat.Data.Initial_aTc{iexp}]+1e-7);
                end
            case 2
                for iexp=1:length(LugDat.Data.expName)
                    y0(iexp,:) = M2_Compute_SteadyState_OverNight(inputs,global_theta_guess,...
                        [LugDat.Data.exp_data{iexp}(1,1), LugDat.Data.exp_data{iexp}(1,2)],...
                        [LugDat.Data.Initial_IPTG{iexp} LugDat.Data.Initial_aTc{iexp}]+1e-7);
                end
            case 3
                for iexp=1:length(LugDat.Data.expName)
                    y0(iexp,:) = M3_Compute_SteadyState_OverNight(inputs,model.par,...
                        [LugDat.Data.exp_data{iexp}(1,1), LugDat.Data.exp_data{iexp}(1,2)],...
                        [LugDat.Data.Initial_IPTG{iexp} LugDat.Data.Initial_aTc{iexp}]+1e-7);
                end
        end
        
        inputs.model.par = global_theta_guess;
        clear exps;

        % Define the exps structure, containing the experimental data to fit
        exps.n_exp = length(LugDat.Data.expName);

        for iexp=1:length(LugDat.Data.expName)
            
            duration = round(LugDat.Data.t_samples{iexp}(end,1));                                                
            inducer = LugDat.Data.input{iexp}+1e-7;
            inputs.meta.init_U{iexp} = [LugDat.Data.Initial_IPTG{iexp} LugDat.Data.Initial_aTc{iexp}]+1e-7;
            exps.exp_type{iexp} = 'fixed'; 

            exps.n_obs{iexp} = 2; 
            exps.obs_names{iexp} = char('RFP_LacI','GFP_TetR');
            exps.obs{iexp} = char('RFP_LacI = L_RFP','GFP_TetR = T_GFP');

            exps.t_f{iexp} = round(LugDat.Data.t_samples{iexp}(end,1));
            exps.n_s{iexp} = LugDat.Data.n_samples{iexp}(1);
            sampl = round(LugDat.Data.t_samples{iexp}(10,1)-LugDat.Data.t_samples{iexp}(9,1));
            exps.t_s{iexp} = 0:sampl:duration; 

            exps.u_interp{iexp} = 'step';
            exps.t_con{iexp} = round(LugDat.Data.t_con{iexp}(:,1))'; 
            exps.n_steps{iexp} = length(LugDat.Data.input{iexp}(:,1));
            exps.u{iexp} = inducer';
            exps.data_type = 'real';
            exps.noise_type = 'homo';

            exps.exp_data{iexp} = LugDat.Data.exp_data{iexp};
            exps.error_data{iexp} = LugDat.Data.standard_dev{iexp};
            exps.exp_y0{iexp} = y0(iexp,:);
            
            
        end


        best_global_theta = global_theta_guess'; 

        % Compile the model
        inputs.model = model;
        inputs.model.par = best_global_theta;
        inputs.exps  = exps;

        inputs.pathd.results_folder = results_folder;                        
        inputs.pathd.short_name     = short_name;
        inputs.pathd.runident       = strcat('pe-',int2str(epcc_exps));


        % GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERMENTS)
        inputs.PEsol.id_global_theta=model.par_names;
        inputs.PEsol.global_theta_guess=(best_global_theta);
        inputs.PEsol.global_theta_max=global_theta_max;  % Maximum allowed values for the paramters
        inputs.PEsol.global_theta_min=global_theta_min;  % Minimum allowed values for the parameters


        %COST FUNCTION RELATED DATA
        inputs.PEsol.PEcost_type='llk';                       % 'lsq' (weighted least squares default) | 'llk' (log likelihood) | 'user_PEcost'
        inputs.PEsol.llk_type='homo_var';                      % [] To be defined for llk function, 'homo' | 'homo_var' | 'hetero'

        %SIMULATION
        inputs.ivpsol.ivpsolver='cvodes';
        inputs.ivpsol.senssolver='cvodes';
        inputs.ivpsol.rtol=1.0D-13;
        inputs.ivpsol.atol=1.0D-13;

        %OPTIMIZATION
        inputs.nlpsol.nlpsolver='eSS';
        inputs.nlpsol.eSS.maxeval = 200000;
        inputs.nlpsol.eSS.maxtime = 5000;
        inputs.nlpsol.eSS.local.solver = 'fmincon'; 
        inputs.nlpsol.eSS.local.finish = 'fmincon'; 

        inputs.plotd.plotlevel='noplot';

        
        AMIGO_Prep(inputs);
        
        pe_start = now;
        pe_inputs = inputs;
        results = AMIGO_PE(inputs);
%         results = AMIGO_PE_UpdateY0(inputs);
        pe_results = results;
        pe_end = now;

        %Save the best theta
        best_global_theta = results.fit.thetabest;

        %Write results to the output file
        fid = fopen(epccOutputResultFileNameBase,'a');
        used_par_names = model.par_names;

        for j=1:size(used_par_names,1)
            fprintf(fid,'PARAM_FIT %s %f\n', used_par_names(j,:), results.fit.thetabest(j));
        end

        %Time in seconds
        fprintf(fid,'PE_TIME %.1f\n', (pe_end-pe_start)*24*60*60);
        fclose(fid);

        save(strcat(epccOutputResultFileNameBase,'PEMultiTest.mat'),'pe_results','exps','pe_inputs','best_global_theta');

out = 1;
end
