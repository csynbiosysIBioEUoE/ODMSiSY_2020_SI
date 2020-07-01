%% Analyse OED results

%% Define path of the results
foldw = 'E:\UNI\D_Drive\PhD\Year_1\2020_03_03_Processes2ToggleModelComparison\Scripts\1_FrequentistAnalysis\OED_ModelSelection\Results\'; % Main directory
fold = [foldw, '\Final\F_M1vsM3\']; % Folder containing the results
tag = 'OED1_M1vsM3_ThetaAMIGOlsq_Rep1-Steps_'; % Common file name for all the runs
steps = [2,4,6,8]; % Number of steps for each run (needed)
trial = 'M1vsM3_Bistability_'; % Tag that will be added in the name of the saved resulting plots

%% RFP case
if ~isfolder([fold,'Plots']) % Make folder to save plots
    mkdir([fold,'Plots'])
end

for i=1:4 % Loop for each experiment (in this case, defined by the number of steps)
    x = load([fold, tag, num2str(steps(i)), '-RFPOED_ModelSelection.mat']); % Load OED result
    % Simulations and input profile
    h = figure('Renderer', 'painters', 'Position', [100 100 900 600]); % Plot OED with fluorescent reporters and input profiles

    subplot(3,1,1)
    title(['RFP optimisation, ', num2str(steps(i)), ' Steps'])
    hold on
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,3), 'r')
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,7), 'b')
    ylabel('RFP')
    legend('M1', 'M3')
    subplot(3,1,2)
    hold on
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,4), 'g')
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,8), 'b')
    ylabel('GFP')
    legend('M1', 'M3')
    subplot(3,1,3)
    yyaxis left
    stairs(x.oed_results.do.t_con, [x.oed_results.do.u(1,:), x.oed_results.do.u(1,end)])
%     ylim([0,1])
    ylabel('IPTG')
    yyaxis right
    stairs(x.oed_results.do.t_con, [x.oed_results.do.u(2,:), x.oed_results.do.u(2,end)])
%     ylim([0,100])
    ylabel('aTc')
    xlabel('time (min)')
%     saveas(h, [fold,'Plots\',trial,'RFPOptim_',num2str(steps(i)),'Steps.png'])
end

% Cost functions plots
h2 = figure('Renderer', 'painters', 'Position', [50 50 900 600]);
hold on
for i=1:4
    x = load([fold, tag, num2str(steps(i)), '-RFPOED_ModelSelection.mat']);
    stairs([1;x.oed_results.nlpsol.conv_curve2(:,1)], [x.oed_results.nlpsol.fiter1; x.oed_results.nlpsol.conv_curve(:,2)])
    ylabel('CFV')
    xlabel('Iteration')
    title('Convergence curve RFP optimisation')
end
legend('2s', '4s','6s','8s')
saveas(h2, [fold,'Plots\',trial,'RFPOptim_CostFunction.png'])

%% GFP case

for i=1:4
    x = load([fold, tag, num2str(steps(i)), '-GFPOED_ModelSelection.mat']);
    % Simulations and input profile
    h = figure('Renderer', 'painters', 'Position', [100 100 900 600]);

    subplot(3,1,1)
    title(['GFP optimisation, ', num2str(steps(i)), ' Steps'])
    hold on
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,3), 'r')
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,7), 'b')
    ylabel('RFP')
    legend('M1', 'M3')
    subplot(3,1,2)
    hold on
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,4), 'g')
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,8), 'b')
    ylabel('GFP')
    legend('M1', 'M3')
    subplot(3,1,3)
    yyaxis left
    stairs(x.oed_results.do.t_con, [x.oed_results.do.u(1,:), x.oed_results.do.u(1,end)])
    ylabel('IPTG')
    yyaxis right
    stairs(x.oed_results.do.t_con, [x.oed_results.do.u(2,:), x.oed_results.do.u(2,end)])
    ylabel('aTc')
    xlabel('time (min)')
    saveas(h, [fold,'Plots\',trial,'GFPOptim_',num2str(steps(i)),'Steps.png'])
end

% Cost functions
h2 = figure('Renderer', 'painters', 'Position', [50 50 900 600]);
hold on
for i=1:4
    x = load([fold, tag, num2str(steps(i)), '-GFPOED_ModelSelection.mat']);
    stairs([1;x.oed_results.nlpsol.conv_curve2(:,1)], [x.oed_results.nlpsol.fiter1; x.oed_results.nlpsol.conv_curve(:,2)])
    ylabel('CFV')
    xlabel('Iteration')
    title('Convergence curve GFP optimisation')
end
legend('2s', '4s','6s','8s')
saveas(h2, [fold,'Plots\',trial,'GFPOptim_CostFunction.png'])



%% Average case


for i=1:4
    x = load([fold, tag, num2str(steps(i)), '-AverOED_ModelSelection.mat']);
    % Simulations and input profile
    h = figure('Renderer', 'painters', 'Position', [100 100 900 600]);

    subplot(3,1,1)
    title(['Average optimisation, ', num2str(steps(i)), ' Steps'])
    hold on
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,3), 'r')
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,7), 'b')
    ylabel('RFP')
    legend('M1', 'M3')
    subplot(3,1,2)
    hold on
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,4), 'g')
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,8), 'b')
    ylabel('GFP')
    legend('M1', 'M3')
    subplot(3,1,3)
    yyaxis left
    stairs(x.oed_results.do.t_con, [x.oed_results.do.u(1,:), x.oed_results.do.u(1,end)])
    ylabel('IPTG')
    yyaxis right
    stairs(x.oed_results.do.t_con, [x.oed_results.do.u(2,:), x.oed_results.do.u(2,end)])
    ylabel('aTc')
    xlabel('time (min)')
    saveas(h, [fold,'Plots\',trial,'AverageOptim_',num2str(steps(i)),'Steps.png'])
end

% Cost functions
h2 = figure('Renderer', 'painters', 'Position', [50 50 900 600]);
hold on
for i=1:4
    x = load([fold, tag, num2str(steps(i)), '-AverOED_ModelSelection.mat']);
    stairs([1;x.oed_results.nlpsol.conv_curve2(:,1)], [x.oed_results.nlpsol.fiter1; x.oed_results.nlpsol.conv_curve(:,2)])
    ylabel('CFV')
    xlabel('Iteration')
    title('Convergence curve Average optimisation')
end
legend('2s', '4s','6s','8s')
saveas(h2, [fold,'Plots\',trial,'AverageOptim_CostFunction.png'])




%% Multiplicative case

for i=1:4
    x = load([fold, tag, num2str(steps(i)), '-MultOED_ModelSelection.mat']);
    % Simulations and input profile
    h = figure('Renderer', 'painters', 'Position', [100 100 900 600]);

    subplot(3,1,1)
    title(['Multiplicative optimisation, ', num2str(steps(i)), ' Steps'])
    hold on
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,3), 'r')
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,7), 'b')
    ylabel('RFP')
    legend('M1', 'M3')
    subplot(3,1,2)
    hold on
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,4), 'g')
    plot(x.oed_results.sim.tsim{1}, x.oed_results.sim.states{1}(:,8), 'b')
    ylabel('GFP')
    legend('M1', 'M3')
    subplot(3,1,3)
    yyaxis left
    stairs(x.oed_results.do.t_con, [x.oed_results.do.u(1,:), x.oed_results.do.u(1,end)])
    ylabel('IPTG')
    yyaxis right
    stairs(x.oed_results.do.t_con, [x.oed_results.do.u(2,:), x.oed_results.do.u(2,end)])
    ylabel('aTc')
    xlabel('time (min)')
    saveas(h, [fold,'Plots\',trial,'MultipOptim_',num2str(steps(i)),'Steps.png'])
end

% Cost functions
h2 = figure('Renderer', 'painters', 'Position', [50 50 900 600]);
hold on
for i=1:4
    x = load([fold, tag, num2str(steps(i)), '-MultOED_ModelSelection.mat']);
    stairs([1;x.oed_results.nlpsol.conv_curve2(:,1)], [x.oed_results.nlpsol.fiter1; x.oed_results.nlpsol.conv_curve(:,2)])
    ylabel('CFV')
    xlabel('Iteration')
    title('Convergence curve Multiplicative optimisation')
end
legend('2s', '4s','6s','8s')
saveas(h2, [fold,'Plots\',trial,'MultipOptim_CostFunction.png'])




%% Multi-Objective case

% Plot Pareto-Fronts

h1 = figure('Renderer', 'painters', 'Position', [50 50 900 600]);
for i=1:4
    x = load([fold, tag, num2str(steps(i)), '-MultiObjOED_ModelSelection.mat']);
    if i==1
    subplot(2,2,1)
    title('2 Steps MOOptimisation Pareto front')
    xlabel('RFP')
    ylabel('GFP')
    hold on
    scatter(x.oed_results.nlpsol.pareto_obj(:,1), x.oed_results.nlpsol.pareto_obj(:,2))
    scatter(x.oed_results.nlpsol.pareto_obj(10,1), x.oed_results.nlpsol.pareto_obj(10,2), 'r', 'MarkerFaceColor', 'r') % Point of the pareto front that will be used to simulate the system and visualise the resultant trajectories (randomly selected, can be changed, but then needs to be changed in the plot section)
    end
    
    if i==2
    subplot(2,2,2)
    title('4 Steps MOOptimisation Pareto front')
    xlabel('RFP')
    ylabel('GFP')
    hold on
    scatter(x.oed_results.nlpsol.pareto_obj(:,1), x.oed_results.nlpsol.pareto_obj(:,2))
    scatter(x.oed_results.nlpsol.pareto_obj(15,1), x.oed_results.nlpsol.pareto_obj(15,2), 'r', 'MarkerFaceColor', 'r')% Point of the pareto front that will be used to simulate the system and visualise the resultant trajectories (randomly selected, can be changed, but then needs to be changed in the plot section)
    end
    
    if i==3
    subplot(2,2,3)
    title('6 Steps MOOptimisation Pareto front')
    xlabel('RFP')
    ylabel('GFP')
    hold on
    scatter(x.oed_results.nlpsol.pareto_obj(:,1), x.oed_results.nlpsol.pareto_obj(:,2))
    scatter(x.oed_results.nlpsol.pareto_obj(33,1), x.oed_results.nlpsol.pareto_obj(33,2), 'r', 'MarkerFaceColor', 'r')% Point of the pareto front that will be used to simulate the system and visualise the resultant trajectories (randomly selected, can be changed, but then needs to be changed in the plot section)
    end
    
    if i==4
    subplot(2,2,4)
    title('8 Steps MOOptimisation Pareto front')
    xlabel('RFP')
    ylabel('GFP')
    hold on 
    scatter(x.oed_results.nlpsol.pareto_obj(:,1), x.oed_results.nlpsol.pareto_obj(:,2))
    scatter(x.oed_results.nlpsol.pareto_obj(33,1), x.oed_results.nlpsol.pareto_obj(33,2), 'r', 'MarkerFaceColor', 'r')% Point of the pareto front that will be used to simulate the system and visualise the resultant trajectories (randomly selected, can be changed, but then needs to be changed in the plot section)
    end
    
end
saveas(h1, [fold,'Plots\',trial,'MultiOpOptim_ParetoFronts.png'])


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
exps.n_exp=0;
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
newExps.n_exp = 1;                                         % Number of experiments 
newExps.n_obs{1}=2;                                        % Number of observables per experiment        
newExps.obs_names{1} = char('LacI_M1','TetR_M1', 'LacI_M3','TetR_M3');
newExps.obs{1} = char('LacI_M1 = L_RFP','TetR_M1 = T_GFP', 'LacI_M3 = L_RFP2','TetR_M3 = T_GFP2');% Name of the observables 
newExps.exp_y0{1}=y0;                                      % Initial condition for the experiment    
newExps.t_f{1}=duration;                                   % Experiment duration
newExps.n_s{1}=duration/5 + 1;                             % Number of sampling times
newExps.t_s{1}=0:5:duration ;                              % Times of samples
newExps.u_interp{1}='step';                                % Interpolating function for the input
inputs.model = model;
inputs.exps = newExps;
inputs.exps.data_type='pseudo';
inputs.exps.noise_type='hetero_proportional';
inputs.exps.std_dev{1}=[0.0 0.0];
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='cvodes';
inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;
inputs.plotd.plotlevel='noplot';

AMIGO_Prep(inputs);

% simDa = AMIGO_SData(inputs);



for i=1:4
    x = load([fold, tag, num2str(steps(i)), '-MultiObjOED_ModelSelection.mat']);
    if i==1
        inputs.exps.n_steps{1} = steps(i);
        inputs.exps.u{1} = x.oed_results.do.u{10};
        inputs.exps.t_con{1} = x.oed_results.do.t_con{10};
        
        
        simMo = AMIGO_SModel(inputs);
        
        h1 = figure('Renderer', 'painters', 'Position', [200 200 900 600]);    
        subplot(3,1,1)
        title(['Multi-Objective optimisation, ', num2str(steps(i)), ' Steps'])
        hold on
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,3), 'r')
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,7), 'b')
        ylabel('RFP')
        legend('M1', 'M3')

        subplot(3,1,2)
        hold on
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,4), 'g')
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,8), 'b')
        ylabel('GFP')
        legend('M1', 'M3')

        subplot(3,1,3)
        yyaxis left
        stairs(x.oed_results.do.t_con{10}, [x.oed_results.do.u{10}(1,:), x.oed_results.do.u{10}(1,end)])
        ylabel('IPTG')
        yyaxis right
        stairs(x.oed_results.do.t_con{10}, [x.oed_results.do.u{10}(2,:), x.oed_results.do.u{10}(2,end)])
        ylabel('aTc')
        xlabel('time (min)')
        saveas(h1, [fold,'Plots\',trial,'MultiOpOptim_Simulations_Step_',num2str(steps(i)),'.png'])

    end
    
    if i==2
        inputs.exps.n_steps{1} = steps(i);
        inputs.exps.u{1} = x.oed_results.do.u{15};
        inputs.exps.t_con{1} = x.oed_results.do.t_con{15};
        
        
        simMo = AMIGO_SModel(inputs);
        
        h1 = figure('Renderer', 'painters', 'Position', [200 200 900 600]);    
        subplot(3,1,1)
        title(['Multi-Objective optimisation, ', num2str(steps(i)), ' Steps'])
        hold on
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,3), 'r')
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,7), 'b')
        ylabel('RFP')
        legend('M1', 'M3')

        subplot(3,1,2)
        hold on
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,4), 'g')
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,8), 'b')
        ylabel('GFP')
        legend('M1', 'M3')

        subplot(3,1,3)
        yyaxis left
        stairs(x.oed_results.do.t_con{15}, [x.oed_results.do.u{15}(1,:), x.oed_results.do.u{15}(1,end)])
        ylabel('IPTG')
        yyaxis right
        stairs(x.oed_results.do.t_con{15}, [x.oed_results.do.u{15}(2,:), x.oed_results.do.u{15}(2,end)])
        ylabel('aTc')
        xlabel('time (min)')
        saveas(h1, [fold,'Plots\',trial,'MultiOpOptim_Simulations_Step_',num2str(steps(i)),'.png'])
    end
    
    if i==3
        inputs.exps.n_steps{1} = steps(i);
        inputs.exps.u{1} = x.oed_results.do.u{33};
        inputs.exps.t_con{1} = x.oed_results.do.t_con{33};
        
        
        simMo = AMIGO_SModel(inputs);
        
        h1 = figure('Renderer', 'painters', 'Position', [200 200 900 600]);    
        subplot(3,1,1)
        title(['Multi-Objective optimisation, ', num2str(steps(i)), ' Steps'])
        hold on
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,3), 'r')
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,7), 'b')
        ylabel('RFP')
        legend('M1', 'M3')

        subplot(3,1,2)
        hold on
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,4), 'g')
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,8), 'b')
        ylabel('GFP')
        legend('M1', 'M3')

        subplot(3,1,3)
        yyaxis left
        stairs(x.oed_results.do.t_con{33}, [x.oed_results.do.u{33}(1,:), x.oed_results.do.u{33}(1,end)])
        ylabel('IPTG')
        yyaxis right
        stairs(x.oed_results.do.t_con{33}, [x.oed_results.do.u{33}(2,:), x.oed_results.do.u{33}(2,end)])
        ylabel('aTc')
        xlabel('time (min)')
        saveas(h1, [fold,'Plots\',trial,'MultiOpOptim_Simulations_Step_',num2str(steps(i)),'.png'])
    end
    
    if i==4
        inputs.exps.n_steps{1} = steps(i);
        inputs.exps.u{1} = x.oed_results.do.u{33};
        inputs.exps.t_con{1} = x.oed_results.do.t_con{33};
        
        
        simMo = AMIGO_SModel(inputs);
        
        h1 = figure('Renderer', 'painters', 'Position', [200 200 900 600]);    
        subplot(3,1,1)
        title(['Multi-Objective optimisation, ', num2str(steps(i)), ' Steps'])
        hold on
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,3), 'r')
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,7), 'b')
        ylabel('RFP')
        legend('M1', 'M3')

        subplot(3,1,2)
        hold on
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,4), 'g')
        plot(simMo.sim.tsim{1}, simMo.sim.states{1}(:,8), 'b')
        ylabel('GFP')
        legend('M1', 'M3')

        subplot(3,1,3)
        yyaxis left
        stairs(x.oed_results.do.t_con{33}, [x.oed_results.do.u{33}(1,:), x.oed_results.do.u{33}(1,end)])
        ylabel('IPTG')
        yyaxis right
        stairs(x.oed_results.do.t_con{33}, [x.oed_results.do.u{33}(2,:), x.oed_results.do.u{33}(2,end)])
        ylabel('aTc')
        xlabel('time (min)')
        saveas(h1, [fold,'Plots\',trial,'MultiOpOptim_Simulations_Step_',num2str(steps(i)),'.png'])
    end
    
end

















