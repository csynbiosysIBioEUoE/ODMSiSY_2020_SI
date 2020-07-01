%% Run scripts for PE
% This scripts allows you to run the OED for mdoel selection as did in the
% paper. 

% This script does NOT consider crossvalidation


% Inputs: 
%       - resultBase: String to be attached to the save results as an
%       identifier for the run
%       - numExperiments: Number of runs with different initial guesses for
%       theta
%       - mdl: Integer between 1 and 3 indicating the Toggle Switch model
%       for which the parameter bounds are desired. 


function [] = MasterRun_fit_to_AllModels(resultBase, numExperiments, mdl)

% Create a matrix of initial guesses for the parameters, having as many
% rows as the number of PE iterations (numExperiments)
% Each vector is passed as input to the computing function

% [~, theta_max, theta_min] = getThetaBounds(mdl)
%               
%               
% M_norm = lhsdesign(numExperiments,length(theta_min));
% M = zeros(size(M_norm));
% for c=1:size(M_norm,2)
%     for r=1:size(M_norm,1)
%         M(r,c) = 10^(M_norm(r,c)*(log10(theta_max(1,c))-log10(theta_min(1,c)))+log10(theta_min(1,c))); % log exploration
%     end
% end 
% 
% % %check the location of the parameters that are fixed
% ParFull = M;  
% save(['MatrixParameters_Model',num2str(mdl),'.mat'],'ParFull');

ParFull = load(['MatrixParameters_Model',num2str(mdl),'.mat']);

% This is to run the for loop for initial guesses that have not been
% considered yet (to solve issues with the parfor stopping because it lost
% connection to a worker)
Folder='.';
filePattern = fullfile(Folder, strcat(resultBase,'-','*PEMultiTest.mat'));
Files = dir(filePattern);
fi = strings(1,length(Files));
for i=1:length(Files)
    fi(i) = Files(i).name;
end


parfor epcc_exps=1:numExperiments
    m = [resultBase,'-',num2str(epcc_exps),'.mat'];
    
    if ~ismember(m,fi)
        
        try
            global_theta_guess = ParFull.ParFull(epcc_exps,:);
            epccOutputResultFileNameBase = [resultBase,'-',num2str(epcc_exps)];
            [out] = Masterfit_to_AllModels(epccOutputResultFileNameBase,epcc_exps,global_theta_guess, mdl);
        catch err
            %open file
            errorFile = [resultBase,'-',num2str(epcc_exps),'.errorLog'];
            fid = fopen(errorFile,'a+');
            fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
            % close file
            fclose(fid);
        end
    end
end

