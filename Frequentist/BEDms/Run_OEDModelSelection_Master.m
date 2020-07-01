%% Run scripts for OED for model selection
% This scripts allows you to run the OED for mdoel selection as did in the
% paper. 

% Inputs: 
%       - resultBase: String to be attached to the save results as an
%       identifier for the run


function [] = Run_OEDModelSelection_Master(resultBase)

% Add necessary folders into the path
addpath AMIGOChanged
addpath CostFunctions

% Vector with all the number of steps for the experiment to be optimised
% that we consider
stepmain = [2,4,6,8];
initvals = cell(1,4);

% Random initial guess for the input profile (it gets saved in the results
% so the same can be used)
for i=1:4
    initvals{1,i} = [randsample(0:0.01:1,stepmain(i), true); randsample(0:100,stepmain(i), true)]+1e-7;
end

% Identifiers for the type of combination of the 2 reporters
cases = {'RFP','GFP','Aver','Mult','MultiObj'};

% This is to run the for loop for initial guesses that have not been
% considered yet (to solve issues with the parfor stopping because it lost
% connection to a worker)
Folder='.';
filePattern = fullfile(Folder, strcat(resultBase,'-','*OED_ModelSelection.mat'));
Files = dir(filePattern);
fi = strings(1,length(Files));
for i=1:length(Files)
    fi(i) = Files(i).name;
end


parfor epcc_exps=1:numExperiments
    m = [resultBase,'-',num2str(epcc_exps),'.mat'];
    
    if ~ismember(m,fi)
        
        try
            steps = stepmain(epcc_exps);
            initGuess = initvals{1,epcc_exps};
            
            for j = 1:length(cases)
                try
                    casee = cases{j};
                    epccOutputResultFileNameBase = [resultBase,'-',num2str(steps),'-',casee];
                    [out] = OEDModelSelection_Master(epccOutputResultFileNameBase,steps,initGuess, casee);
    %                 [out] = OEDModelSelection_Master_ShortIter(epccOutputResultFileNameBase,steps,initGuess, casee);
                catch
                    errorFile = [resultBase,'-',casee,'_',num2str(steps),'.errorLog'];
                    fid = fopen(errorFile,'a+');
                    fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
                    % close file
                    fclose(fid);
                end
            end
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

