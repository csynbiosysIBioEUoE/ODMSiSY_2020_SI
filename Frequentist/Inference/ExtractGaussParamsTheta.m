%% Gaussian approximation of Posterior Distributions for theta
% This script computes the mean and covariance matrix for each of the
% posterior distributions obtained in the RStan inference

% Inputs:
%       - modl: Integer between 1 and 3 indicating the Toggle Switch model
%       for which the parameter bounds are desired. 

function [Mean, Covs, thetaNames] = ExtractGaussParamsTheta(modl)

% Read csv file with all the RStan samples obtained
posterior = readtable(['draws_ALL_Model',num2str(modl),'.csv']);

% Extract data
thetaNames = posterior.Properties.VariableNames;
thetaMatrix = posterior.Variables;

% Compute mean and covariance matrix
Mean = mean(thetaMatrix);
Covs = cov(thetaMatrix);

% Display mean for quick check
disp(thetaNames)
disp(Mean)

end





























