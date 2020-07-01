
% Change the name of the AMIGO2 function AMIGO_PE to
% AMIGO_PE_UpdateY0. 

% Find the call to AMIGO_check_model and change it to
% AMIGO_check_model_NoVer

% Find the call to AMIGO_check_exps and change it to AMIGO_check_exps_NoVer

% Find the call to AMIGO_call_OPTsolver and change it to
% AMIGO_call_OPTsolver_UpdateY0

% Find the calls to AMIGO_PEcost and change it to AMIGO_PEcost_UpdateY0

% Add the call to AMIGO_transform_Y0_TS after line 248 (after definition of
% privstruct.theta

% This modification is needed if you want to update Y0 at each function
% evaluation during PE since these will deppend on the parameter vector
% selected. 


