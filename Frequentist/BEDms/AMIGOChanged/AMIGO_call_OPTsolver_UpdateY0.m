
% Change the name of the AMIGO2 function AMIGO_call_OPTsolver to
% AMIGO_call_OPTsolver_UpdateY0. 

% Find all the strings "AMIGO_PEcost" and change them to
% "AMIGO_PEcost_UpdateY0".

% This modification is needed if you want to update Y0 at each function
% evaluation during PE since these will deppend on the parameter vector
% selected. 

