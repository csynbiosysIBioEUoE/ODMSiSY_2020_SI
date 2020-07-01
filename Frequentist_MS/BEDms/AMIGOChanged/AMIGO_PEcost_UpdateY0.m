
% Change the name of the AMIGO2 function AMIGO_PE to
% AMIGO_PE_UpdateY0. 

% Find the call to AMIGO_ivpsol and change it to
% AMIGO_ivpsol_NoVer


% Add the call to AMIGO_transform_Y0_TS after line 100 (after calling
% AMIGO_transform_theta). 

% This modification is needed if you want to update Y0 at each function
% evaluation during PE since these will deppend on the parameter vector
% selected. 


