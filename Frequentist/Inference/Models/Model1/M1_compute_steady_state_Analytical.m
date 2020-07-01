function [ res ] = M1_compute_steady_state_Analytical(theta, InitialStates_AU, initial_u)
% ToggleSwitch_M1_Compute_SteadyState computes the steady state of the MToggleSwitch model for the given values of
% theta and the inputs u_IPTG and u_aTc.

% Model 3
k_in_iptg = theta(1);
k_out_iptg = theta(2);
k_in_aTc = theta(3);
k_out_aTc = theta(4);

kL_p_m0 = theta(5);
kL_p_m = theta(6);
theta_T = theta(7);
theta_aTc = theta(8);
n_aTc = theta(9);
n_T = theta(10);
kT_p_m0 = theta(11);
kT_p_m = theta(12);
theta_L = theta(13);
theta_IPTG = theta(14);
n_IPTG = theta(15);
n_L = theta(16);


preRFP = InitialStates_AU(1);
preGFP = InitialStates_AU(2);

u_IPTG = initial_u(1);
u_aTc = initial_u(2);

%% Steady state equation

L_RFP = ((1/0.1386)*(kL_p_m0 + (kL_p_m/(1+((preGFP/theta_T)*(1/(1+(u_aTc/theta_aTc)^n_aTc)))^n_T))))/0.0165;
          
T_GFP = ((1/0.1386)*(kT_p_m0 + (kT_p_m/(1+((preRFP/theta_L)*(1/(1+(u_IPTG/theta_IPTG)^n_IPTG)))^n_L))))/0.0165;
          

res = [u_IPTG u_aTc L_RFP T_GFP];

end
