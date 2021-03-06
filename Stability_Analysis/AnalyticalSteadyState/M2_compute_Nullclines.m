function [ res ] = M2_compute_Nullclines(theta, InitialStates_AU, initial_u)
% ToggleSwitch_M1_Compute_SteadyState computes the analytical steady state of the MToggleSwitch model 2 for the given values of
% theta and the inputs u_IPTG and u_aTc.

% Model 3
k_in_iptg = theta(1);
k_in_aTc = theta(2);

kL_p_m0 = theta(3);
kL_p_m = theta(4);
theta_T = theta(5);
theta_aTc = theta(6);
n_aTc = theta(7);
n_T = theta(8);
kT_p_m0 = theta(9);
kT_p_m = theta(10);
theta_L = theta(11);
theta_IPTG = theta(12);
n_IPTG = theta(13);
n_L = theta(14);


preRFP = InitialStates_AU(:,1);
preGFP = InitialStates_AU(:,2);

u_IPTG = initial_u(1);
u_aTc = initial_u(2);

%% Steady state equation

L_RFP = ((1/0.1386).*(kL_p_m0 + (kL_p_m./(1+((preGFP./theta_T).*(1./(1+(u_aTc./theta_aTc).^n_aTc))).^n_T))))/0.0165;
          
T_GFP = ((1/0.1386).*(kT_p_m0 + (kT_p_m./(1+((preRFP./theta_L).*(1./(1+(u_IPTG./theta_IPTG).^n_IPTG))).^n_L))))/0.0165;
          

res = [L_RFP T_GFP];

end
