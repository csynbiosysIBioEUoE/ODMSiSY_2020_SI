%% Cost function for Average Case
% This script computs the cost function as the product of euclidean distance
% between the RFP and GFP simulations of the 2 models. 

% Used in script ComparisonBayesianVsFrequentistAverage.m


% Inputs: 
%       - od: Input profile to be tested
%       - inputs,results,privstruct: AMIGO structures

function [f,g,h] = OEDMSCostGFP_ComparisonBF(sims)

    y = sims;
    
    % Extract each simualtion vector
	IPTGi =y(:,1);
	aTci  =y(:,2);
	L_RFP =y(:,3);
	T_GFP =y(:,4);
	IPTGi2=y(:,5);
	aTci2 =y(:,6);
	L_RFP2=y(:,7);
	T_GFP2=y(:,8);
% 	u_IPTG=[	 od(1)	 od(2)	];
% 	u_aTc =[	 od(3)	 od(4)	];

% Definition of the cost function to be minimised (average euclidean
% distance)

    % Euclidean distance for RFP simulations
%     subs = L_RFP-L_RFP2;
%     sqr = zeros(length(subs),1);
%     for i=1:length(sqr)
%         sqr(i,1) = subs(i,1)^2;
%     end
% 
% 	f1 = sqrt(sum(sqr));
    
    % Euclidean distance for GFP simulations
    subs2 = T_GFP-T_GFP2;
    sqr2 = zeros(length(subs2),1);
    for i=1:length(sqr2)
        sqr2(i,1) = subs2(i,1)^2;
    end

	f2 = sqrt(sum(sqr2));
    
    % Cost function
    f = -(f2);
	 h(1)=0;
	 g(1)=0;

return


















end


