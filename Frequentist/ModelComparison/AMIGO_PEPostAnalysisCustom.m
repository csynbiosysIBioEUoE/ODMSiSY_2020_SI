
% Change the name of the AMIGO2 function AMIGO_PEPostAnalysis to
% AMIGO_PEPostAnalysisCustom. 

% Modify lines 71 to 79 to: 
%     for iexp = 1:nexp
%         postAnalysis.residuals{iexp} = PEresults.sim.states{iexp}(:,3:4) - PEresults.sim.exp_data.Data.exp_data{iexp};
%         postAnalysis.data{iexp} = PEresults.sim.exp_data.Data.exp_data{iexp};
%         postAnalysis.std_dev{iexp} = PEresults.sim.exp_data.Data.standard_dev{iexp};
%         postAnalysis.n_globalParameters = PEresults.model.n_par;
%         postAnalysis.n_localParameters{iexp} = PEresults.model.n_par;
%         postAnalysis.obs_names{iexp} = PEinputs.exps.obs_names{iexp};
%        % postAnalysis.costType = PEinputs.PEsol.PEcost_type;
%     end

% Add between lines 262 and 263: 
%     [aic, ~] = AMIGO_aic(Rexp{iexp}(:),Dexp{iexp}(:),sexp{iexp}(:),npars,'homo_var');
%     [bic] = AMIGO_bic(Rexp{iexp}(:),Dexp{iexp}(:),sexp{iexp}(:),npars,'homo_var');
%     stats(iexp).AkaikeInformationCrit.AIC = aic;
%     stats(iexp).BayesianInformationCrit = bic;










