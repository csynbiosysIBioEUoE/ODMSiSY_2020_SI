%% Compute mono/bi-stable regions for Bayesian results
% The script computes the monostable and bistable regions for the 3 models
% given a set of inducers. To asses mono/bi-stability the code computes the
% nullclines for all the cases and assesses the number of intersections
% between these. 

%% Range of inducers
IPTG = 0:0.01:1;
aTc = 0:100;

%% Check if matrix already exist
% If not generate an empty one for pre-allocation
try
    load('BayesBiStabilityCheck\BayesBiStabilityMatrixModel1.mat', 'BStab1', 'Bvals1')
catch
    BStab1 = nan(101,101); % Matrix containing a tag for monostability and bi-stability
    Bvals1 = cell(101,101); % Structure containing the intersection points and RFP/GFP values
end 
try
    load('BayesBiStabilityCheck\BayesBiStabilityMatrixModel2.mat', 'BStab2', 'Bvals2')
catch
    BStab2 = nan(101,101); % Matrix containing a tag for monostability and bi-stability
    Bvals2 = cell(101,101); % Structure containing the intersection points and RFP/GFP values
end
try 
    load('BayesBiStabilityCheck\BayesBiStabilityMatrixModel3.mat', 'BStab3', 'Bvals3')
catch
    BStab3 = nan(101,101); % Matrix containing a tag for monostability and bi-stability
    Bvals3 = cell(101,101); % Structure containing the intersection points and RFP/GFP values
end

%% Model 1
% Read csv file with all the RStan samples obtained
posterior = readtable(['drawsRedT_ALL_Model',num2str(1),'.csv']);
% Extract data
thetaMatrix = posterior.Variables;

for i=1:length(IPTG)
    for j=1:length(aTc)
        Bvtemp = cell(400,1);
        temp1 = nan(1,400);
        if isnan(BStab1(i,j))
            parfor h=1:400
                if isempty(Bvtemp{h,1})
                    tht = thetaMatrix(h,:);
                    m1 = M1_compute_Nullclines(tht, [(0:3000)', (0:3000)'], [IPTG(i),aTc(j)]);
                    i1 = InterX([m1(:,1), (0:3000)']',[(0:3000)', m1(:,2)]');
                    Bvtemp{h,1} = i1;
                    if length(i1(1,:))==3
                        temp1(1,h)=1;
                    else
                        temp1(1,h)=0;
                    end
                end
            end
            Bvals1{i,j} = Bvtemp;
            BStab1{i,j} = mean(temp1, 'omitnan');
        end
    end
end
save('BayesBiStabilityCheck\BayesBiStabilityMatrixModel1.mat', 'BStab1', 'Bvals1')

h = figure;
pcolor(aTc, IPTG, BStab1)
ylabel('IPTG (nM)')
xlabel('aTc (ng/ml)')
title('Bayesian Stability Analysis Model 1') 
colorbar
saveas(h, ['BayesBiStabilityCheck\','BayesModel1Check.png'])


%% Model 2

% Read csv file with all the RStan samples obtained
posterior = readtable(['drawsRedT_ALL_Model',num2str(2),'.csv']);
% Extract data
thetaMatrix = posterior.Variables;

for i=1:length(IPTG)
    for j=1:length(aTc)
        Bvtemp = cell(400,1);
        temp1 = nan(1,400);
        if isnan(BStab2(i,j))
            parfor h=1:400
                if isempty(Bvtemp{h,1})
                    tht = thetaMatrix(h,:);
                    m1 = M2_compute_Nullclines(tht, [(0:3000)', (0:3000)'], [IPTG(i),aTc(j)]);
                    i1 = InterX([m1(:,1), (0:3000)']',[(0:3000)', m1(:,2)]');
                    Bvtemp{h,1} = i1;
                    if length(i1(1,:))==3
                        temp1(1,h)=1;
                    else
                        temp1(1,h)=0;
                    end
                end
            end
            Bvals2{i,j} = Bvtemp;
            BStab2(i,j) = mean(temp1, 'omitnan');
        end
    end
end
save('BayesBiStabilityCheck\BayesBiStabilityMatrixModel2.mat', 'BStab2', 'Bvals2')

h = figure;
pcolor(aTc, IPTG, BStab2)
ylabel('IPTG (nM)')
xlabel('aTc (ng/ml)')
title('Bayesian Stability Analysis Model 2') 
colorbar
saveas(h, ['BayesBiStabilityCheck\','BayesModel2Check.png'])



%% Model 3

% Read csv file with all the RStan samples obtained
posterior = readtable(['drawsRedT_ALL_Model',num2str(3),'.csv']);
% Extract data
thetaMatrix = posterior.Variables;

for i=1:length(IPTG)
    for j=1:length(aTc)
        Bvtemp = cell(400,1);
        temp1 = nan(1,400);
        if isnan(BStab3(i,j))
            parfor h=1:400
                if isempty(Bvtemp{h,1})
                    tht = thetaMatrix(h,:);
                    m1 = M3_compute_Nullclines(tht, [(0:3000)', (0:3000)'], [IPTG(i),aTc(j)]);
                    i1 = InterX([m1(:,1), (0:3000)']',[(0:3000)', m1(:,2)]');
                    Bvtemp{h,1} = i1;
                    if length(i1(1,:))==3
                        temp1(1,h)=1;
                    else
                        temp1(1,h)=0;
                    end
                end
            end
            Bvals3{i,j} = Bvtemp;
            BStab3{i,j} = mean(temp1, 'omitnan');
        end
    end
end
save('BayesBiStabilityCheck\BayesBiStabilityMatrixModel3.mat', 'BStab3', 'Bvals3')

h = figure;
pcolor(aTc, IPTG, BStab3)
ylabel('IPTG (nM)')
xlabel('aTc (ng/ml)')
title('Bayesian Stability Analysis Model 3') 
colorbar
saveas(h, ['BayesBiStabilityCheck\','BayesModel3Check.png'])















































