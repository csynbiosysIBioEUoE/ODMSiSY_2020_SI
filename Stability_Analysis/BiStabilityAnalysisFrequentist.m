%% Compute mono/bi-stable regions for Frequentist results
% The script computes the monostable and bistable regions for the 3 models
% given a set of inducers. To asses mono/bi-stability the code computes the
% nullclines for all the cases and assesses the number of intersections
% between these. 

%% Load models
M1 = ToggleSwitch_load_model_M1;
M2 = ToggleSwitch_load_model_M2;
M3 = ToggleSwitch_load_model_M3;

%% Range of inducers
IPTG = 0:0.01:1;
aTc = 0:100;

%% Check if matrix already exist
% If not generate an empty one for pre-allocation
try
    load('BiStabilityCheck\BiStabilityMatrixModel1.mat', 'Stab1', 'vals1')
catch
    Stab1 = nan(101,101); % Matrix containing a tag for monostability and bi-stability
    vals1 = cell(101,101); % Structure containing the intersection points and RFP/GFP values
end 
try
    load('BiStabilityCheck\BiStabilityMatrixModel2.mat', 'Stab2', 'vals2')
catch
    Stab2 = nan(101,101);% Matrix containing a tag for monostability and bi-stability
    vals2 = cell(101,101);% Structure containing the intersection points and RFP/GFP values
end
try 
    load('BiStabilityCheck\BiStabilityMatrixModel3.mat', 'Stab3', 'vals3')
catch
    Stab3 = nan(101,101);% Matrix containing a tag for monostability and bi-stability
    vals3 = cell(101,101);% Structure containing the intersection points and RFP/GFP values
end

%% Loop Model 1
for i=1:length(IPTG)
    for j=1:length(aTc)
        if isnan(Stab1(i,j))
            tht = M1.par;
            m1 = M1_compute_Nullclines(tht, [(0:3000)', (0:3000)'], [IPTG(i),aTc(j)]);
            i1 = InterX([m1(:,1), (0:3000)']',[(0:3000)', m1(:,2)]');
            vals1{i,j} = i1;
            if length(i1(1,:))==3
                Stab1(i,j)=1; % Bistable
            else
                Stab1(i,j)=0; % Monostable
            end
        end
    end
end
save('BiStabilityCheck\BiStabilityMatrixModel1.mat', 'Stab1', 'vals1') % Save resuts

% Plot results
h = figure;
pcolor(aTc, IPTG, Stab1)
ylabel('IPTG (nM)')
xlabel('aTc (ng/ml)')
title('Stability Analysis Model 1') 
saveas(h, ['BiStabilityCheck\','Model1Check.png'])


%% Loop Model 2
for i=1:length(IPTG)
    for j=1:length(aTc)
        if isnan(Stab2(i,j))
            tht = M2.par;
            m1 = M2_compute_Nullclines(tht, [(0:3000)', (0:3000)'], [IPTG(i),aTc(j)]);
            i1 = InterX([m1(:,1), (0:3000)']',[(0:3000)', m1(:,2)]');
            vals2{i,j} = i1;
            if length(i1(1,:))==3
                Stab2(i,j)=1;
            else
                Stab2(i,j)=0;
            end
        end
    end
end
save('BiStabilityCheck\BiStabilityMatrixModel2.mat', 'Stab2', 'vals2')

h = figure;
pcolor(aTc, IPTG, Stab2)
ylabel('IPTG (nM)')
xlabel('aTc (ng/ml)')
title('Stability Analysis Model 2') 
saveas(h, ['BiStabilityCheck\','Model2Check.png'])


%% Loop Model 3
for i=1:length(IPTG)
    for j=1:length(aTc)
        if isnan(Stab3(i,j))
            tht = M3.par;
            m1 = M3_compute_Nullclines(tht, [(0:3000)', (0:3000)'], [IPTG(i),aTc(j)]);
            i1 = InterX([m1(:,1), (0:3000)']',[(0:3000)', m1(:,2)]');
            vals3{i,j} = i1;
            if length(i1(1,:))==3
                Stab3(i,j)=1;
            else
                Stab3(i,j)=0;
            end
        end
    end
end
save('BiStabilityCheck\BiStabilityMatrixModel3.mat', 'Stab3', 'vals3')

h = figure;
pcolor(aTc, IPTG, Stab3)
ylabel('IPTG (nM)')
xlabel('aTc (ng/ml)')
title('Stability Analysis Model 3') 
saveas(h, ['BiStabilityCheck\','Model3Check.png'])






