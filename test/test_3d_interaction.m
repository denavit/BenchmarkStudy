clear all; close all; clc;

%% Input and Definitions
H  = 20;
B  = 12;
t  = 0.75;
Fy = 50;
fc = 4;
L  = 15*H;
units = 'US';
section_type = 'RCFT';

% Define Section
section = RCFT(H,B,t,Fy,fc,units);
section.Lx = L;
section.Ly = L;

% Define Frame Data
data = struct;
data.section        = section;
data.section_type   = section_type;
data.L              = L;
data.frame_typeZ    = 'Sidesway_Inhibited';
data.frame_typeY    = 'Sidesway_Inhibited';
data.Delta0Z        = 0.0;
data.Delta0Y        = 0.0;
data.delta0Z        = L/1000;
data.delta0Y        = 0.0;
data.MzTop          = 1.0;
data.MzBot          = 0;
data.MyTop          = 0;
data.MyBot          = -1.0;

% Set Design Options
alpha = 1.5;
StiffnessReduction = 0.8;
tau = 0.8; 

% Create base scratch path
base_scratch_path = fullfile(pathOf.scratch,'test_3d_interaction');
if ~exist(base_scratch_path,'dir')
    mkdir(base_scratch_path);
end

%% Curve 1
numPoints = 9;

% Define Analysis Options
analysisOptions = struct;
analysisOptions.controlled_dof_override = [];
analysisOptions.includeInitialGeometricImperfections = true;
analysisOptions.reportMaximumLimitpoint = true;  

% Define Fiber Section Definition
opts = struct;
opts.includePackageDefinition = true;
opts.nf1 = 30;
opts.nf2 = 30;
switch data.section_type
    case 'RCFT'
        opts.SteelMaterialType                  = 'ModifiedAbdelRahman';
        opts.AbdelRahmanResidualStressParameter = 0.75;
        opts.AbdelRahmanHardeningRatio          = 0.0;
        opts.ConcreteMaterialType               = 'ProposedForDesign';

    case 'SRC'
        opts.SteelMaterialType                  = 'ElasticPP';
        opts.ConcreteMaterialType               = 'ProposedForDesign';                    

    case 'WF'
        opts.SteelMaterialType                  = 'ElasticPP';

    otherwise 
        error('Unknown section type: %s',data.section_type)
end
FiberSectionDefinitonOptions = opts;
FiberSectionDefinition_SectionID        = 1;
FiberSectionDefinition_StartMatID       = 1;
section_def = FiberSectionDefinition(data.section,'3d',...
    FiberSectionDefinition_SectionID,...
    FiberSectionDefinition_StartMatID,...
    FiberSectionDefinitonOptions);

% Initilize Output
P1 = nan(numPoints,1);
X1 = nan(numPoints,1);
status1 = cell(numPoints,1);

% Run Axial Only Analysis
BA = BenchmarkAnalysis3d_OpenSees(data,section_def,analysisOptions);
BA.scratchPath = fullfile(base_scratch_path,'Axial');

results = BA.runAnalysis('LimitPoint_Proportional',0,[],1);
if ~results.limitPoint.good
    results = BA.runAnalysis('LimitPoint_Proportional',0,[],2);
end
if ~results.limitPoint.good
    results = BA.runAnalysis('LimitPoint_Proportional',0,[],3);
end                    
if ~results.limitPoint.good
    error('Limit Point Not Obtained: Axial Only');
end

P1(1) = results.limitPoint.P1;
X1(1) = results.limitPoint.X1;
status1{1} = results.limitPoint.limit_type;

% Run Non-proportional Analyses
Ps = linspace(P1(1),0,numPoints);
parfor i = 2:length(Ps)
    BA = BenchmarkAnalysis3d_OpenSees(data,section_def,analysisOptions);
    BA.runZeroAxialLoadAsSection = false;
    BA.scratchPath = fullfile(base_scratch_path,sprintf('Lateral_%i',i));
    
    results = BA.runAnalysis('LimitPoint_NonProportional',Ps(i),0,1);
    if ~results.limitPoint.good
        results = BA.runAnalysis('LimitPoint_NonProportional',Ps(i),0,2);
    end
    if ~results.limitPoint.good
        error('Limit Point Not Obtained: Lateral Loading');
    end
    
    P1(i) = results.limitPoint.P1;
    X1(i) = results.limitPoint.X1;
    status1{i} = results.limitPoint.limit_type;
end


%% Curve 4 - General Calculations

% Modify Data For Elastic Analysis
elastic_data = data;
elastic_data.sectionType = 'elastic';
if strcmp(data.frame_typeZ,'Sidesway_Uninhibited')
    elastic_data.kqtopZ = StiffnessReduction * data.kqtopZ; 
    elastic_data.kqbotZ = StiffnessReduction * data.kqbotZ;              
end
if strcmp(data.frame_typeY,'Sidesway_Uninhibited')
    elastic_data.kqtopY = StiffnessReduction * data.kqtopY; 
    elastic_data.kqbotY = StiffnessReduction * data.kqbotY;
end
elastic_data.Delta0Z        = 0.0;
elastic_data.Delta0Y        = 0.0;
elastic_data.delta0Z        = 0.0;
elastic_data.delta0Y        = 0.0;

% Define Analysis Options
analysisOptions = struct;
analysisOptions.controlled_dof_override = 1;
analysisOptions.includeInitialGeometricImperfections = false;
analysisOptions.reportMaximumLimitpoint = false;  
analysisOptions.geomTransfType = 'Corotational';

% Define Fiber Section Definition
[E,A,Iz,Iy,G,J] = elastic_data.section.sectionPropertiesForElasticAnalysis3d('ColumnStrength');
section_def = sprintf('section Elastic 1 %g %g %g %g %g %g \n',...
    StiffnessReduction*E,A,tau*Iz,tau*Iy,StiffnessReduction*G,J);

% Determine Interaction Strength
[Pa,~]   = section.pointA;
[Pc,Mcz] = section.pointC('strong');
[ ~,Mcy] = section.pointC('weak');
xi  = section.stabilityReduction('min',section.Pnco);
xPa = xi*abs(Pa);
xPc = xi*abs(Pc);

% Determine Critical Load
% @todo - can this be estimated from simple equations (e.g., pi^2EI/KL^2)
BA = BenchmarkAnalysis3d_OpenSees(elastic_data,section_def,analysisOptions);
BA.scratchPath = fullfile(base_scratch_path,'Elastic_Critical_Load');  
results = BA.runAnalysis('LimitPoint_Proportional',0,[],1);
Peb  = results.limitPoint.P1;

%% Curve 4 - cross section by cross section evaluation
numPoints = 9;

% Initilize Output
P4 = nan(numPoints,1);
X4 = nan(numPoints,1);

% Run Axial Only Analysis
P0 = -max(1.1*xPa,abs(Peb));
BA.scratchPath = fullfile(base_scratch_path,'Elastic_Axial');  
results = BA.runAnalysis('TargetForce_Proportional',P0,0,1);

Pr  = abs(results.path.P2);
Mrz = abs(results.path.M2z);
Mry = abs(results.path.M2y);
P4_loc = nan(1,size(Pr,2)); 
for i = 1:size(Pr,2)
    % Calculate strength ratio a each point along the length of the column 
    R = Interaction_Check_ACB(Pr(:,i),Mrz(:,i),Mry(:,i),xPa,xPc,Mcz,Mcy,alpha);
    [ind,x] = find_limit_point_in_vector(R,1.0);
    P4_loc(i) = ind+x;
end
ind = floor(min(P4_loc));
x = min(P4_loc) - ind;
P4(1) = interpolate_vector(results.path.P1,ind,x);
X4(1) = 0;

% Run Non-proportional Analyses
Ps = linspace(P4(1),0,numPoints)';
for i = 2:numPoints
    ALR = 1.1 * max([Mcz,Mcy]); 
    BA.scratchPath = fullfile(base_scratch_path,'Elastic_Lateral');
    results = BA.runAnalysis('TargetForce_NonProportional',Ps(i),ALR,1);
    
    Pr  = abs(results.path.P2);
    Mrz = abs(results.path.M2z);
    Mry = abs(results.path.M2y);
    X4_loc = nan(1,size(Pr,2));
    for j = 1:size(Pr,2)
        R = Interaction_Check_ACB(Pr(:,j),Mrz(:,j),Mry(:,j),xPa,xPc,Mcz,Mcy,alpha);
        [ind,x] = find_limit_point_in_vector(R,1.0);
        if isempty(ind)
            % @todo - should we raise an error here?
            ind = nan;
            x = nan;
        end
        X4_loc(j) = ind+x;
    end
    ind = floor(min(X4_loc));
    x = min(X4_loc) - ind;  
    X4(i) = interpolate_vector(results.path.time,ind,x);
    P4(i) = interpolate_vector(results.path.P1,ind,x);
end

%% Curve 4 - maximum along the length evaluation
numPoints = 9;

% Initilize Output
P4m = nan(numPoints,1);
X4m = nan(numPoints,1);

% Run Axial Only Analysis
% @todo - this is the same analysis as for the cross section by cross
% section analysis - it could be eliminted. 
P0 = -max(1.1*xPa,abs(Peb));
BA.scratchPath = fullfile(base_scratch_path,'Elastic_Axial');  
results = BA.runAnalysis('TargetForce_Proportional',P0,0,1);

Pr  = max(abs(results.path.P2),[],2);
Mrz = max(abs(results.path.M2z),[],2);
Mry = max(abs(results.path.M2y),[],2);
R = Interaction_Check_ACB(Pr,Mrz,Mry,xPa,xPc,Mcz,Mcy,alpha);
[ind,x] = find_limit_point_in_vector(R,1.0);
P4m(1) = interpolate_vector(results.path.P1,ind,x);
X4m(1) = 0;

% Run Non-proportional Analyses
P0 = linspace(P4m(1),0,numPoints)';
for i = 2:numPoints
    ALR = 1.1 * max([Mcz,Mcy]);
    BA.scratchPath = fullfile(base_scratch_path,'Elastic_Lateral');
    results = BA.runAnalysis('TargetForce_NonProportional',Ps(i),ALR,1);

    Pr  = max(abs(results.path.P2),[],2);
    Mrz = max(abs(results.path.M2z),[],2);
    Mry = max(abs(results.path.M2y),[],2);
    R   = Interaction_Check_ACB(Pr,Mrz,Mry,xPa,xPc,Mcz,Mcy,alpha);
    [ind,x] = find_limit_point_in_vector(R,1.0);
    X4m(i) = interpolate_vector(results.path.time,ind,x);
    P4m(i) = interpolate_vector(results.path.P1,ind,x);
end


%% Make Plots
fs = figureStyle('Display');

fs.figure(5,5)
fs.axes;

plot(X1, -P1, '-xc','LineWidth',1)
plot(X4, -P4, '-sg','LineWidth',1)
plot(X4m,-P4m,'-or','LineWidth',1)

xlabel('X')
ylabel('Axial Compression (kips)')

legend('Curve 1','Curve 4','Curve 4 max','Location','NE')

axis_limits('Margin',0.1)
axis_limits('Quadrant',1)


%% Remove Base Scratch Path
rmdir(base_scratch_path,'s')