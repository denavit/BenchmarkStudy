clear all; close all; clc;

%% Plot Options
fs = figureStyle('Main');
calcUnits = unitSystem('US');
plotUnits = unitSystem('SI');
plotUnits.forceUnits = 'MN';
plotUnits.momentUnits = 'kN-m';
max_M = 850;
max_P = 10;
color_CS    = [0 0 0];
color_MMa   = [ 68 119 170]/255;
color_MMb   = [102 204 238]/255;
color_MMc   = [ 43 136  51]/255;
color_SOE   = [204 187  68]/255; 
color_SOI   = [238 102 119]/255;
color_extra = [170  51 119]/255;
color_grey  = [187 187 187]/255;

%% Analysis Options 
% it is assumed that include_strength_reduction = false;
% it is assumed that include_stiffness_reduction = false;
use_stored_inelastic_results = true;
compare_to_main_study = false;

%% Input
axis = 'x';

% Material Properties
fc  = 6;
Ec  = 57*sqrt(1000*fc);
fy  = 60;
fu  = 90;
fyt = 60;
Es  = 29000;

% Cross-Sectional Geometric Properties
H = 20;
B = 20;
conc_cross_section = Rectangle_Shape(H,B);

% Define transverse reinforcement
transverse_bar_size = reinf_bar_lookup('#3');
s = 16;
transverse_reinf_type = 'ties';
nLegX = 3;
nLegY = 3;

% Define longitudinal reinforcement
longitudinal_bar_size = reinf_bar_lookup('#8');
cover = 1.5;
nBarX = 3;
nBarY = 3;
Bc = B - 2*cover - 2*transverse_bar_size.diameter - longitudinal_bar_size.diameter;
Hc = H - 2*cover - 2*transverse_bar_size.diameter - longitudinal_bar_size.diameter;
longitudinal_config = sprintf('%ix-%iy',nBarX,nBarY);
reinforcement = reinf_rect(Bc,Hc,0,0,nBarX,nBarY,longitudinal_bar_size.area,longitudinal_bar_size.diameter);

% Define cross section
section = RC(fc,fy,conc_cross_section,reinforcement,'US');
section.transverse_reinf_type = transverse_reinf_type;
section.Ec = Ec;

% Computed Cross Section Properties
Po  = section.Po();
Ast = section.Asr();
Ag  = section.Ag();
h   = section.depth(axis);
Ig  = section.Ig(axis);

% Frame Properties
beta = 1.0;
L_over_H = 20;
L = L_over_H*H;

% Define data structure 
data = struct;
data.units = 'US';                            
data.section_type = 'RC';
data.shape = 'square';
data.B = B;
data.H = H;
data.cover = cover;
data.fc = fc;
data.Ec = Ec;
data.L = L;
data.nBarX = nBarX;
data.nBarY = nBarY;
data.Es = Es;
data.fy = fy;
data.fu = fu;
data.db = longitudinal_bar_size.diameter;
data.Ab = longitudinal_bar_size.area;
data.fyt = fyt;
data.dbt = transverse_bar_size.diameter;
data.Abt = transverse_bar_size.area;
data.nLegX = nLegX;
data.nLegY = nLegY;
data.s = s;
data.longitudinal_config = longitudinal_config;
%data.transverse_config = transverse_config;
data.longitudinal_bar_size = longitudinal_bar_size;
data.transverse_bar_size = transverse_bar_size;
data.section = section;
data.axis = axis;
data.frame_type = 'Sidesway_Inhibited';
data.beta = beta;
data.delta0 = L/1000;

%% Cross Section
hf = fs.picture(2.5,2.5);
section.plotSection();
axis_limits('Margin',0.1);

%% Benchmark Analysis Object
BA_ACI = BenchmarkAnalysis2d_ACI(data);
BA_ACI.include_strength_reduction = false;
BA_ACI.include_stiffness_reduction = false;
P_cs = BA_ACI.id2d_nom.idY;
M_cs = BA_ACI.id2d_nom.idX;
Mn = BA_ACI.id2d_nom.findXgivenY(0,'Pos');

%% Moment Magnification - EI options (a) and (b)
num_points = 50;
BA_ACI.EIeff_type = 'a';
BA_ACI.second_order_moment_ratio_limit = Inf;
[P_MMa_nl,M1_MMa_nl,~,M2_MMa_nl] = BA_ACI.design_interaction(num_points);
BA_ACI.second_order_moment_ratio_limit = 1.4;
[P_MMa,M1_MMa,~,M2_MMa] = BA_ACI.design_interaction(num_points);
BA_ACI.EIeff_type = 'b';
BA_ACI.second_order_moment_ratio_limit = Inf;
[P_MMb_nl,M1_MMb_nl,~,M2_MMb_nl] = BA_ACI.design_interaction(num_points);
BA_ACI.second_order_moment_ratio_limit = 1.4;
[P_MMb,M1_MMb,~,M2_MMb] = BA_ACI.design_interaction(num_points);

hf = fs.figure(3.25,3.25);
ha = fs.axes([0.10 0.10 0.85 0.87]);
plot(unitConvert('moment',M_cs,calcUnits,plotUnits),...
    -unitConvert('force',P_cs,calcUnits,plotUnits),...
    'Color',color_CS)

plot(unitConvert('moment',M1_MMa_nl,calcUnits,plotUnits),...
    -unitConvert('force',P_MMa_nl,calcUnits,plotUnits),...
    '--','Color',color_MMa)
% plot(unitConvert('moment',M2_MMa_nl,calcUnits,plotUnits),...
%     -unitConvert('force',P_MMa_nl,calcUnits,plotUnits),...
%     '-.','Color',color_MMa)
plot(unitConvert('moment',M1_MMa,calcUnits,plotUnits),...
    -unitConvert('force',P_MMa,calcUnits,plotUnits),...
    'Color',color_MMa)
% plot(unitConvert('moment',M2_MMa,calcUnits,plotUnits),...
%     -unitConvert('force',P_MMa,calcUnits,plotUnits),...
%     '-.','Color',color_MMa)

plot(unitConvert('moment',M1_MMb_nl,calcUnits,plotUnits),...
    -unitConvert('force',P_MMb_nl,calcUnits,plotUnits),...
    '--','Color',color_MMb)
% plot(unitConvert('moment',M2_MMb_nl,calcUnits,plotUnits),...
%     -unitConvert('force',P_MMb_nl,calcUnits,plotUnits),...
%     '-.','Color',color_MMb)
plot(unitConvert('moment',M1_MMb,calcUnits,plotUnits),...
    -unitConvert('force',P_MMb,calcUnits,plotUnits),...
    'Color',color_MMb)
% plot(unitConvert('moment',M2_MMb,calcUnits,plotUnits),...
%     -unitConvert('force',P_MMb,calcUnits,plotUnits),...
%     '-.','Color',color_MMb)

xlim([0 max_M])
ylim([0 max_P])
xlabel(sprintf('Bending Moment (%s)',plotUnits.momentUnits))
ylabel(sprintf('Axial Compression (%s)',plotUnits.forceUnits))

save_figure(hf,'example_figure_1','svg')

%% Moment Magnification - EI option (c)
num_points = 100;
BA_ACI.EIeff_type = 'c';
BA_ACI.second_order_moment_ratio_limit = Inf;
[P_MMc_nl,M1_MMc_nl,~,M2_MMc_nl] = BA_ACI.design_interaction(num_points);
BA_ACI.second_order_moment_ratio_limit = 1.4;
[P_MMc,M1_MMc,~,M2_MMc] = BA_ACI.design_interaction(num_points);

hf = fs.figure(3.25,3.25);
ha = fs.axes([0.10 0.10 0.85 0.87]);
plot(unitConvert('moment',M1_MMc_nl,calcUnits,plotUnits),...
    -unitConvert('force',P_MMc_nl,calcUnits,plotUnits),...
    '--','Color',color_MMc)
plot(unitConvert('moment',M2_MMc_nl,calcUnits,plotUnits),...
    -unitConvert('force',P_MMc_nl,calcUnits,plotUnits),...
    '-.','Color',color_grey)
plot(unitConvert('moment',M1_MMc,calcUnits,plotUnits),...
    -unitConvert('force',P_MMc,calcUnits,plotUnits),...
    'Color',color_MMc)
plot(unitConvert('moment',M2_MMc,calcUnits,plotUnits),...
    -unitConvert('force',P_MMc,calcUnits,plotUnits),...
    '-.','Color',color_grey)
plot(unitConvert('moment',M_cs,calcUnits,plotUnits),...
    -unitConvert('force',P_cs,calcUnits,plotUnits),...
    'Color',color_CS)
xlim([0 max_M])
ylim([0 max_P])
xlabel(sprintf('Bending Moment (%s)',plotUnits.momentUnits))
ylabel(sprintf('Axial Compression (%s)',plotUnits.forceUnits))

save_figure(hf,'example_figure_2','svg')

% Behavior a specific load
P  = unitConvert('force',2.4,'MN',calcUnits);
nM = 100;
M2 = linspace(0,unitConvert('moment',max_M,plotUnits,calcUnits),nM);
M1 = nan(1,nM);
I_over_Ig = nan(1,nM);
Cm = 1; 

for i = 1:nM   
    I = (0.80+25*Ast/Ag)*(1-M2(i)/(P*h)-0.5*P/Po)*Ig;
    if I < 0.35*Ig
        I = 0.35*Ig;
    end
    if I > 0.875*Ig
        I = 0.875*Ig;
    end
    Pc = pi^2*Ec*I/L^2;
    I_over_Ig(i) = I/Ig;
    if P < Pc
        delta = max(Cm./(1-P/Pc),1);
        M1(i) = M2(i)/delta;
    end
end

hf = fs.figure(3.25,2.25);
ha = fs.axes([0.13 0.17 0.82 0.80]);
max_M1 = 500;

plot([0 max_M],[0 max_M/1.4],'--r')

M2cs = BA_ACI.id2d_nom.findXgivenY(-P,'Pos');
plot(unitConvert('moment',M2cs,calcUnits,plotUnits)*[1 1],[0 max_M1],'--r')

plot(unitConvert('moment',M2,calcUnits,plotUnits),...
    unitConvert('moment',M1,calcUnits,plotUnits),...
    'Color','k')
xlim([0 max_M])
ylim([0 max_M1])
xlabel(sprintf('Maximum Internal Bending Moment (%s)',plotUnits.momentUnits))
ylabel(sprintf('Applied Bending Moment (%s)',plotUnits.momentUnits))

save_figure(hf,'example_figure_3','svg')

%% Second-Order Elastic Analysis
num_points = 100;

elastic_frame = BenchmarkAnalysis2d_Elastic(data,'RC_Study_Without_Stiffness_Reduction');
elastic_frame.includeInitialGeometricImperfections = false;
notionalLoadObject = notional_load(0.000,0.000,Inf);
designStrengthType = 'ACI';

[P_SOE_nl,M1_SOE_nl,~,M2_SOE_nl] = elastic_frame.designInteraction(...
    section,axis,designStrengthType,num_points,...
    notionalLoadObject,'zero',1.0,1.0,'none');

P_at_PeakMomentRatio = elastic_frame.determineLoadForMomentRatio(1.4,1.0,1.0,'none',nan);
P_SOE = max(P_SOE_nl,-P_at_PeakMomentRatio);
M1_SOE = M1_SOE_nl;
M2_SOE = M2_SOE_nl;

hf = fs.figure(3.25,3.25);
ha = fs.axes([0.10 0.10 0.85 0.87]);
plot(unitConvert('moment',M1_SOE_nl,calcUnits,plotUnits),...
    -unitConvert('force',P_SOE_nl,calcUnits,plotUnits),...
    '--','Color',color_SOE)
% plot(unitConvert('moment',M2_SOE_nl,calcUnits,plotUnits),...
%     -unitConvert('force',P_SOE_nl,calcUnits,plotUnits),...
%     '-.','Color',color_grey)
plot(unitConvert('moment',M1_SOE,calcUnits,plotUnits),...
    -unitConvert('force',P_SOE,calcUnits,plotUnits),...
    'Color',color_SOE)
% plot(unitConvert('moment',M2_SOE,calcUnits,plotUnits),...
%     -unitConvert('force',P_SOE,calcUnits,plotUnits),...
%     '-.','Color',color_grey)
plot(unitConvert('moment',M_cs,calcUnits,plotUnits),...
    -unitConvert('force',P_cs,calcUnits,plotUnits),...
    'Color',color_CS)
xlim([0 max_M])
ylim([0 max_P])
xlabel(sprintf('Bending Moment (%s)',plotUnits.momentUnits))
ylabel(sprintf('Axial Compression (%s)',plotUnits.forceUnits))

save_figure(hf,'example_figure_4','svg')

%% Second-Order Inelastic Analysis
if use_stored_inelastic_results
    % read data from file
    [~,S] = csvread2('saved_example_results.csv');
    P_SOI  = S.P_SOI;
    M1_SOI = S.M1_SOI;
    M2_SOI = S.M2_SOI;
else
    num_points = 40;

    % Initilize results
    M1_SOI = nan(num_points,1);
    P1_SOI = nan(num_points,1);
    M2_SOI = nan(num_points,1);
    P2_SOI = nan(num_points,1);

    % Define Analysis Object
    fiber_section_definition_options = struct;
    fiber_section_definition_options.nf1 = 30;
    fiber_section_definition_options.includePackageDefinition = true;
    fiber_section_definition_options.conc_material = 'Concrete04';
    fiber_section_definition_options.steel_material = 'ElasticPP';    
    
    section_definition = FiberSectionDefinition(data,axis,1,1,fiber_section_definition_options);
    
    analysis_options = struct;    
    analysis_options.absoluteStrainLimit = 0.05;
    analysis_options.store_extra_data = false;
    
    BA_OPS = BenchmarkAnalysis2d_OpenSees(data,section_definition,analysis_options);
    BA_OPS.scratchPath = pathOf.scratch;

    % Run axial only analysis to get Pn
    iResults = BA_OPS.runAnalysis('LimitPoint_Proportional',0,[],1);

    if ~iResults.limitPoint.good
        iResults = BA_OPS.runAnalysis('LimitPoint_Proportional',0,[],2);
    end

    if ~iResults.limitPoint.good
        error('Limit Point Not Obtained: Case %i, Axial Only',iData);
    end

    M1_SOI(1) = iResults.limitPoint.M1;
    P1_SOI(1) = iResults.limitPoint.P1;
    M2_SOI(1) = iResults.limitPoint.M2;
    P2_SOI(1) = iResults.limitPoint.P2;

    % Run nonproportional analyses at different axial loads
    Ps = linspace(P1_SOI(1),0,num_points);
    for i = 2:length(Ps)
        if Ps(i) ~= 0
            iResults = BA_OPS.runAnalysis('LimitPoint_NonProportional',Ps(i),[],1);

            if ~iResults.limitPoint.good
                iResults = BA_OPS.runAnalysis('LimitPoint_NonProportional',Ps(i),[],2);
            end

            if ~iResults.limitPoint.good
                iResults = BA_OPS.runAnalysis('LimitPoint_NonProportional',Ps(i),[],3);
            end

            if ~iResults.limitPoint.good
                error('Limit Point Not Obtained: Axial Load Level %i (%g)',i,Ps(i));
            end

            M1_SOI(i) = iResults.limitPoint.M1;
            P1_SOI(i) = iResults.limitPoint.P1;
            M2_SOI(i) = iResults.limitPoint.M2;
            P2_SOI(i) = iResults.limitPoint.P2;
        else
            iResults = BA_OPS.runSectionAnalysis('LimitPoint_NonProportional',0,[],1);

            if ~iResults.limitPoint.good
                iResults = BA_OPS.runSectionAnalysis('LimitPoint_NonProportional',0,[],2);
            end

            if ~iResults.limitPoint.good
                error('Limit Point Not Obtained: Axial Load Level %i (0)',i);
            end                

            M1_SOI(i) = iResults.limitPoint.M1;
            P1_SOI(i) = 0;
            M2_SOI(i) = iResults.limitPoint.M2;
            P2_SOI(i) = 0;
        end
    end  
    
    P_SOI = P1_SOI;
    
    %fprintf('P,M1,M2\n')
    %for i = 1:num_points
    %    fprintf('%.5f,%.5f,%.5f\n',P_SOI(i),M1_SOI(i),M2_SOI(i))
    %end
end 
  
hf = fs.figure(3.25,3.25);
ha = fs.axes([0.10 0.10 0.85 0.87]);
plot(unitConvert('moment',M1_SOI,calcUnits,plotUnits),...
    -unitConvert('force',P_SOI,calcUnits,plotUnits),...
    'Color',color_SOI)
plot(unitConvert('moment',M2_SOI,calcUnits,plotUnits),...
    -unitConvert('force',P_SOI,calcUnits,plotUnits),...
    '-.','Color',color_grey)
plot(unitConvert('moment',M_cs,calcUnits,plotUnits),...
    -unitConvert('force',P_cs,calcUnits,plotUnits),...
    'Color',color_CS)
xlim([0 max_M])
ylim([0 max_P])
xlabel(sprintf('Bending Moment (%s)',plotUnits.momentUnits))
ylabel(sprintf('Axial Compression (%s)',plotUnits.forceUnits))

save_figure(hf,'example_figure_5','svg')

%% Appliled Load Comparison
hf = fs.figure(3.25,3.25);
ha = fs.axes([0.15 0.12 0.80 0.85]);

plot(M1_MMa/Mn,-P_MMa/Po,'Color',color_MMa)
plot(M1_MMb/Mn,-P_MMb/Po,'Color',color_MMb)
plot(M1_MMc/Mn,-P_MMc/Po,'Color',color_MMc)
plot(M1_SOE/Mn,-P_SOE/Po,'Color',color_SOE)
plot(M1_SOI/Mn,-P_SOI/Po,'Color',color_SOI)

xlim([0 1.7])
ylim([0 0.5])
xlabel('Normalized Bending Moment (M/M_n)')
ylabel('Normalized Axial Compression (P/P_o)')

save_figure(hf,'example_figure_6','svg')

%% Error vs. Second-Order Inelastic
id2d_MMa = interactionDiagram2d(M1_MMa/Mn,-P_MMa/Po);
id2d_MMb = interactionDiagram2d(M1_MMb/Mn,-P_MMb/Po);
id2d_MMc = interactionDiagram2d(M1_MMc/Mn,-P_MMc/Po);
id2d_SOE = interactionDiagram2d(M1_SOE/Mn,-P_SOE/Po);
id2d_MMa_nl = interactionDiagram2d(M1_MMa_nl/Mn,-P_MMa_nl/Po);
id2d_MMb_nl = interactionDiagram2d(M1_MMb_nl/Mn,-P_MMb_nl/Po);
id2d_MMc_nl = interactionDiagram2d(M1_MMc_nl/Mn,-P_MMc_nl/Po);
id2d_SOE_nl = interactionDiagram2d(M1_SOE_nl/Mn,-P_SOE_nl/Po);
id2d_SOI = interactionDiagram2d(M1_SOI/Mn,-P_SOI/Po);
angles = linspace(0,pi/2,300);
errors_MMa = id2d_SOI.compareTwo(id2d_MMa,angles);
errors_MMb = id2d_SOI.compareTwo(id2d_MMb,angles);
errors_MMc = id2d_SOI.compareTwo(id2d_MMc,angles);
errors_SOE = id2d_SOI.compareTwo(id2d_SOE,angles);
errors_MMa_nl = id2d_SOI.compareTwo(id2d_MMa_nl,angles);
errors_MMb_nl = id2d_SOI.compareTwo(id2d_MMb_nl,angles);
errors_MMc_nl = id2d_SOI.compareTwo(id2d_MMc_nl,angles);
errors_SOE_nl = id2d_SOI.compareTwo(id2d_SOE_nl,angles);


hf = fs.figure(3.25,2.25);
ha = fs.axes([0.12 0.15 0.83 0.82]);

plot(rad2deg(angles),errors_MMa_nl,'--','Color',color_MMa)
plot(rad2deg(angles),errors_MMa,'Color',color_MMa)
plot(rad2deg(angles),errors_MMb_nl,'--','Color',color_MMb)
plot(rad2deg(angles),errors_MMb,'Color',color_MMb)
plot(rad2deg(angles),errors_MMc_nl,'--','Color',color_MMc)
plot(rad2deg(angles),errors_MMc,'Color',color_MMc)
plot(rad2deg(angles),errors_SOE_nl,'--','Color',color_SOE)
plot(rad2deg(angles),errors_SOE,'Color',color_SOE)

set(gca,'XTick',[0 15 30 45 60 75 90])
xlim([0 90])
%ylim([0 0.5])
xlabel('Angle in Normalized Interaction Diagram (deg)')
ylabel('Error')

save_figure(hf,'example_figure_7','svg')

%% Compare to results from main study
if compare_to_main_study
    i = 3852;

    study_path = fullfile(pathOf.BenchmarkStudyResults,'RC');
    study = BenchmarkStudy('open',study_path);
    fs2 = figureStyle('Display');

    load(study.path_of_results('ACI Interaction (a - no limit)'));
    hf = fs2.figure(5,5);
    ha = fs2.axes();
    plot(M_cs,-P_cs,'Color',color_CS)
    plot(M1_MMa_nl,-P_MMa_nl,'--','Color',color_MMa)
    plot(results(i).M1,-results(i).P1,'o-')
    plot(results(i).M2,-results(i).P2,'o-')
    xlabel(sprintf('Bending Moment (%s)',calcUnits.momentUnits))
    ylabel(sprintf('Axial Compression (%s)',calcUnits.forceUnits))

    load(study.path_of_results('ACI Interaction (a)'));
    hf = fs2.figure(5,5);
    ha = fs2.axes();
    plot(M_cs,-P_cs,'Color',color_CS)
    plot(M1_MMa,-P_MMa,'Color',color_MMa)
    plot(results(i).M1,-results(i).P1,'o--')
    plot(results(i).M2,-results(i).P2,'o--')
    xlabel(sprintf('Bending Moment (%s)',calcUnits.momentUnits))
    ylabel(sprintf('Axial Compression (%s)',calcUnits.forceUnits))

    load(study.path_of_results('ACI Interaction (c)'));
    hf = fs2.figure(5,5);
    ha = fs2.axes();
    plot(M_cs,-P_cs,'Color',color_CS)
    plot(M1_MMc,-P_MMc,'--','Color',color_MMa)
    plot(M2_MMc,-P_MMc,'-','Color',color_MMa)
    plot(results(i).M1,-results(i).P1,'o-')
    plot(results(i).M2,-results(i).P2,'o-')
    xlabel(sprintf('Bending Moment (%s)',calcUnits.momentUnits))
    ylabel(sprintf('Axial Compression (%s)',calcUnits.forceUnits))

    load(study.path_of_results('Analysis Interaction'));
    hf = fs2.figure(5,5);
    ha = fs2.axes();
    plot(M1_SOI,-P_SOI,'x-','Color',color_SOI)
    plot(M2_SOI,-P_SOI,'x-','Color',color_SOI)
    plot(results(i).M1,-results(i).P1,'o--')
    plot(results(i).M2,-results(i).P2,'o--')
    xlabel(sprintf('Bending Moment (%s)',calcUnits.momentUnits))
    ylabel(sprintf('Axial Compression (%s)',calcUnits.forceUnits))
end


