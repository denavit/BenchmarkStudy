clear all; close all; clc;

fs = figureStyle('Display');

study_path = fullfile(pathOf.BenchmarkStudyResults,'RC');
study = BenchmarkStudy('open',study_path);
load(study.path_of_data,'data')

%% Select Data
iData = 39;
%iData = 1012;

%% Moment Magnification
num_points = 40;

BA = BenchmarkAnalysis2d_ACI(data(iData));

BA.EIeff_type = 'a';
[P_a,M1_a,~,M2_a] = BA.design_interaction(num_points);

BA.EIeff_type = 'b';
[P_b,M1_b,~,M2_b] = BA.design_interaction(num_points);

BA.EIeff_type = 'c';
[P_c,M1_c,~,M2_c] = BA.design_interaction(num_points);

%% Second-Order Elastic Analysis
num_points = 40;

EI = 0.7*data(iData).section.Ec*data(iData).section.Ig(data(iData).axis);
elastic_frame = BenchmarkAnalysis2d_Elastic(data(iData),EI);
elastic_frame.includeInitialGeometricImperfections = false;

notionalLoadObject = notional_load(0.000,0.000,Inf);

[P1_SOEA,M1_SOEA,P2_SOEA,M2_SOEA] = elastic_frame.designInteraction(...
    data(iData).section,data(iData).axis,'ACI',num_points,...
    notionalLoadObject,'zero',1.0,1.0,'none');

P_at_PeakMomentRatio = elastic_frame.determineLoadForMomentRatio(1.4,1.0,1.0,'none',nan);
P1_SOEA = max(P1_SOEA,-P_at_PeakMomentRatio);
P2_SOEA = max(P2_SOEA,-P_at_PeakMomentRatio);

%% Second-Order Inelastic Analysis (OpenSees)


%% Make Plot
my_colors = lines(5);

hf = fs.figure(8,6);
ha = fs.axes([0.10 0.10 0.85 0.85]);
plot(BA.InteractionDiagram_M,-BA.InteractionDiagram_P,'Color',my_colors(1,:))

plot(M1_a,-P_a, '-','Color',my_colors(2,:))
plot(M2_a,-P_a,'--','Color',my_colors(2,:))

plot(M1_b,-P_b, '-','Color',my_colors(3,:))
plot(M2_b,-P_b,'--','Color',my_colors(3,:))

plot(M1_c,-P_c, '-','Color',my_colors(4,:))
plot(M2_c,-P_c,'--','Color',my_colors(4,:))

plot(M1_SOEA,-P1_SOEA, '-','Color',my_colors(5,:))
plot(M2_SOEA,-P2_SOEA,'--','Color',my_colors(5,:))

legend('Cross Section',...
    'M_1 - EI_{eff} Type (a)','M_2 - EI_{eff} Type (a)',...
    'M_1 - EI_{eff} Type (b)','M_2 - EI_{eff} Type (b)',...
    'M_1 - EI_{eff} Type (c)','M_2 - EI_{eff} Type (c)',...
    'M_1 - 2nd-Order','M_2 - 2nd-Order',...
    'Location','NE')
legend('BoxOff')

axis_limits('Quadrant',1)
xlabel('Bending Moment (kip-in)')
ylabel('Axial Compression (kips)')

