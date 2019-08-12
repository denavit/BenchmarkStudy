clear all; close all; clc;

%% Input
fc      = 4;
fy      = 60;
H       = 40;
B       = 40;
cover   = 0.15*H;
nbB     = 2;
nbH     = 5;
rhosr   = 0.08;
axis    = 'x';
units   = 'US';

num_points = 40;
fs = figureStyle('display');

% Define cross section
conc_cross_section = Rectangle_Shape(H,B);
Ab = H*B*rhosr/(2*nbB+2*nbH-4);
reinforcement = reinf_rect(B-2*cover,H-2*cover,0,0,nbB,nbH,Ab);
section = RC(fc,fy,conc_cross_section,reinforcement,units);
section.transverse_reinf_type = 'ties';

%% Plot Cross Section
fs.picture(6,6);
section.plotSection();
axis_limits('Margin',0.1)

%% Nonsway Frame Analysis
L       = 10*H;
beta    = 1;

BA = BenchmarkAnalysis2d_ACI_Nonsway_Frame(section,axis,L,beta);
% BA.second_order_moment_ratio_limit = [];

% Determine Design Interaction
BA.EIeff_type = 'a';
[P_a,M1_a,~,M2_a] = BA.design_interaction(num_points);

BA.EIeff_type = 'b';
[P_b,M1_b,~,M2_b] = BA.design_interaction(num_points);

BA.EIeff_type = 'c';
[P_c,M1_c,~,M2_c] = BA.design_interaction(num_points);

% Make Plot
hf = fs.figure(8,6);
ha = fs.axes([0.15 0.15 0.80 0.80]);
plot(BA.InteractionDiagram_M,-BA.InteractionDiagram_P)

plot(M1_a,-P_a,'r-')
plot(M2_a,-P_a,'r--')

plot(M1_b,-P_b,'k-')
plot(M2_b,-P_b,'k--')

plot(M1_c,-P_c,'g-')
plot(M2_c,-P_c,'g--')

legend('Cross Section',...
    'M_1 - EI_{eff} Type (a)','M_2 - EI_{eff} Type (a)',...
    'M_1 - EI_{eff} Type (b)','M_2 - EI_{eff} Type (b)',...
    'M_1 - EI_{eff} Type (c)','M_2 - EI_{eff} Type (c)',...
    'Location','NE')

axis_limits('Quadrant',1)
xlabel('Bending Moment (kip-in)')
ylabel('Axial Compression (kips)')


%% Sway Frame Analysis
L       = 10*H;
Ec      = section.Ec; % Use same material for column and beam
Igb     = 100000;
Lb      = 400;

BA = BenchmarkAnalysis2d_ACI_Sway_Frame(section,axis,L,Ec*Igb/Lb);
% BA.second_order_moment_ratio_limit = [];

% Determine Design Interaction
BA.EIeff_type = 'a';
[P_a,M1_a,~,M2_a] = BA.design_interaction(num_points);

BA.EIeff_type = 'b';
[P_b,M1_b,~,M2_b] = BA.design_interaction(num_points);

BA.EIeff_type = 'c';
[P_c,M1_c,~,M2_c] = BA.design_interaction(num_points);

% Make Plot
hf = fs.figure(8,6);
ha = fs.axes([0.15 0.15 0.80 0.80]);
plot(BA.InteractionDiagram_M,-BA.InteractionDiagram_P)

plot(M1_a,-P_a,'r-')
plot(M2_a,-P_a,'r--')

plot(M1_b,-P_b,'k-')
plot(M2_b,-P_b,'k--')

plot(M1_c,-P_c,'g-')
plot(M2_c,-P_c,'g--')

legend('Cross Section',...
    'M_1 - EI_{eff} Type (a)','M_2 - EI_{eff} Type (a)',...
    'M_1 - EI_{eff} Type (b)','M_2 - EI_{eff} Type (b)',...
    'M_1 - EI_{eff} Type (c)','M_2 - EI_{eff} Type (c)',...
    'Location','NE')

axis_limits('Quadrant',1)
xlabel('Bending Moment (kip-in)')
ylabel('Axial Compression (kips)')
