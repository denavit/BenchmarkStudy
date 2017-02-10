clear all; close all; clc;

data = struct;
data.EI      = 1000;
data.EA      = 100000;
data.L       = 100;
data.axis    = 'strong';
data.section = general_elastic_section(1,data.EA,data.EI);
data.type    = 'Elastic';

% data.frame_type = 'Sidesway_Inhibited';
% data.beta    = 1;
% data.delta0  = data.L/100;

data.frame_type = 'Sidesway_Uninhibited';
data.kqtop   = 30*data.EI/data.L; %Inf;
data.kqbot   = 50*data.EI/data.L; %Inf;
data.gamma   = 2;
data.delta0  = data.L/1000;
data.Delta0  = data.L/500;

%%
analysisOptions = struct;
analysisOptions.numElements                             = 20;
analysisOptions.includeInitialGeometricImperfections    = true;
analysisOptions.extraResultsOutput                      = true;
analysisOptions.deleteFilesAfterAnalysis                = true;
analysisOptions.scratchPath                             = '.';


%% Elastic Analysis
BA_Elastic = BenchmarkAnalysis2d_Elastic(data);

Pcr = BA_Elastic.eulerLoad();

P = -0.3*Pcr;
load = 0.005; % M or H

x = linspace(0,data.L,100);
BA_Elastic_V1      = BA_Elastic.firstOrderDisplacement(P,load,x);
BA_Elastic_V2      = BA_Elastic.secondOrderDisplacement(P,load,x);
BA_Elastic_M1      = BA_Elastic.firstOrderMoment(P,load,x); 
BA_Elastic_M2      = BA_Elastic.secondOrderMoment(P,load,x);

BA_Elastic.includeInitialGeometricImperfections = false;
BA_Elastic_V1_woii = BA_Elastic.firstOrderDisplacement(P,load,x);
BA_Elastic_V2_woii = BA_Elastic.secondOrderDisplacement(P,load,x);
BA_Elastic_M1_woii = BA_Elastic.firstOrderMoment(P,load,x);
BA_Elastic_M2_woii = BA_Elastic.secondOrderMoment(P,load,x);


%%
sectionDef = FiberSectionDefinition(data.section,'2dStrong',1,1,struct);

BA_OpenSees = BenchmarkAnalysis2d_OpenSees(data,sectionDef,analysisOptions);

BA_OpenSees.geomTransfType = 'Corotational';
BA_OpenSees.includeInitialGeometricImperfections = true;
results2      = BA_OpenSees.runAnalysis('TargetForce_NonProportional',P,load,2);

if ~strcmp(results2.exitStatus,'Loads Reached')
    fprintf(results2.textOutput)
end

BA_OpenSees.geomTransfType = 'Corotational';
BA_OpenSees.includeInitialGeometricImperfections = false;
results2_woii = BA_OpenSees.runAnalysis('TargetForce_NonProportional',P,load,2);

if ~strcmp(results2_woii.exitStatus,'Loads Reached')
    fprintf(results2_woii.textOutput)
end

BA_OpenSees.geomTransfType = 'Linear';
BA_OpenSees.includeInitialGeometricImperfections = true;
results1      = BA_OpenSees.runAnalysis('TargetForce_NonProportional',P,load,2);

if ~strcmp(results1.exitStatus,'Loads Reached')
    fprintf(results1.textOutput)
end

BA_OpenSees.geomTransfType = 'Linear';
BA_OpenSees.includeInitialGeometricImperfections = false;
results1_woii = BA_OpenSees.runAnalysis('TargetForce_NonProportional',P,load,2);

if ~strcmp(results1_woii.exitStatus,'Loads Reached')
    fprintf(results1_woii.textOutput)
end


%%
my_colors = lines(2);

figure
hold all
plot(BA_Elastic_V1     ,x, '-','Color',my_colors(1,:))
plot(BA_Elastic_V1_woii,x,'--','Color',my_colors(1,:))
plot(BA_Elastic_V2     ,x, '-','Color',my_colors(2,:))
plot(BA_Elastic_V2_woii,x,'--','Color',my_colors(2,:))
plot(     results1.path.def(end)*[1 1],[0 data.L], '-','Color',my_colors(1,:))
plot(results1_woii.path.def(end)*[1 1],[0 data.L],'--','Color',my_colors(1,:))
plot(     results2.path.def(end)*[1 1],[0 data.L], '-','Color',my_colors(2,:))
plot(results2_woii.path.def(end)*[1 1],[0 data.L],'--','Color',my_colors(2,:))

xlabel('Lateral Displacement') 
axis_limits('Quadrant',1)
legend('V1 (w/ ii)','V1 (w/o ii)','V2 (w/ ii)','V2 (w/o ii)',...
    'Location','EastOutside')


figure
hold all
plot(BA_Elastic_M1     ,x, '-','Color',my_colors(1,:))
plot(BA_Elastic_M1_woii,x,'--','Color',my_colors(1,:))
plot(BA_Elastic_M2     ,x, '-','Color',my_colors(2,:))
plot(BA_Elastic_M2_woii,x,'--','Color',my_colors(2,:))
plot(     -results1.distE.M(end,:),     results1.distE.xi(end,:),'^','Color',my_colors(1,:))
plot(-results1_woii.distE.M(end,:),results1_woii.distE.xi(end,:),'o','Color',my_colors(1,:))
plot(     -results2.distE.M(end,:),     results2.distE.xi(end,:),'^','Color',my_colors(2,:))
plot(-results2_woii.distE.M(end,:),results2_woii.distE.xi(end,:),'o','Color',my_colors(2,:))

xlabel('Bendind Moment') 
axis_limits('BalanceX')
legend('M1 (w/ ii)','M1 (w/o ii)','M2 (w/ ii)','M2 (w/o ii)',...
    'M1 (w/ ii)','M1 (w/o ii)','M2 (w/ ii)','M2 (w/o ii)',...
    'Location','EastOutside')



%% Critical Load Check

load_path = 0.001;
P_path       = linspace(0,-1.1*Pcr,200);
M1_path      = nan(size(P_path));
M2_path      = nan(size(P_path));
M1_path_woii = nan(size(P_path));
M2_path_woii = nan(size(P_path));
for i = 1:length(P_path)
    BA_Elastic.includeInitialGeometricImperfections = true;
    M1_path(i) = BA_Elastic.maxFirstOrderMoment(P_path(i),load_path);
    M2_path(i) = BA_Elastic.maxSecondOrderMoment(P_path(i),load_path);
    
    BA_Elastic.includeInitialGeometricImperfections = false;
    M1_path_woii(i) = BA_Elastic.maxFirstOrderMoment(P_path(i),load_path);
    M2_path_woii(i) = BA_Elastic.maxSecondOrderMoment(P_path(i),load_path);    
end

figure
hold all
xlim([0 5])
plot(M1_path     ,-P_path, '-','Color',my_colors(1,:))
plot(M1_path_woii,-P_path,'--','Color',my_colors(1,:))
plot(M2_path     ,-P_path, '-','Color',my_colors(2,:))
plot(M2_path_woii,-P_path,'--','Color',my_colors(2,:))
plot(xlim,Pcr*[1 1],'k--')

ylabel('Axial Compression')
xlabel('Bendind Moment')
axis_limits('Quadrant',1)
legend('M1 (w/ ii)','M1 (w/o ii)','M2 (w/ ii)','M2 (w/o ii)',...
    'Location','EastOutside')