clear all; close all; clc;

fs = figureStyle('Display');

study_path = fullfile(pathOf.BenchmarkStudyResults,'RC');
study = BenchmarkStudy('open',study_path);
load(study.path_of_data,'data')

%% Select Data
%iData = 4;
iData = 1012;

%% Define SectionAnalysis Object
fs_options = struct;
fs_options.includePackageDefinition = true;
fs_options.conc_material = 'Concrete04';
fs_options.steel_material = 'ElasticPP';
[definition,numMat] = FiberSectionDefinition(data(iData),'3d',1,1,fs_options);
SA = SectionAnalysis(definition,1);
SA.fiber_discretization_method = 'legacy';

%% Table of section data
SA.printMaterialInfo();

%% Discretization
fs.figure(4,4)
fs.axes()
data(iData).section.plotSection();
SA.plotDiscretization();

%% Concrete Stress Strain
ss_results_tens_1 = SA.getStressStrain([0  0.005],1,'Steps',100);
ss_results_comp_1 = SA.getStressStrain([0 -0.02],1,'Steps',100);
ss_results_tens_2 = SA.getStressStrain([0  0.005],2,'Steps',100);
ss_results_comp_2 = SA.getStressStrain([0 -0.02],2,'Steps',100);

fs.figure(6,4)
fs.axes()
plot(ss_results_tens_1.disp,ss_results_tens_1.force)
plot(ss_results_comp_1.disp,ss_results_comp_1.force)
plot(ss_results_tens_2.disp,ss_results_tens_2.force)
plot(ss_results_comp_2.disp,ss_results_comp_2.force)
xlabel('Strain')
xlabel('Stress (ksi)')

%% Steel Stress Strain
ss_results_tens_3 = SA.getStressStrain([0  0.05],3,'Steps',100);
ss_results_comp_3 = SA.getStressStrain([0 -0.05],3,'Steps',100);

fs.figure(6,4)
fs.axes()
plot(ss_results_tens_3.disp,ss_results_tens_3.force)
plot(ss_results_comp_3.disp,ss_results_comp_3.force)
xlabel('Strain')
xlabel('Stress (ksi)')

