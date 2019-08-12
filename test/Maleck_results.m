clear all; close all; clc;

study_name = 'Maleck_SA';
%study_name = 'Maleck_WA';

print_to_file = true;

% Load data and results
study_path = fullfile(pathOf.BenchmarkStudyResults,study_name);
study = BenchmarkStudy('open',study_path);
load(study.path_of_data);
numData = length(data);
load(study.path_of_results('Analysis Interaction'));
results_FN = results;
load(study.path_of_results('Design Interaction (Maleck (EL))'));
results_EL = results;
load(study.path_of_results('Design Interaction (Maleck (DA))'));
results_DA = results;
clear results;

% Print table of results
if print_to_file
    fid = fopen(sprintf('%s_results.txt',study_name),'w');
else
    fid = 1; 
end
fprintf(fid,'                             EL             DA    \n');
fprintf(fid,'  Frame Name   lambda     e-     e+      e-     e+\n');
for iData = 1:numData
    lambda = sqrt(AISC_inverse_column_curve(-results_FN(iData).Pn_ana/results_FN(iData).Pno));
    EL_en  = 100*min([0 min(results_EL(iData).errors)]);
    EL_ep  = 100*max([0 max(results_EL(iData).errors)]);
    DA_en  = 100*min([0 min(results_DA(iData).errors)]);
    DA_ep  = 100*max([0 max(results_DA(iData).errors)]);
    fprintf(fid,'%13s %7.3f %6.1f %6.1f  %6.1f %6.1f\n',data(iData).frame_name,lambda,EL_en,EL_ep,DA_en,DA_ep);
end
if print_to_file
    fclose(fid);
end
