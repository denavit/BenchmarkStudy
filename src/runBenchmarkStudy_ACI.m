function runBenchmarkStudy_ACI(study,tag,EIeff_type,include_strength_reduction,include_stiffness_reduction,second_order_moment_ratio_limit)

%% Load Data
% Load Sections
if iscell(study)
    selectedData = study{2};
    study = study{1};
    load(study.path_of_data);
    numData = length(data);
    newStudy = false;
else
    load(study.path_of_data);
    numData = length(data);
    selectedData = 1:numData;
    newStudy = true;    
end

% Number of Points in Interaction
numPoints = study.numPoints;

%% Initilize Results Structure
if newStudy
    results(numData) = struct;
    if ~isempty(study.check_results_tag(tag))
        error('Results file already exists, will not be able to save results')
    end    
else
    load(study.path_of_results(tag),'results');
end

%% Run Study
fprintf('\nACI Interaction Study (%s)\n',study.study_name);
study_name = strrep(study.study_name,'_','\_');
studystr = sprintf('ACI Interaction Study (%s)',study_name);
hwait = waitbar(0,studystr);
for iData = selectedData
    str = {'',studystr,sprintf('Case %i of %i',iData,numData)};
    waitbar(iData/numData,hwait,str);

    % ACI Analysis Object
    BA = BenchmarkAnalysis2d_ACI(data(iData),EIeff_type);
    BA.include_strength_reduction = include_strength_reduction;
    BA.include_stiffness_reduction = include_stiffness_reduction;
    BA.second_order_moment_ratio_limit = second_order_moment_ratio_limit;
    
    % Run Analysis 
    [P1,M1,P2,M2] = BA.design_interaction(numPoints);
    
    % Store Analysi Results
    results(iData).study_type = 'ACI Interaction';
    results(iData).P1               = P1;
    results(iData).M1               = M1;
    results(iData).P2               = P2;
    results(iData).M2               = M2;
end
fprintf('  analyses complete\n');
close(hwait);

%% Save Data
if ~newStudy
    study.remove_results_file(tag);
end
path = study.add_results_file(tag);
save(path,'results','-v7.3');
fprintf('  results saved\n');

end