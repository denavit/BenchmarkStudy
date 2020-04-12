function runBenchmarkStudy_ACI_Section(study,tag,option)

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
%numPoints = study.numPoints;

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
fprintf('\nACI Interaction (Section) Study (%s)\n',study.study_name);
study_name = strrep(study.study_name,'_','\_');
studystr = sprintf('ACI Interaction (Section) Study (%s)',study_name);
hwait = waitbar(0,studystr);
for iData = selectedData
    str = {'',studystr,sprintf('Case %i of %i',iData,numData)};
    waitbar(iData/numData,hwait,str);

        
    num_points = 50;
    sc = data(iData).section.strainCompatibilityAciObject;
    switch lower(data(iData).axis)
        case {'x','z'}
            [P,M,~,et] = sc.interactionSweep(0,num_points);
        case 'y'
            [P,~,M,et] = sc.interactionSweep(pi/2,num_points);
        otherwise
            error('Bad axis');
    end
    
    switch option
        case 'nominal'
            
        case 'reduced'
            phi = obj.section.phi(et);
            P = phi.*P;
            M = phi.*M;
            
        otherwise
            error('Unknown option: %s',option);
    end
    
    % Store Analysis Results
    results(iData).study_type = 'ACI Interaction (Section)';
    results(iData).P1               = P;
    results(iData).M1               = M;
    results(iData).P2               = P;
    results(iData).M2               = M;
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