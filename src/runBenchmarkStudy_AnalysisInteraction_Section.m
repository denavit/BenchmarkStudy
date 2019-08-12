function runBenchmarkStudy_AnalysisInteraction_Section(study,tag,fiberSectionDefinitionOptions,analysisOptions)

%% Load Data
% Load Sections
load(study.path_of_data);
numData = length(data);

% Number of Points in Interaction
numPoints = study.numPoints;

% Set the working path
analysisOptions.scratchPath = study.scratch_path;

%% Initilize Results Structure
if iscell(tag)
    newStudy = isempty(study.check_results_tag(tag{1}));
    selectedData = tag{2};
    tag = tag{1};
else
    newStudy = true;
    selectedData = 1:numData;
    assert(isempty(study.check_results_tag(tag)),'tag is not unique')    
end

if newStudy
    results(numData) = struct;
else
    load(study.path_of_results(tag));
end

%% Run Study
fprintf('\nSection Analysis Interaction Study (%s)\n',study.study_name);
study_name = strrep(study.study_name,'_','\_');
studystr = sprintf('Section Analysis Interaction Study (%s)',study_name);
hwait = waitbar(0,studystr);
for iData = selectedData
    str = {studystr,sprintf('Case %i of %i',iData,numData)};
    waitbar(iData/numData,hwait,str);

    % Get cross section definition
    if length(fiberSectionDefinitionOptions) == 1
        fsDefOpts = fiberSectionDefinitionOptions;
    elseif length(fiberSectionDefinitionOptions) == numData
        fsDefOpts = fiberSectionDefinitionOptions(iData);
    else
        error('Bad size for fiberSectionDefinitionOptions')
    end
    
    sectionDef = FiberSectionDefinition(data(iData).section,data(iData).axis,1,1,fsDefOpts);
    
    % Create Analysis Object
    ba = BenchmarkAnalysis2d_OpenSees(data(iData),sectionDef,analysisOptions);

    % Initilize results
    limit_type = cell(numPoints,1);
    M1         = nan(numPoints,1);
    P1         = nan(numPoints,1);    

    % Set Data
    results(iData).study_type = 'AnalysisInteraction_Section';
    
    % Run axial only analysis to get Pn
    iResults = ba.runSectionAnalysis('LimitPoint_Axial',0,[],1);
    limit_type{1} = iResults.limitPoint.limit_type;
    M1(1)         = iResults.limitPoint.M1;
    P1(1)         = iResults.limitPoint.P1;
       
    % Run nonproportional analyses at different axial loads
    Ps = linspace(P1(1),0,numPoints);
    Ps = [Ps(1) 0.995*Ps(1) Ps(2:end)];
    for i = 2:length(Ps)
        iP = Ps(i);
        iResults = ba.runSectionAnalysis('LimitPoint_NonProportional',iP,[],1);
        limit_type{i} = iResults.limitPoint.limit_type;
        M1(i)         = iResults.limitPoint.M1;
        P1(i)         = iResults.limitPoint.P1;
    end
    
    % Store Results
    results(iData).limit_type = limit_type;
    results(iData).P1         = P1;
    results(iData).M1         = M1;
    
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









