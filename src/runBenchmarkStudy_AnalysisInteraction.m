function runBenchmarkStudy_AnalysisInteraction(study,tag,fiberSectionDefinitionOptions,analysisOptions)

%% Load Data
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

% Set the working path
analysisOptions.scratchPath = study.scratch_path;

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
fprintf('\nAnalysis Interaction Study (%s)\n',study.study_name);
study_name = strrep(study.study_name,'_','\_');
studystr = sprintf('Analysis Interaction Study (%s)',study_name);
hwait = waitbar(0,studystr);
set(hwait,'Position',[100 100 275 75])
for iData = selectedData
    str = {'',studystr,sprintf('Case %i of %i',iData,numData)};
    waitbar(iData/numData,hwait,str);
    
    % Get cross section definition
    sectionDef = FiberSectionDefinition(data(iData).section,data(iData).axis,1,1,fiberSectionDefinitionOptions);
    
    % Create Analysis Object
    ba = BenchmarkAnalysis2d_OpenSees(data(iData),sectionDef,analysisOptions);
    % ba.base_tolerance_force = 0.00001; % Adjust force tolerance if necessary during re-run

    % Initilize results
    limit_type  = cell(numPoints,1);
    M1          = nan(numPoints,1);
    P1          = nan(numPoints,1);
    M2          = nan(numPoints,1);
    P2          = nan(numPoints,1);
    def         = nan(numPoints,1);

    % Set Data
    results(iData).study_type = 'AnalysisInteraction';

    
    % Run axial only analysis to get Pn
    try
        iResults = ba.runAnalysis('LimitPoint_Proportional',0,[],1);

        if ~iResults.limitPoint.good
            iResults = ba.runAnalysis('LimitPoint_Proportional',0,[],2);
        end

        if ~iResults.limitPoint.good
            fprintf('Limit Point Not Obtained: Case %i, Axial Only\n',iData);
        end

        limit_type{1} = iResults.limitPoint.limit_type;
        M1(1)         = iResults.limitPoint.M1;
        P1(1)         = iResults.limitPoint.P1;
        M2(1)         = iResults.limitPoint.M2;
        P2(1)         = iResults.limitPoint.P2;
        def(1)        = iResults.limitPoint.def;
    catch err
        fprintf('%s\n',err.message)
        fprintf('Axial Only Analysis Failed: Case %i\n',iData);
        continue
    end


    % Run nonproportional analyses at different axial loads
    Ps = linspace(P1(1),0,numPoints);
    for i = 2:length(Ps)
        P = Ps(i);
        if P ~= 0
            try
                iResults = ba.runAnalysis('LimitPoint_NonProportional',P,[],1);

                if ~iResults.limitPoint.good
                    iResults = ba.runAnalysis('LimitPoint_NonProportional',P,[],2);
                end

                if ~iResults.limitPoint.good
                    iResults = ba.runAnalysis('LimitPoint_NonProportional',P,[],3);
                end

                if ~iResults.limitPoint.good
                    fprintf('Limit Point Not Obtained: Case %i, Axial Load Level %i (%g)\n',iData,i,P);
                end

                limit_type{i} = iResults.limitPoint.limit_type;
                M1(i)         = iResults.limitPoint.M1;
                P1(i)         = iResults.limitPoint.P1;
                M2(i)         = iResults.limitPoint.M2;
                P2(i)         = iResults.limitPoint.P2;
                def(i)        = iResults.limitPoint.def;
            catch err
                fprintf('%s\n',err.message)
                fprintf('Non-Proportional Analysis Failed: Case %i, Axial Load Level %i (%g)\n',iData,i,P);
                continue
            end
        else
            % Run cross section analysis for P=0
            try
                iResults = ba.runSectionAnalysis('LimitPoint_NonProportional',0,[],1);
                
                limit_type{i} = iResults.limitPoint.limit_type;
                M1(i)         = iResults.limitPoint.M1;
                P1(i)         = iResults.limitPoint.P1;
                M2(i)         = iResults.limitPoint.M2;
                P2(i)         = iResults.limitPoint.P2;
                def(i)        = NaN;
            catch err
                fprintf('%s\n',err.message)
                fprintf('Cross Section Analysis Failed: Case: %i\n',iData);
                continue
            end
        end
    end

    % Store Analysis Results
    results(iData).limit_type       = limit_type;
    results(iData).P1               = P1;
    results(iData).M1               = M1;
    results(iData).P2               = P2;
    results(iData).M2               = M2;
    results(iData).def              = def;
    
    % Store Extra Data
    results(iData).Pn_ana           = results(iData).P1(1);
    results(iData).Pno              = data(iData).section.Pnco;
    results(iData).Pn               = data(iData).section.Pnc(data(iData).axis);
    results(iData).Mno              = data(iData).section.Mno(data(iData).axis);
    results(iData).EIeff            = data(iData).section.EI(data(iData).axis,'ColumnStrength');
    results(iData).EIgross          = data(iData).section.EI(data(iData).axis,'Gross');
    results(iData).EsIs             = data(iData).section.getSectionData('GrossSteelFlexuralRigidity',data(iData).axis);
    results(iData).EsIsr            = data(iData).section.getSectionData('GrossReinforcingFlexuralRigidity',data(iData).axis);
    results(iData).EcIc             = data(iData).section.getSectionData('GrossConcreteFlexuralRigidity',data(iData).axis);
    results(iData).AsFy             = data(iData).section.getSectionData('AsFy',data(iData).axis);
    results(iData).Acfc             = data(iData).section.getSectionData('Acfc',data(iData).axis);
    results(iData).Fy               = data(iData).section.getSectionData('SteelStrength');
    results(iData).fc               = data(iData).section.getSectionData('ConcreteStrength');
    results(iData).rhos             = data(iData).section.getSectionData('SteelRatio');
    results(iData).rhosr            = data(iData).section.getSectionData('ReinforcingRatio');
    results(iData).ShapeFactor      = data(iData).section.getSectionData('ShapeFactor',data(iData).axis);
    
    PnOverPo = -results(iData).P1(1)/results(iData).Pno;
    results(iData).lambdaoe_equiv   = sqrt(AISC_inverse_column_curve(PnOverPo));
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