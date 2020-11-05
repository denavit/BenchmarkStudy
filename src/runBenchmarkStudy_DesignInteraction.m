function runBenchmarkStudy_DesignInteraction(study,tag,tag_Analysis,option)

%% Basic Information
angles = linspace(0,pi/2,300)';

%% Load Data
if iscell(study)
    selectedData = study{2};
    study = study{1};
    load(study.path_of_data,'data');
    numData = length(data);
    newStudy = false;
else
    load(study.path_of_data,'data');
    numData = length(data);
    selectedData = 1:numData;
    newStudy = true;    
end

% Load Results from Analysis Interaction Study
load(study.path_of_results(tag_Analysis),'results');
results_FN = results;
clear results;

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
fprintf('\nDesign Interaction Study (%s)\n',study.study_name);
study_name = strrep(study.study_name,'_','\_');
studystr = sprintf('Design Interaction Study (%s)',study_name);
hwait = waitbar(0,studystr);
set(hwait,'Position',[100 100 275 75])
for iData = selectedData
    str = {studystr,sprintf('Case %i of %i',iData,numData)};
    waitbar(iData/numData,hwait,str);

    results(iData).study_type = 'DesignInteraction';
    
    % Retreive Fully Nonlinear Results
    P1_FN = results_FN(iData).P1;
    M1_FN = results_FN(iData).M1;
    if abs(P1_FN(end)) < 1e-3
        P1_FN(end) = 0;
    end

    Py_tau = data(iData).section.Pnco;

    %%%%%%%% Construct Interaction Diagrams %%%%%%%%
    switch option
        case 'AISC 2016 (DA)'
            notionalLoadObject          = notional_load(0.000,0.002,1.7);
            effectiveLengthFactorType   = 'one';
            columnStiffnessReduction    = 0.8;
            beamStiffnessReduction      = 0.8; 
            tauType                     = 'AISC';
            peakMomentRatio             = [];
            neglectInitialImperf        = true;
            elasticStiffnessType        = 'ColumnStrength';
            designStrengthType          = 'AISC/in_plane';
            % Modifications for composite
            if data(iData).section.hasConcrete
                tauType = 'Composite';
                data(iData).section.option_EI = 'AISC2016';
            end
           
        case 'AISC 2016 (DA) - Trial-ACDB'
            assert(data(iData).section.hasConcrete,'Section should be composite')
            data(iData).section.option_EI = 'AISC2016';
            
            notionalLoadObject          = notional_load(0.000,0.002,1.7);
            effectiveLengthFactorType   = 'one';
            columnStiffnessReduction    = 0.8;
            beamStiffnessReduction      = 0.8; 
            tauType                     = 'Composite';
            peakMomentRatio             = [];
            neglectInitialImperf        = true;
            elasticStiffnessType        = 'ColumnStrength';
            designStrengthType          = 'Trial-ACDB';         
            
        case 'AISC 2016 (DMMI)'
            notionalLoadObject          = notional_load(0.000,0.000,Inf);
            effectiveLengthFactorType   = 'zero';
            columnStiffnessReduction    = 0.8;
            beamStiffnessReduction      = 0.8; 
            tauType                     = 'AISC';
            peakMomentRatio             = [];
            neglectInitialImperf        = false;
            elasticStiffnessType        = 'ColumnStrength';
            designStrengthType          = 'AISC/in_plane';
            % Modifications for composite
            if data(iData).section.hasConcrete
                tauType = 'Composite';
                data(iData).section.option_EI = 'AISC2016';
            end
            
        case 'AISC 2016 (EL)'
            notionalLoadObject          = notional_load(0.000,0.002,1.7);
            effectiveLengthFactorType   = 'StoryBasedK';
            columnStiffnessReduction    = 1.0;
            beamStiffnessReduction      = 1.0; 
            tauType                     = 'none';
            peakMomentRatio             = [];
            neglectInitialImperf        = true;
            elasticStiffnessType        = 'ColumnStrength';
            designStrengthType          = 'AISC/in_plane';
            % Modifications for composite
            if data(iData).section.hasConcrete
                tauType = 'Composite';
                data(iData).section.option_EI = 'AISC2016';
            end
            
        case 'ACI 2019'
            notionalLoadObject          = notional_load(0.000,0.000,Inf);
            effectiveLengthFactorType   = 'zero';
            columnStiffnessReduction    = 1.0;
            beamStiffnessReduction      = 1.0;
            tauType                     = 'none';
            peakMomentRatio             = 1.4;
            neglectInitialImperf        = true;
            elasticStiffnessType        = '0.7EcIg';
            designStrengthType          = 'ACI';
            
        case 'ACI 2019 (No Moment Ratio Limit)'
            notionalLoadObject          = notional_load(0.000,0.000,Inf);
            effectiveLengthFactorType   = 'zero';
            columnStiffnessReduction    = 1.0;
            beamStiffnessReduction      = 1.0;
            tauType                     = 'none';
            peakMomentRatio             = [];                
            neglectInitialImperf        = true;
            elasticStiffnessType        = '0.7EcIg';
            designStrengthType          = 'ACI';
            
        case 'Maleck (DA)'
            switch data(iData).axis
                case 'strong'
                    stiffnessReduction = 0.9;
                case 'weak'
                    stiffnessReduction = 0.8;
                otherwise
                    error('Bad axis')
            end
            notionalLoadObject          = notional_load(0.002,0.000,Inf);
            effectiveLengthFactorType   = 'one';
            columnStiffnessReduction    = stiffnessReduction;
            beamStiffnessReduction      = stiffnessReduction;
            tauType                     = 'maleck';
            peakMomentRatio             = [];
            neglectInitialImperf        = true;
            elasticStiffnessType        = 'ColumnStrength';
            designStrengthType          = 'AISC';
            
        case 'Maleck (EL)'
            notionalLoadObject          = notional_load(0.000,0.000,Inf);
            effectiveLengthFactorType   = 'StoryBasedK';
            columnStiffnessReduction    = 1.0;
            beamStiffnessReduction      = 1.0;
            tauType                     = 'none';
            peakMomentRatio             = [];
            neglectInitialImperf        = true;
            elasticStiffnessType        = 'ColumnStrength';
            designStrengthType          = 'AISC';
          
        case 'Scratch'
            tauType                     = 'none';
            notionalLoadObject          = notional_load(0.000,0.002,1.7);
            effectiveLengthFactorType   = 'one';
            columnStiffnessReduction    = 0.8;
            beamStiffnessReduction      = 0.8;
            peakMomentRatio             = [];
            neglectInitialImperf        = true;
            elasticStiffnessType        = 'ColumnStrength';
            designStrengthType          = 'AISC';            
            
        otherwise
            error('Unknown Design Type')
    end

    % Get Elastic Frame
    data(iData).EI = data(iData).section.EI(data(iData).axis,elasticStiffnessType);
    if neglectInitialImperf 
        data(iData).delta0 = 0;
        data(iData).Delta0 = 0;
    end
    BA_Elastic = BenchmarkAnalysis2d_Elastic(data(iData));    
    if neglectInitialImperf 
        BA_Elastic.includeInitialGeometricImperfections = false;
    end
    
    % Store Interaction Results
    [P1,M1,P2,M2] = BA_Elastic.designInteraction(...
        data(iData).section,data(iData).axis,designStrengthType,numPoints,...
        notionalLoadObject,effectiveLengthFactorType,...
        columnStiffnessReduction,beamStiffnessReduction,tauType);

    % Moment Ratios
    if ~isempty(peakMomentRatio)
        loadAtPeakMomentRatio = BA_Elastic.determineLoadForMomentRatio(peakMomentRatio,...
            columnStiffnessReduction,beamStiffnessReduction,tauType,Py_tau);

        P1 = max(P1,-loadAtPeakMomentRatio);
        P2 = max(P2,-loadAtPeakMomentRatio);
    end        

    results(iData).P1 = P1;
    results(iData).M1 = M1;
    results(iData).P2 = P2;
    results(iData).M2 = M2; 

    % Pmax with Phi
    [Pmax,M2atPmax] = BA_Elastic.maximumLoad(...
        data(iData).section,data(iData).axis,designStrengthType,...
        notionalLoadObject,...
        effectiveLengthFactorType,...            
        columnStiffnessReduction,beamStiffnessReduction,tauType);

    results(iData).PmaxWithPhi        = -Pmax;        
    results(iData).M2atPmaxWithPhi    = M2atPmax; 

    % Store Implied Interaction Diagram 
    [P2e,M2e] = BA_Elastic.impliedInteraction(P1_FN,M1_FN,...
        notionalLoadObject,...
        columnStiffnessReduction,beamStiffnessReduction,tauType,Py_tau);        

    results(iData).P2e = P2e;
    results(iData).M2e = M2e;        

    % Axial Strength 
    sectionKL = data(iData).section;
    sectionKL.Lx = BA_Elastic.L;
    sectionKL.Ly = BA_Elastic.L;
    sectionKL.Kx = BA_Elastic.K;
    sectionKL.Ky = BA_Elastic.K;
    results(iData).PnK = sectionKL.Pnc(data(iData).axis);


    %%%%%%%% Compute Errors %%%%%%%%        
    % Py = sections(iSection).section.Pnco;
    % My = sections(iSection).section.Mn(sections(iSection).axis);
    Py = P1_FN(1);
    My = M1_FN(end);        
    base_id = interactionDiagram2d(M1_FN/My,P1_FN/Py);
    test_id = interactionDiagram2d(M1/My,P1/Py);
    errors  = base_id.compareTwo(test_id,angles);

    % Store Error Results
    results(iData).errors = errors;
    results(iData).angles = angles;


    %%%%%%%% Compute Key Angles %%%%%%%%
    FN_id = interactionDiagram2d(M1_FN,P1_FN);
    % Determine angle for high axial
    P1_HP = 0.8*P1_FN(1);
    M1_HP = FN_id.findXgivenY(P1_HP,'pos');
    if isempty(M1_HP)
        angle_HP = NaN;
    else
        angle_HP = atan((P1_HP/Py)/(M1_HP/My));
    end

    % Determine angle for low axial
    P1_LP = 0.2*P1_FN(1);
    M1_LP = FN_id.findXgivenY(P1_LP,'pos');
    if isempty(M1_LP)
        angle_LP = NaN;
    else
        angle_LP = atan((P1_LP/Py)/(M1_LP/My));
    end

    % Determine angle for high moment
    M1_HM = 0.8*M1_FN(end);
    P1_HM = FN_id.findYgivenX(M1_HM,'neg');
    if isempty(P1_HM)
        angle_HM = NaN;
    else
        angle_HM = atan((P1_HM/Py)/(M1_HM/My));
    end

    % Determine angle for low moment
    M1_LM = 0.2*M1_FN(end);
    P1_LM = FN_id.findYgivenX(M1_LM,'neg');
    if isempty(P1_LM)
        angle_LM = NaN;
    else
        angle_LM = atan((P1_LM/Py)/(M1_LM/My));
    end

    % Store Results
    results(iData).angle_HP = angle_HP;
    results(iData).angle_LP = angle_LP;
    results(iData).angle_HM = angle_HM;
    results(iData).angle_LM = angle_LM;


    %%%%%%%% Compute Drifts and Moment Ratios for Each Angle %%%%%%%%    
    id    = interactionDiagram2d(M1/My,-P1/Py);
    radii = id.radial_distance(angles);
    iM = My*cos(angles).*radii;
    iP = Py*sin(angles).*radii;

    [drift_s1,drift_s2,~]           = BA_Elastic.computeDrifts(iP/1.6,iM/1.6);
    [drift_u1,drift_u2,driftRatios] = BA_Elastic.computeDrifts(iP,iM);

    % Store Results
    results(iData).drift_s1       = drift_s1;
    results(iData).drift_s2       = drift_s2;
    results(iData).drift_u1       = drift_u1;
    results(iData).drift_u2       = drift_u2;
    results(iData).driftRatios    = driftRatios;
    
    
    % Computed Data
    results(iData).errors       = results(iData).errors;
    results(iData).angles       = results(iData).angles;

    results(iData).drift_s1     = results(iData).drift_s1;
    results(iData).drift_s2     = results(iData).drift_s2;
    results(iData).drift_u1     = results(iData).drift_u1;
    results(iData).drift_u2     = results(iData).drift_u2;
    results(iData).driftRatios  = results(iData).driftRatios;

    results(iData).highAxial    = results(iData).angles >= results(iData).angle_HP;
    results(iData).lowAxial     = results(iData).angles <= results(iData).angle_LP;
    results(iData).highMoment   = results(iData).angles <= results(iData).angle_HM;
    results(iData).lowMoment    = results(iData).angles >= results(iData).angle_LM;

    results(iData).Pno          = data(iData).section.Pnco;
    results(iData).Pn           = data(iData).section.Pnc(data(iData).axis);
    results(iData).Mno          = data(iData).section.Mno(data(iData).axis);
    results(iData).EIeff        = data(iData).section.EI(data(iData).axis,'ColumnStrength');
    results(iData).EIgross      = data(iData).section.EI(data(iData).axis,'Gross');
    results(iData).EsIs         = data(iData).section.getSectionData('GrossSteelFlexuralRigidity',data(iData).axis);
    results(iData).EcIc         = data(iData).section.getSectionData('GrossConcreteFlexuralRigidity',data(iData).axis);
    results(iData).Fy           = data(iData).section.getSectionData('SteelStrength');
    results(iData).fc           = data(iData).section.getSectionData('ConcreteStrength');
    results(iData).rhos         = data(iData).section.getSectionData('SteelRatio');
    results(iData).rhosr        = data(iData).section.getSectionData('ReinforcingRatio');
    results(iData).ShapeFactor  = data(iData).section.getSectionData('ShapeFactor',data(iData).axis);

    results(iData).K            = BA_Elastic.K;
    results(iData).StoryBasedK  = BA_Elastic.storyBasedK;
    results(iData).L            = data(iData).L;
    results(iData).sidesway     = BA_Elastic.type2; % Returns "U" or "I"
    %results(iData).endMomentRatio(ind)    = BA_Elastic.endMomentRatio;

    results(iData).Pn_ana      = results_FN(iData).P1(1);
     
    results(iData).lambdaoe1 = (results(iData).L/pi)*sqrt(results(iData).Pno/results(iData).EIeff);
    results(iData).lambdaoe  = results(iData).K*results(iData).lambdaoe1;
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