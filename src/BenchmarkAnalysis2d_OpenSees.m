classdef BenchmarkAnalysis2d_OpenSees < OpenSeesAnalysis
    
    properties
        
        % Cross Section Properties
        section
        axis 
        sectionDefinition
        sectionTag = 1;
        
        % Frame Properties
        L
        frame_type
        kqtop  = [];
        kqbot  = [];
        gamma  = [];
        beta   = [];
        Delta0 = 0;
        delta0 = 0;
        
        % Analysis Options
        analysisTclFile = fullfile(pathOf.BenchmarkStudySource,'BenchmarkAnalysis2d.tcl');
        numElements                             = 6;
        numIntegrationPoints                    = 3;
        geomTransfType                          = 'Corotational';
        includeInitialGeometricImperfections    = true;
        eigenvalueLimitPointTolerance           = 1e-3;
        absoluteStrainLimit                     = 0.05;
        concreteCompressiveStrainLimit          = [];
        compressiveForceLimitType               = '';
        numStepsGravity                         = 30;
        numStepsLateral                         = 100;  
        limitPointWarning                       = false;
        runZeroAxialLoadAsSection               = true;
        extraResultsOutput                      = false;
        base_tolerance_force                    = 0.1;
        determine_M1_with_initial_imperfections = false;
        
    end
    
    methods
        function obj = BenchmarkAnalysis2d_OpenSees(data,sectionDefinition,analysisOptions)

            % Cross Section Properties
            obj.section = data.section;
            obj.axis = data.axis;
            obj.sectionDefinition = sectionDefinition;
            
            % Frame Properties
            obj.L = data.L;
            obj.frame_type = data.frame_type;
            
            switch obj.frame_type
                case 'Section'
                    assert(obj.L == 0);
                case 'Sidesway_Inhibited'
                    obj.beta   = data.beta;
                    obj.delta0 = data.delta0;
                case 'Sidesway_Uninhibited'
                    obj.kqtop  = data.kqtop;
                    obj.kqbot  = data.kqbot;
                    obj.gamma  = data.gamma;
                    obj.Delta0 = data.Delta0;
                    obj.delta0 = data.delta0;
                otherwise
                    error('Unknown frame type: %s',obj.frame_type)
            end
            
            % Analysis Options
            if nargin > 2
                names = fieldnames(analysisOptions);
                pnames = properties(obj);
                for i = 1:length(names)
                    if ismember(names{i},pnames)
                        obj.(names{i}) = analysisOptions.(names{i});
                    elseif strcmp(names{i},'store_extra_data')
                        % do nothing
                    else
                        warning('Unknown option: %s',names{i})
                    end
                end
            end
        end
        
        function results = runAnalysis(obj,analysisType,X1,X2,tryNumber)
            % Compression load is negative
            
            if nargin < 5
                tryNumber = 1;
            end
            switch analysisType
                case 'LimitPoint_Proportional'
                    x_over_P = X1;
                    checkForLimitPoint = true;
                case 'LimitPoint_NonProportional'
                    P = X1;
                    checkForLimitPoint = true;
                    if P == 0 && obj.runZeroAxialLoadAsSection
                        results = obj.runSectionAnalysis(analysisType,X1,X2,tryNumber);
                        return
                    end
                case 'TargetForce_Proportional'
                    P        = X1;
                    x_over_P = X2;
                    checkForLimitPoint = false;
                case 'TargetForce_NonProportional'
                    P = X1;
                    x = X2;
                    checkForLimitPoint = false;                  
                % case 'TargetDisplacement_Proportional'
                case 'TargetDrift_NonProportional'
                    P = X1;
                    d = X2;
                    checkForLimitPoint = false;                  
                otherwise
                    error('Unknown analysis type: %s',analysisType)
            end
            
            % Filenames
            inputFilename                     = obj.scratchFile('BenchmarkAnalysis2d_Driver.tcl');
            nodeDisplacementFilename          = obj.scratchFile('BenchmarkAnalysis2d_NodeDisp.out');
            elementForceFilename              = obj.scratchFile('BenchmarkAnalysis2d_EleForce.out');
            elementSectionDeformationFilename = obj.scratchFile('BenchmarkAnalysis2d_EleSectionDef.out');
            elementSectionStiffnessFilename   = obj.scratchFile('BenchmarkAnalysis2d_EleSectionStiff.out');
            eigenFilename                     = obj.scratchFile('BenchmarkAnalysis2d_EigVals.out');
                                   
            %% Write Input File
            fid = fopen(inputFilename,'w');
            fprintf(fid,'package require OpenSeesComposite \n');
            fprintf(fid,'namespace import OpenSeesComposite::* \n');
            fprintf(fid,'model basic -ndm 2 -ndf 3 \n');
            
            % Column Length
            fprintf(fid,'set Lc %g \n',obj.L);
            
            % Section Definition
            fprintf(fid,'set columnSectionTag %i\n',obj.sectionTag);
            if iscell(obj.sectionDefinition)
                for i = 1:length(obj.sectionDefinition)
                    fprintf(fid,'%s \n',obj.sectionDefinition{i});
                end
            else
                fprintf(fid,'%s \n',obj.sectionDefinition);
            end
            
            % Information Specific to the Type of Frame
            if strcmp(obj.frame_type,'Sidesway_Uninhibited')
                % Initial Imperfections
                if obj.includeInitialGeometricImperfections
                    fprintf(fid,'set Delta0 %g \n',obj.Delta0);
                    fprintf(fid,'set delta0 %g \n',obj.delta0);
                else
                    fprintf(fid,'set Delta0 0.0 \n');
                    fprintf(fid,'set delta0 0.0 \n');
                end
                
                % Leaning Column Load
                fprintf(fid,'set gamma %g \n',obj.gamma);
                if (obj.gamma ~= 0)
                    fprintf(fid,'set leaningE %g \n',1.0);
                    fprintf(fid,'set leaningA %g \n',1.0e10);
                    fprintf(fid,'set leaningI %g \n',1.0e10);
                end
                
                % Rotational Stiffness at the Top
                if (obj.kqtop == Inf)
                    fprintf(fid,'set rotStiffTop Fixed \n');
                elseif (obj.kqtop == 0)
                    fprintf(fid,'set rotStiffTop Free \n');
                else
                    fprintf(fid,'set rotStiffTop %g \n',obj.kqtop);
                end
                
                % Rotational Stiffness at the Bottom
                if (obj.kqbot == Inf)
                    fprintf(fid,'set rotStiffBot Fixed \n');
                elseif (obj.kqbot == 0)
                    fprintf(fid,'set rotStiffBot Free \n');
                else
                    fprintf(fid,'set rotStiffBot %g \n',obj.kqbot);
                end
                
                % Loading Steps
                switch analysisType
                    case 'LimitPoint_Proportional'                       
                        if obj.Delta0 == 0
                            direction = 1;
                        else
                            direction = sign(obj.Delta0);
                        end
                        
                        switch tryNumber
                            case 1
                                
                                fprintf(fid,'set coarseDispStepSize %g \n',obj.L/50000*direction);
                                fprintf(fid,'set fineDispStepSize %g \n',obj.L/500000*direction);
                                fprintf(fid,'set maxDisp %g \n',0.10*obj.L);
                                fprintf(fid,'set maxAnalysisTime %i \n',240);
                                fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                                fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);
                            case 2
                                fprintf(fid,'set coarseDispStepSize %g \n',obj.L/2000000*direction);
                                fprintf(fid,'set fineDispStepSize %g \n',obj.L/2000000*direction);
                                fprintf(fid,'set maxDisp %g \n',0.10*obj.L);
                                fprintf(fid,'set maxAnalysisTime %i \n',240);
                                fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                                fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);
                            otherwise
                                error('Undefined for tryNumber: %i',tryNumber)
                        end
                        
                    case 'LimitPoint_NonProportional'
                        switch tryNumber
                            case 1
                                fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                                fprintf(fid,'set coarseDispStepSize %g \n',obj.L/2000);
                                fprintf(fid,'set fineDispStepSize %g \n',obj.L/20000);
                                fprintf(fid,'set maxDisp %g \n',0.15*obj.L);
                                fprintf(fid,'set maxAnalysisTime %i \n',60);
                                fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                                fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);
                            case 2
                                fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                                fprintf(fid,'set coarseDispStepSize %g \n',obj.L/100000);
                                fprintf(fid,'set fineDispStepSize %g \n',obj.L/100000);
                                fprintf(fid,'set maxDisp %g \n',0.15*obj.L);
                                fprintf(fid,'set maxAnalysisTime %i \n',240);
                                fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                                fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);
                            case 3
                                fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                                fprintf(fid,'set coarseDispStepSize %g \n',obj.L/100000);
                                fprintf(fid,'set fineDispStepSize %g \n',obj.L/1000000);
                                fprintf(fid,'set maxDisp %g \n',0.15*obj.L);
                                fprintf(fid,'set maxAnalysisTime %i \n',600);
                                fprintf(fid,'set baseForceTolerance %g \n',10*obj.base_tolerance_force);
                                fprintf(fid,'set baseDisplacementTolerance %g \n',1e-7*obj.L);
                            otherwise
                                error('Undefined for tryNumber: %i',tryNumber)
                        end       
                        
                    case 'TargetForce_Proportional'
                        fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                        fprintf(fid,'set numStepsLateral %i\n',obj.numStepsLateral);
                        fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                        fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);                           

                    case 'TargetForce_NonProportional'
                        fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                        fprintf(fid,'set numStepsLateral %i\n',obj.numStepsLateral);
                        fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                        fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);
                        
                    case 'TargetDrift_NonProportional'
                        fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                        fprintf(fid,'set numStepsLateral %i\n',obj.numStepsLateral);
                        fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                        fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);  
                        
                    otherwise
                        error('Unknown analysis type: %s',analysisType)
                end                
                
            elseif strcmp(obj.frame_type,'Sidesway_Inhibited')
                
                % Initial Imperfections
                if obj.includeInitialGeometricImperfections
                    fprintf(fid,'set delta0 %g \n',obj.delta0);
                else
                    fprintf(fid,'set delta0 0.0 \n');
                end
                
                
                % Moment Ratio
                fprintf(fid,'set beta %g \n',obj.beta);
                
                % Loading Steps
                switch analysisType
                    case 'LimitPoint_Proportional'
                        
                        if obj.delta0 == 0
                            direction = 1;
                        else
                            direction = sign(obj.delta0);
                        end
                        
                        switch tryNumber
                            case 1                        
                                fprintf(fid,'set coarseDispStepSize %g \n',obj.L/200000*direction);
                                fprintf(fid,'set fineDispStepSize %g \n',obj.L/200000*direction);
                                fprintf(fid,'set maxDisp %g \n',0.05*obj.L);
                                fprintf(fid,'set maxAnalysisTime %i \n',240);
                                fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                                fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);
                                
                            case 2
                                fprintf(fid,'set coarseDispStepSize %g \n',obj.L/2000000*direction);
                                fprintf(fid,'set fineDispStepSize %g \n',obj.L/2000000*direction);
                                fprintf(fid,'set maxDisp %g \n',0.05*obj.L);
                                fprintf(fid,'set maxAnalysisTime %i \n',600);
                                fprintf(fid,'set baseForceTolerance %g \n',0.1*obj.base_tolerance_force);
                                fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);

                            case 3
                                fprintf(fid,'set coarseDispStepSize %g \n',obj.L/20000000*direction);
                                fprintf(fid,'set fineDispStepSize %g \n',obj.L/20000000*direction);
                                fprintf(fid,'set maxDisp %g \n',0.05*obj.L);
                                fprintf(fid,'set maxAnalysisTime %i \n',600);
                                fprintf(fid,'set baseForceTolerance %g \n',0.01*obj.base_tolerance_force);
                                fprintf(fid,'set baseDisplacementTolerance %g \n',1e-6*obj.L);
                                                                
                            otherwise
                                error('Undefined for tryNumber: %i',tryNumber)
                        end
                        
                    case 'LimitPoint_NonProportional'
                        switch tryNumber
                            case 1
                                fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                                fprintf(fid,'set coarseDispStepSize %g \n',obj.L/10000);
                                fprintf(fid,'set fineDispStepSize %g \n',obj.L/100000);
                                fprintf(fid,'set maxDisp %g \n',0.10*obj.L);
                                % fprintf(fid,'set coarseDispStepSize %g \n',-1/200000);
                                % fprintf(fid,'set fineDispStepSize %g \n',-1/2000000);
                                % fprintf(fid,'set maxDisp %g \n',0.10);
                                fprintf(fid,'set maxAnalysisTime %i \n',60);
                                fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                                fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);

                            case 2
                                fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                                fprintf(fid,'set coarseDispStepSize %g \n',obj.L/1000000);
                                fprintf(fid,'set fineDispStepSize %g \n',obj.L/1000000);
                                fprintf(fid,'set maxDisp %g \n',0.10*obj.L);
                                % fprintf(fid,'set coarseDispStepSize %g \n',-1/10000000);
                                % fprintf(fid,'set fineDispStepSize %g \n',-1/10000000);
                                % fprintf(fid,'set maxDisp %g \n',0.10);
                                fprintf(fid,'set maxAnalysisTime %i \n',600);
                                fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                                fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);

                            case 3
                                fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                                fprintf(fid,'set coarseDispStepSize %g \n',obj.L/10000000);
                                fprintf(fid,'set fineDispStepSize %g \n',obj.L/10000000);
                                fprintf(fid,'set maxDisp %g \n',0.10*obj.L);
                                fprintf(fid,'set maxAnalysisTime %i \n',6000);
                                fprintf(fid,'set baseForceTolerance %g \n',10*obj.base_tolerance_force);
                                fprintf(fid,'set baseDisplacementTolerance %g \n',10*1e-8*obj.L);                                
                                
                            otherwise
                                error('Undefined for tryNumber: %i',tryNumber)
                        end
                        
                    case 'TargetForce_Proportional'
                        fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                        fprintf(fid,'set numStepsLateral %i\n',obj.numStepsLateral);
                        fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                        fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);
                        
                    case 'TargetForce_NonProportional'
                        fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                        fprintf(fid,'set numStepsLateral %i\n',obj.numStepsLateral);
                        fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                        fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);
                        
                    case 'TargetDrift_NonProportional'
                        %error('TargetDrift_NonProportional not implemented for Sidesway_Inhibited')                        
                        fprintf(fid,'set numStepsGravity %i\n',obj.numStepsGravity);
                        fprintf(fid,'set numStepsLateral %i\n',obj.numStepsLateral);
                        fprintf(fid,'set baseForceTolerance %g \n',obj.base_tolerance_force);
                        fprintf(fid,'set baseDisplacementTolerance %g \n',1e-8*obj.L);  
                        
                    otherwise
                        error('Unknown analysis type: %s',analysisType)
                end                     
                
            else
                error('Unknown frame type: %s',obj.frame_type);
            end
            
            
            % Information Specific to the Type of Loading
            switch analysisType
                case 'LimitPoint_Proportional'
                    fprintf(fid,'set analysisType LimitPoint_Proportional\n');
                    if strcmp(obj.frame_type,'Sidesway_Uninhibited')
                        fprintf(fid,'set H %g\n',x_over_P);
                    elseif strcmp(obj.frame_type,'Sidesway_Inhibited')
                        fprintf(fid,'set M %g\n',x_over_P);
                    else
                        error('Unknown frame type: %s',obj.frame_type);
                    end    
                    
                case 'LimitPoint_NonProportional'
                    fprintf(fid,'set analysisType LimitPoint_NonProportional\n');
                    fprintf(fid,'set P %g\n',P);     
                    
                case 'TargetForce_Proportional'
                    fprintf(fid,'set analysisType TargetForce_Proportional\n');
                    fprintf(fid,'set P %g\n',P);
                    if strcmp(obj.frame_type,'Sidesway_Uninhibited')
                        fprintf(fid,'set H %g\n',x_over_P);
                    elseif strcmp(obj.frame_type,'Sidesway_Inhibited')
                        fprintf(fid,'set M %g\n',x_over_P);
                    else
                        error('Unknown frame type: %s',obj.frame_type);
                    end                    
                    
                case 'TargetForce_NonProportional'
                    fprintf(fid,'set analysisType TargetForce_NonProportional\n');
                    fprintf(fid,'set P %g\n',P);
                    if strcmp(obj.frame_type,'Sidesway_Uninhibited')
                        fprintf(fid,'set H %g\n',x);
                    elseif strcmp(obj.frame_type,'Sidesway_Inhibited')
                        fprintf(fid,'set M %g\n',x);
                    else
                        error('Unknown frame type: %s',obj.frame_type);
                    end  
                    
                case 'TargetDrift_NonProportional'
                    fprintf(fid,'set analysisType TargetDrift_NonProportional\n');
                    fprintf(fid,'set P %g\n',P);
                    fprintf(fid,'set d %g\n',d);
                    
                otherwise
                    error('Unknown analysis type: %s',analysisType)
            end
                       
            
            % Miscellaneous Analysis Parameters
            fprintf(fid,'set numEles %i\n',obj.numElements);
            fprintf(fid,'set midNode %i\n',0.5*obj.numElements+1);
            fprintf(fid,'set numNodes %i\n',obj.numElements+1);
            fprintf(fid,'set numIP %i\n',obj.numIntegrationPoints);
            fprintf(fid,'set geomTransfType %s\n',obj.geomTransfType);
            fprintf(fid,'set nodeDisplacementFilename {%s} \n',nodeDisplacementFilename);
            fprintf(fid,'set elementForceFilename {%s} \n',elementForceFilename);
            fprintf(fid,'set elementSectionDeformationFilename {%s} \n',elementSectionDeformationFilename);
            fprintf(fid,'set elementSectionStiffnessFilename {%s} \n',elementSectionStiffnessFilename);
            fprintf(fid,'set eigenFilename {%s} \n',eigenFilename);
            fprintf(fid,'set tryNumber %i \n',tryNumber);
            fprintf(fid,'set frameType %s \n',obj.frame_type);
            fprintf(fid,'source {%s} \n',obj.path_for_tcl(obj.analysisTclFile));
            
            % Close File
            fclose(fid);
            
            % Run Analysis
            [status, result] = obj.runOpenSees(inputFilename);

            % Store Analysis Data
            results = struct;
            results.member_type = 'Beam-Column';
            results.textOutput = result;
            switch status
                case 1
                    results.exitStatus = 'Peak Obtained';
                case 2
                    results.exitStatus = 'Step Limits Reached';
                case 3
                    fprintf('%s\n',result);
                    error('Analysis failed in gravity loading');
                case 4
                    results.exitStatus = 'Analysis Failed In Main Loading';
                case 5
                    results.exitStatus = 'Loads Reached';
                case 6
                    results.exitStatus = 'Analysis Failed In Lateral Loading';
                    %fprintf('%s\n',result);
                    %warning(results.exitStatus);
                case 7
                    results.exitStatus = 'Deformation Limits Reached';
                case 8
                    results.exitStatus = 'Analysis Timed Out';
                otherwise
                    fprintf('%s\n',result);
                    error('Analysis ended in an unknown manner, exit code = %i',status);
            end
            
            
            %% Read Results
            temp = dlmread(nodeDisplacementFilename);
            nodeDisp.x = temp(:,1:3:end);
            nodeDisp.y = temp(:,2:3:end);
            nodeDisp.q = temp(:,3:3:end);
            temp = dlmread(elementForceFilename);
            eleForce.Pi = temp(:,1:6:end);
            eleForce.Vi = temp(:,2:6:end);
            eleForce.Mi = temp(:,3:6:end);
            eleForce.Pj = temp(:,4:6:end);
            eleForce.Vj = temp(:,5:6:end);
            eleForce.Mj = temp(:,6:6:end);
            temp = dlmread(elementSectionDeformationFilename);
            sectionDeformation.axial     = temp(:,1:2:end);
            sectionDeformation.curvature = temp(:,2:2:end);
            temp = dlmread(elementSectionStiffnessFilename);
            sectionStiffness.EA = temp(:,1:2:end);
            sectionStiffness.EI = temp(:,2:2:end);
            temp = dlmread(eigenFilename);
            time = temp(:,1);
            lowestEigenvalue = temp(:,2);                
            
            
            %% Post-Process Results
                        
            % Applied Loads
            switch analysisType
                case {'LimitPoint_Proportional','TargetForce_Proportional'}
                    ind_grav = [];
                    ind_main = 1:numel(time);

                    P1 = -time;
                    x  = x_over_P*time;
                    
                case {'LimitPoint_NonProportional','TargetForce_NonProportional','TargetDrift_NonProportional'}
                    if P == 0
                        ind_grav = [];
                        ind_main = 1:numel(time);

                        P1 = zeros(size(time));
                        x  = time;
                    else
                        ind_grav = 1:(obj.numStepsGravity+1);
                        ind_main = (obj.numStepsGravity+1):numel(time);

                        time_1 = time(1:(obj.numStepsGravity+1));
                        P1_1 = time_1*P;
                        x_1  = zeros(size(P1_1));
                        
                        time_2 = time((obj.numStepsGravity+2):end);
                        P1_2 = ones(size(time_2))*P;
                        x_2  = time_2;
                        
                        P1 = vertcat(P1_1,P1_2);
                        x  = vertcat(x_1,x_2);
                    end                    
                    
                otherwise
                    error('Unknown analysis type: %s',analysisType)
            end
            
            % First Order Moment
            data = struct;
            data_fields = fields(obj);
            for i = 1:length(data_fields)
                data.(data_fields{i}) = obj.(data_fields{i});
            end
            data.EI = data.section.EI(obj.axis,'Gross');
            BA_Elastic = BenchmarkAnalysis2d_Elastic(data);
            BA_Elastic.includeInitialGeometricImperfections = ...
                obj.determine_M1_with_initial_imperfections;
            
            M1 = nan(size(P1));
            for i = 1:length(M1)
                if P1(i) > 0
                    error('Positive P')
                end
                M1(i) = BA_Elastic.maxFirstOrderMoment(P1(i),x(i));
            end
            
            % Internal Loads
            M2 = max(abs(horzcat(eleForce.Mi,eleForce.Mj)),[],2);
            P2 = min(horzcat(-eleForce.Pi,eleForce.Pj),[],2);
            
            % Deformations
            def = max(abs(nodeDisp.x),[],2);
            
            % Strains
            maxAbsoluteStrain = ...
                obj.section.longitudinalStrain2d(obj.axis,...
                sectionDeformation.axial,sectionDeformation.curvature,...
                'MaxAbsolute');
            maxAbsoluteStrain = max(maxAbsoluteStrain,[],2);
            
            if obj.section.hasConcrete
                maxConcreteCompressiveStrain = ...
                    obj.section.longitudinalStrain2d(obj.axis,...
                    sectionDeformation.axial,sectionDeformation.curvature,...
                    'MaxConcreteCompressive');
                maxConcreteCompressiveStrain = max(maxConcreteCompressiveStrain,[],2);
            end
            
            
            % Store Path Results
            results.path.eigen = lowestEigenvalue;
            results.path.P1    = P1;
            results.path.M1    = M1;
            results.path.P2    = P2;
            results.path.M2    = M2;
            results.path.def   = def;
            results.path.maxAbsoluteStrain = maxAbsoluteStrain;
            if obj.section.hasConcrete
                results.path.maxConcreteCompressiveStrain = maxConcreteCompressiveStrain;
            else
                results.path.maxConcreteCompressiveStrain = [];
            end
            results.path.topDisplacement = nodeDisp.x(:,end);
            results.path.botRotation = -nodeDisp.q(:,1);
            
            if ~isempty(ind_main)
                results.mainPath.eigen = lowestEigenvalue(ind_main);
                results.mainPath.P1    = P1(ind_main);
                results.mainPath.M1    = M1(ind_main);
                results.mainPath.P2    = P2(ind_main);
                results.mainPath.M2    = M2(ind_main);
                results.mainPath.def   = def(ind_main);
                results.mainPath.maxAbsoluteStrain = maxAbsoluteStrain(ind_main);
                if obj.section.hasConcrete
                    results.mainPath.maxConcreteCompressiveStrain = maxConcreteCompressiveStrain(ind_main);
                else
                    results.mainPath.maxConcreteCompressiveStrain = [];
                end
                results.mainPath.topDisplacement = nodeDisp.x(ind_main,end);
                results.mainPath.botRotation = -nodeDisp.q(ind_main,1);
            end
            
            if ~isempty(ind_grav)
                results.gravPath.eigen = lowestEigenvalue(ind_grav);
                results.gravPath.P1    = P1(ind_grav);
                results.gravPath.M1    = M1(ind_grav);
                results.gravPath.P2    = P2(ind_grav);
                results.gravPath.M2    = M2(ind_grav);
                results.gravPath.def   = def(ind_grav);
                results.gravPath.maxAbsoluteStrain = maxAbsoluteStrain(ind_grav);
                if obj.section.hasConcrete
                    results.gravPath.maxConcreteCompressiveStrain = maxConcreteCompressiveStrain(ind_grav);
                else
                    results.gravPath.maxConcreteCompressiveStrain = [];
                end
                results.gravPath.topDisplacement = nodeDisp.x(ind_grav,end);
                results.gravPath.botRotation = -nodeDisp.q(ind_grav,1);
            end
            
            
            % Limit Point
            if checkForLimitPoint
                results = obj.findLimitPoint(results);
            end
            
            % Extra Results
            if obj.extraResultsOutput
                nEL = obj.numElements;
                nIP = obj.numIntegrationPoints;
                Lele = obj.L/nEL;
                
                % Distribution along length by elements
                xi = nan(1,2*nEL);
                temp = linspace(0,obj.L,nEL+1);
                xi(1:2:2*nEL) = temp(1:(end-1));
                xi(2:2:2*nEL) = temp(2:end);
                P = nan(length(time),2*nEL);
                M = nan(length(time),2*nEL);
                V = nan(length(time),2*nEL);
                P(:,1:2:2*nEL) =  eleForce.Pi;
                P(:,2:2:2*nEL) = -eleForce.Pj;
                M(:,1:2:2*nEL) =  eleForce.Mi;
                M(:,2:2:2*nEL) = -eleForce.Mj;
                V(:,1:2:2*nEL) =  eleForce.Vi;
                V(:,2:2:2*nEL) = -eleForce.Vj;
                
                % Store results
                results.distE.xi = xi;
                results.distE.P = P;
                results.distE.M = M;
                results.distE.V = V;
                               
                
                % Distribution along length by sections
                [ixi,iwt] = lobatto_quadrature(nIP,Lele);
                xi = nan(1,nEL*nIP);
                wt = nan(1,nEL*nIP);
                for i = 1:nEL
                    xi(1+(i-1)*nIP:i*nIP) = ixi + (i-1)*Lele;
                    wt(1+(i-1)*nIP:i*nIP) = iwt;
                end
           
                % Store Results
                results.distS.xi = xi;
                results.distS.wt = wt;
                results.distS.axialStrain = sectionDeformation.axial;
                results.distS.curvature = sectionDeformation.curvature;
                results.distS.EA = sectionStiffness.EA;
                results.distS.EI = sectionStiffness.EI;
            end
            
            % Clean Folder
            if obj.deleteFilesAfterAnalysis
                delete(inputFilename,nodeDisplacementFilename,...
                    elementForceFilename,elementSectionDeformationFilename,...
                    elementSectionStiffnessFilename,eigenFilename)
            end
            
        end
        
        function results = runSectionAnalysis(obj,analysisType,X1,X2,tryNumber)
            if nargin < 5
                tryNumber = 1;
            end
            switch analysisType
                case 'LimitPoint_Axial'
                    checkForLimitPoint = true;
                case 'LimitPoint_Proportional'
                    e = X1;
                    checkForLimitPoint = true;
                case 'LimitPoint_Proportional2'
                    e = X1;
                    checkForLimitPoint = true;                    
                case 'LimitPoint_NonProportional'
                    P = X1;
                    checkForLimitPoint = true;
                case 'TargetForce_Proportional'
                    P = X1;
                    M = X2;
                    checkForLimitPoint = false;
                case 'TargetForce_NonProportional'
                    P = X1;
                    M = X2;
                    checkForLimitPoint = false;                  
                % case 'TargetDisplacement_Proportional'
                % case 'TargetDisplacement_NonProportional'
                case 'TargetCurvature_NonProportional'
                    P = X1;
                    kappa = X2;
                    checkForLimitPoint = false;
                otherwise
                    error('Unknown analysis type: %s',analysisType)
            end
            
            % Run Analysis
            si = SectionAnalysis(obj.sectionDefinition);
            si.scratchPath              = obj.scratchPath;
            si.deleteFilesAfterAnalysis = obj.deleteFilesAfterAnalysis;
            
            % Step Sizes
            switch tryNumber
                case 1
                    deformationStep_Axial = -0.00005;
                    maxNumSteps_Axial = 1000;
                    deformationStep_Bending = 0.0001/obj.section.depth(obj.axis);
                    maxNumSteps_Bending = 3000;
                case 2
                    deformationStep_Axial = -0.00005;
                    maxNumSteps_Axial = 1000;
                    deformationStep_Bending = 0.00001/obj.section.depth(obj.axis);
                    maxNumSteps_Bending = 30000;
                case 3
                    deformationStep_Axial = -0.00005;
                    maxNumSteps_Axial = 1000;
                    deformationStep_Bending = 0.000001/obj.section.depth(obj.axis);
                    maxNumSteps_Bending = 300000;
                otherwise
                    error('Undefined for tryNumber: %i',tryNumber)
            end
            
            switch analysisType
                case 'LimitPoint_Axial'
                    siResults = si.runSectionToPeak2d('AxialOnly',...
                        deformationStep_Axial,maxNumSteps_Axial);
                case 'LimitPoint_Proportional'
                    siResults = si.runSectionToPeak2d('Proportional',...
                        deformationStep_Axial,maxNumSteps_Axial,e);
                case 'LimitPoint_Proportional2'
                    siResults = si.runSectionToPeak2d('Proportional2',...
                        deformationStep_Bending,maxNumSteps_Bending,e);
                case 'LimitPoint_NonProportional'
                    siResults = si.runSectionToPeak2d('NonProportional',...
                        deformationStep_Bending,maxNumSteps_Bending,P,10);
                case 'TargetForce_Proportional';
                    siResults = si.runSectionAnalysis([2 0 2],[0 0; P M],[P/10 M/200]);
                case 'TargetForce_NonProportional';
                    siResults = si.runSectionAnalysis([2 0 2],[0 0; P 0; P M],[P/10 M/200]);
                case 'TargetCurvature_NonProportional'
                    siResults = si.runSectionAnalysis([2 0 1],[0 0; P 0; P kappa],[P/10 kappa/200]);
                otherwise
                    error('Unknown loadingType');
            end
                    
            
            % Store Analysis Data
            results = struct;
            results.member_type = 'Section';
            results.textOutput  = siResults.textOutput;
            results.exitStatus  = siResults.status;
            
            % Store Path Results
            results.path.P1             = siResults.axialForce;
            results.path.M1             = siResults.moment;
            results.path.P2             = results.path.P1;
            results.path.M2             = results.path.M1;
            results.path.def            = nan(size(results.path.P1));
            results.path.axialStrain    = siResults.axialStrain;
            results.path.curvature      = siResults.curvature;
            
            
            results.path.maxAbsoluteStrain = obj.section.longitudinalStrain2d(...
                obj.axis,...
                results.path.axialStrain,results.path.curvature,...
                'MaxAbsolute');
            
            if obj.section.hasConcrete
                results.path.maxConcreteCompressiveStrain = obj.section.longitudinalStrain2d(...
                    obj.axis,...
                    results.path.axialStrain,results.path.curvature,...
                    'MaxConcreteCompressive');
            end
            
            
            % Limit Point
            if checkForLimitPoint
                results.path.eigen = siResults.eigen;
                results = obj.findLimitPoint(results);
            end 
            
        end
        
        function results = findLimitPoint(obj,results)
            
            % Find Stability Limit Point
            [ind,x] = find_limit_point_in_vector(results.path.eigen,0.0);
            
            if ~isempty(ind)
                stabilityLimitPoint.limit_type = 'Stability Limit (Reached)';
                stabilityLimitPoint.good = true;
                stabilityLimitPoint.time = ind+x;
            else
                [ind,x] = find_limit_point_in_vector(...
                    results.path.eigen/results.path.eigen(1),...
                    obj.eigenvalueLimitPointTolerance);
                
                if isempty(ind)
                    stabilityLimitPoint.limit_type = 'Stability Limit (Not Reached)';
                    stabilityLimitPoint.good = false;
                    stabilityLimitPoint.time = Inf;
                else
                    stabilityLimitPoint.limit_type = 'Stability Limit (Tolerance Reached)';
                    stabilityLimitPoint.good = true;
                    stabilityLimitPoint.time = ind+x;
                end
            end
            
            % Find Strain Limit
            [ind,x] = find_limit_point_in_vector(...
                results.path.maxAbsoluteStrain,obj.absoluteStrainLimit);
            
            if isempty(ind)
                strainLimitPoint.limit_type = 'Strain Limit (Not Reached)';
                strainLimitPoint.good = false;
                strainLimitPoint.time = Inf;
            else
                strainLimitPoint.limit_type = 'Strain Limit (Reached)';
                strainLimitPoint.good = true;
                strainLimitPoint.time = ind+x;
            end
            
            % Find Concrete Compressive Strain Limit Point
            if ~isempty(obj.concreteCompressiveStrainLimit) && obj.section.hasConcrete
                [ind,x] = find_limit_point_in_vector(...
                    results.path.maxConcreteCompressiveStrain,...
                    obj.concreteCompressiveStrainLimit);
                
                if isempty(ind)
                    concreteCompressiveStrainLimitPoint.limit_type = 'Concrete Compressive Strain Limit (Not Reached)';
                    concreteCompressiveStrainLimitPoint.good = false;
                    concreteCompressiveStrainLimitPoint.time = Inf;
                else
                    concreteCompressiveStrainLimitPoint.limit_type = 'Concrete Compressive Strain Limit (Reached)';
                    concreteCompressiveStrainLimitPoint.good = true;
                    concreteCompressiveStrainLimitPoint.time = ind+x;
                end
            else
                concreteCompressiveStrainLimitPoint.limit_type = 'Concrete Compressive Strain Limit (Not Checked)';
                concreteCompressiveStrainLimitPoint.good = false;
                concreteCompressiveStrainLimitPoint.time = Inf;
            end
            
            % Find Compressive Force Limit Point
            if ~isempty(obj.compressiveForceLimitType) && ~strcmpi(obj.compressiveForceLimitType,'none')
                
                switch lower(obj.compressiveForceLimitType)
                    case 'pnco'
                        compressiveForceLimit = obj.section.Pnco;
                    otherwise
                        error('Unknown compressiveForceLimitType')
                end
                
                [ind,x] = find_limit_point_in_vector(results.path.P2,compressiveForceLimit);
                
                if isempty(ind)
                    compressiveForceLimitPoint.limit_type = 'Compressive Force Limit (Not Reached)';
                    compressiveForceLimitPoint.good = false;
                    compressiveForceLimitPoint.time = Inf;
                else
                    compressiveForceLimitPoint.limit_type = 'Compressive Force Limit (Reached)';
                    compressiveForceLimitPoint.good = true;
                    compressiveForceLimitPoint.time = ind+x;
                end
            else
                compressiveForceLimitPoint.limit_type = 'Compressive Force Limit (Not Checked)';
                compressiveForceLimitPoint.good = false;
                compressiveForceLimitPoint.time = Inf;
            end
            
            % Determine Controling Limit
            limitPointTimes = [
                stabilityLimitPoint.time
                strainLimitPoint.time
                concreteCompressiveStrainLimitPoint.time
                compressiveForceLimitPoint.time ];
            controllingLimitPointTime = min(limitPointTimes);
            
            if stabilityLimitPoint.time == controllingLimitPointTime
                results.limitPoint = stabilityLimitPoint;
            elseif strainLimitPoint.time == controllingLimitPointTime
                results.limitPoint = strainLimitPoint;
            elseif concreteCompressiveStrainLimitPoint.time == controllingLimitPointTime
                results.limitPoint = concreteCompressiveStrainLimitPoint;
            elseif compressiveForceLimitPoint.time == controllingLimitPointTime
                results.limitPoint = compressiveForceLimitPoint;
            else
                error('error')
            end
            
            % Warn if limit point could not be determined
            if obj.limitPointWarning && ~results.limitPoint.good
                warning('Limit point could not be determined, \n   Exit Status: %s \n   Eigenvalue Ratio: %g (limit: %g)\n   Max Strain: %g (limit: %g)',...
                    results.exitStatus,...
                    min(results.path.eigen/results.path.eigen(1)),obj.eigenvalueLimitPointTolerance,...
                    max(results.path.maxAbsoluteStrain),obj.absoluteStrainLimit);
            end
            
            % Compute and Store Limit Point Data
            if results.limitPoint.time ~= Inf
                ind = floor(results.limitPoint.time);
                x = results.limitPoint.time - ind;
                results.limitPoint.P1  = interpolate_vector(results.path.P1,ind,x);
                results.limitPoint.M1  = interpolate_vector(results.path.M1,ind,x);
                results.limitPoint.P2  = interpolate_vector(results.path.P2,ind,x);
                results.limitPoint.M2  = interpolate_vector(results.path.M2,ind,x);
                results.limitPoint.def = interpolate_vector(results.path.def,ind,x);
            else
                results.limitPoint.P1  = 0.0;
                results.limitPoint.M1  = 0.0;
                results.limitPoint.P2  = 0.0;
                results.limitPoint.M2  = 0.0;
                results.limitPoint.def = 0.0;
            end
            
        end
        
    end
    
end

