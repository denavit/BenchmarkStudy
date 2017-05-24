classdef BenchmarkAnalysis2d_Elastic_Base
        
    properties (Abstract = true, SetAccess = immutable)
        EI      % Elastic flexural stiffness
        L       % Length of the member
        K       % Effective Length Factor
    end
    
    methods (Abstract,Static)
        t = type()
        t = type2()
    end
  
    methods (Abstract)
        x = endMomentRatio(obj)
        x = driftRatio(obj,P)
        v = firstOrderDisplacement(obj,H,x)
        v = secondOrderDisplacement(obj,P,H,x)
        vmax = maxFirstOrderDisplacement(obj,H)
        vmax = maxSecondOrderDisplacement(obj,P,H)
        m = firstOrderMoment(obj,H,x)
        m = secondOrderMoment(obj,P,H,x)
        mmax = maxFirstOrderMoment(obj,H)
        mmax = maxSecondOrderMoment(obj,P,H)
        m2m1 = secondOrderMomentRatio(obj,P)
        M1 = determineAppliedMoment(obj,P,M2)
        [Pmax,M2atPmax] = maximumLoad(obj,section,...
            designStrengthType,...
            notionalLoadObject,...
            effectiveLengthFactorType,...
            columnStiffnessReduction,beamStiffnessReduction,tauType)
        [P1,M1,P2,M2] = designInteraction(obj,section,...
            designStrengthType,numPoints,...
            notionalLoadObject,...
            effectiveLengthFactorType,...
            columnStiffnessReduction,beamStiffnessReduction,tauType)
        [P2e,M2e] = impliedInteraction(obj,P1_FN,M1_FN,...
            notionalLoadObject,...
            columnStiffnessReduction,beamStiffnessReduction,tauType,Py)
        [drift_1,drift_2,driftRatio] = computeDrifts(obj,P1,M1)
        p = determineLoadForMomentRatio(obj,momentRatio,...
            columnStiffnessReduction,beamStiffnessReduction,tauType,Py_tau);
    end
    
    methods
        function x = position(obj)
            numPoints = 100;
            x = linspace(0,obj.L,numPoints);
        end
        function p = eulerLoad(obj)
            p = pi^2*obj.EI/(obj.K*obj.L)^2;
        end        
        function x = displacementRatio(obj,P,x)
            x = obj.secondOrderDisplacement(P,1,x)/obj.firstOrderDisplacement(1,x);
        end        
        function x = determineAppliedLoadFromM1(obj,P,M1)
            x = M1/obj.maxFirstOrderMoment(P,1);
        end
        function EIeff = determineElasticStiffness(obj,P1,M1,vmax,type)
            
            if nargin < 5
                type = 'd';
            end
            
            assert(isequal(size(P1),size(M1),size(vmax)),...
                'P1, M1, and vmax should be arrays of the same size');
            
            EIeff = zeros(size(vmax));
            
            for i = 1:numel(EIeff)
                if M1(i) == 0 || vmax(i) == 0
                    EIeff(i) = NaN;
                else
                    x1 = obj.determineAppliedLoadFromM1(P1(i),M1(i));
                    
                    % Use fsolve to find effective stiffness
                    options.TolFun = 1.0e-12*vmax(i);
                    %options.TolX = abs(1.0e-9*obj.EI);
                    options.Display = 'off';
                    
                    switch type
                        case 'd'
                            [EIeff(i),~,exitflag] =  fsolve(...
                                @(x)obj.errorEIeff_d(P1(i),x1,vmax(i),x),...
                                obj.EI,options);
                            if exitflag <= 0
                                [EIeff(i),~,exitflag] =  fsolve(...
                                    @(x)obj.errorEIeff_d(P1(i),x1,vmax(i),x),...
                                    0.5*obj.EI,options);
                                if exitflag <= 0
                                    [EIeff(i),~,exitflag] =  fsolve(...
                                        @(x)obj.errorEIeff_d(P1(i),x1,vmax(i),x),...
                                        1.1*obj.EI,options);
                                    if exitflag <= 0
                                        warning('fsolve could not find solution');
                                        EIeff(i) = NaN;
                                    end
                                end
                            end
                        case 'f'
                            [EIeff(i),~,exitflag] =  fsolve(...
                                @(x)obj.errorEIeff_f(P1(i),x1,vmax(i),x),...
                                obj.EI,options);
                            if exitflag <= 0
                                [EIeff(i),~,exitflag] =  fsolve(...
                                    @(x)obj.errorEIeff_f(P1(i),x1,vmax(i),x),...
                                    0.5*obj.EI,options);
                                if exitflag <= 0
                                    [EIeff(i),~,exitflag] =  fsolve(...
                                        @(x)obj.errorEIeff_f(P1(i),x1,vmax(i),x),...
                                        1.1*obj.EI,options);
                                    if exitflag <= 0
                                        warning('fsolve could not find solution');
                                        EIeff(i) = NaN;
                                    end
                                end
                            end                            
                        otherwise
                            error('')
                    end
                end
            end
        end
    end 
    
    methods (Static)
        function obj = BenchmarkAnalysis2d_Elastic(data)
            % Define Frame
            switch data.frame_type  
                case 'Sidesway_Uninhibited'
                    obj = BenchmarkAnalysis2d_Elastic_Sidesway_Uninhibited(data);
                case 'Sidesway_Inhibited'
                    obj = BenchmarkAnalysis2d_Elastic_Sidesway_Inhibited(data);
                otherwise
                    error('Unknown frame type: %s',data.frame_type);
            end            
        end
        function elasticFrame = getFrame(section,frame,elasticStiffness,Lcol)
            EIgross  = section.section.EIgross(section.axis);
            Pnogross = section.section.getSectionData('GrossSectionCompressionStrength');                      
            
            % Column Length
            if nargin < 4
                Lcol = frame.lambdaoe1g*pi*sqrt(EIgross/Pnogross);
            end
            
            % Elastic Flexural Stiffness
            if isnumeric(elasticStiffness) && isscalar(elasticStiffness)
                EIelastic = elasticStiffness;
            elseif ischar(elasticStiffness)
                [E,~,I] = section.section.sectionPropertiesForElasticAnalysis2d(...
                    section.axis,elasticStiffness);
                EIelastic = E*I;
            else
                error('Bad elasticStiffness')
            end
            
            % Define Frame
            switch frame.type  
                case 'sideswayUninhibited'
                    kqtop = (6*EIgross)/(frame.Ggtop*Lcol);
                    kqbot = (6*EIgross)/(frame.Ggbot*Lcol);
                    elasticFrame = BenchmarkAnalysis2d_Elastic_Sidesway_Uninhibited(...
                        EIelastic,Lcol,kqtop,kqbot,frame.gamma);
                case 'sideswayInhibited'
                    elasticFrame = BenchmarkAnalysis2d_Elastic_Sidesway_Inhibited(...
                        EIelastic,Lcol,frame.beta);
                otherwise
                    error('Unknown frame type');
            end
        end
        function l = lambdaoe(section,frame,EIeff)
            Pno = section.section.Pnco; 
            elasticFrame = BenchmarkAnalysis2d_Elastic.getFrame(section,frame,'Gross');
            if isnumeric(EIeff) && isscalar(EIeff)
                EI = EIeff;
            elseif ischar(EIeff)
                [E,~,I] = section.section.sectionPropertiesForElasticAnalysis2d(...
                    section.axis,EIeff);
                EI = E*I;
            else
                error('Bad EIeff')
            end
            l = (elasticFrame.K*elasticFrame.L/pi)*sqrt(Pno/EI);
        end        
    end
end
 
