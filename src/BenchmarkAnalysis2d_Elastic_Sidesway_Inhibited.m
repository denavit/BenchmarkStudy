classdef BenchmarkAnalysis2d_Elastic_Sidesway_Inhibited < BenchmarkAnalysis2d_Elastic_Base

    properties (SetAccess = immutable)
        EI      % Elastic flexural stiffness
        L       % Length of the member
        K       % Effective Length Factor        
        beta    % End moment ratio
        delta0  % Initial out-of-straightness
    end
    
    properties
        includeInitialGeometricImperfections = true;
    end

    methods
        function obj = BenchmarkAnalysis2d_Elastic_Sidesway_Inhibited(EI,L,beta,delta0)
            if nargin == 1
                data = EI;
                obj.EI      = data.EI;
                obj.L       = data.L;
                obj.beta    = data.beta;
                obj.delta0  = data.delta0;
            else
                obj.EI      = EI;
                obj.L       = L;
                obj.beta    = beta;
                obj.delta0  = delta0;
            end
            
            if obj.beta == -1
                obj.K   = 0.5;
            else
                obj.K   = 1;            
            end
        end
        function x = endMomentRatio(obj)
            if abs(obj.beta) > 1
                x = -1/obj.beta;
            else
                x = -obj.beta;
            end
        end
        function x = driftRatio(obj,P)
            x = 1;
        end  
        function k = storyBasedK(obj)
            if obj.beta == -1
                k = 0.5;
            else
                k = 1;            
            end
        end
        function v = firstOrderDisplacement(obj,P,M,x)
            if nargin < 4
                x = obj.position;
            end
            
            % Check for zero axial load case
            if P > 0
                error('Function only valid for compressive loads');
            else
                % Compressive loads should be positive
                P = -P;
            end                
            
            alpha = sqrt(P/obj.EI);
            if obj.includeInitialGeometricImperfections
                v = M*(obj.L-x).*x.*(x*(-1+obj.beta)+obj.L*(obj.beta+2))/...
                    (6*obj.EI*obj.L) + alpha^2*obj.delta0*obj.L^2*sin(pi*x/obj.L)/pi^2;
            else
                v = M*(obj.L-x).*x.*(x*(-1+obj.beta)+obj.L*(obj.beta+2))/...
                    (6*obj.EI*obj.L);
            end
        end
        function v = secondOrderDisplacement(obj,P,M,x)
            if nargin < 4
                x = obj.position;
            end
            
            % Check for zero axial load case
            if P == 0
                v = obj.firstOrderDisplacement(P,M,x);
                return
            elseif P > 0
                error('Function only valid for compressive loads');
            else
                % Compressive loads should be positive
                P = -P;
            end

            alpha = sqrt(P/obj.EI);
            if obj.includeInitialGeometricImperfections
                v = M*((x-obj.L) - x*obj.beta + obj.L*cos(x*alpha) + ...
                    obj.L*csc(obj.L*alpha)*(obj.beta-cos(obj.L*alpha)).*sin(x*alpha))/...
                    (obj.EI*obj.L*alpha^2) + ...
                    (alpha^2*obj.delta0*obj.L^2*sin(pi*x/obj.L))/(pi^2-alpha^2*obj.L^2);
            else
                v = M*((x-obj.L) - x*obj.beta + obj.L*cos(x*alpha) + ...
                    obj.L*csc(obj.L*alpha)*(obj.beta-cos(obj.L*alpha)).*sin(x*alpha))/...
                    (obj.EI*obj.L*alpha^2);
            end
        end
        function vmax = maxFirstOrderDisplacement(obj,P,M)
            vmax = max(abs(obj.firstOrderDisplacement(P,M)));
        end
        function vmax = maxSecondOrderDisplacement(obj,P,M)
            vmax = max(abs(obj.secondOrderDisplacement(P,M)));
        end
        function m = firstOrderMoment(obj,P,M,x)
            if nargin < 4
                x = obj.position;
            end
            
            % Check for zero axial load case
            if P > 0
                error('Function only valid for compressive loads');
            else
                % Compressive loads should be positive for this equation
                P = -P;
            end            
            
            alpha = sqrt(P/obj.EI);
            if obj.includeInitialGeometricImperfections
                m = M*(1+(obj.beta-1)*x/obj.L) + alpha^2*obj.delta0*obj.L^2*sin(pi*x/obj.L)/pi^2;
            else
                m = M*(1+(obj.beta-1)*x/obj.L);
            end
        end
        function m = secondOrderMoment(obj,P,M,x)
            if nargin < 4
                x = obj.position;
            end
            
            % Check for zero axial load case
            if P == 0
                m = obj.firstOrderMoment(P,M);
                return
            elseif P > 0
                error('Function only valid for compressive loads');
            else
                % Compressive loads should be positive
                P = -P;
            end

            alpha = sqrt(P/obj.EI);
            if obj.includeInitialGeometricImperfections
                m = M*(cos(x*alpha)+(-cot(obj.L*alpha)+obj.beta*csc(obj.L*alpha))*sin(x*alpha)) + ...
                    (alpha^2*obj.delta0*obj.EI*pi^2*sin(pi*x/obj.L))/(pi^2-alpha^2*obj.L^2);
            else
                m = M*(cos(x*alpha)+(-cot(obj.L*alpha)+obj.beta*csc(obj.L*alpha))*sin(x*alpha));
            end
        end
        function mmax = maxFirstOrderMoment(obj,P,M)
            m = obj.firstOrderMoment(P,M);
            mmax = max(abs(m));
        end
        function mmax = maxSecondOrderMoment(obj,P,M)
            m = obj.secondOrderMoment(P,M);
            mmax = max(abs(m));
        end
        function m2m1 = secondOrderMomentRatio(obj,P)
            m2m1 = max(abs(obj.secondOrderMoment(P,1)))/max(abs(obj.firstOrderMoment(P,1)));
            % @todo - this will be weird
        end   
        function M1 = determineAppliedMoment(obj,P,M2)
            %options.TolFun = 1.0e-6*M2;
            %options.TolX = 1.0e-6*M2;
            options = struct;
            options.Display = 'off';
            [M1,~,exitflag] =  fsolve(...
                @(M)obj.maxSecondOrderMoment(P,M)-M2,...
                0,options);
            if exitflag <= 0
                error('fsolve could not find solution');
            end
        end
        function [P,M2] = determinePeakLoad(obj,designM,designP,M1,Py,tauType)
            if M1 == 0
                % No lateral loads
                id = interactionDiagram2d(designM,designP);
                Pcr_design = id.findYgivenX(0,'negative');
                assert(~isempty(Pcr_design),'Could not find Pcr_design')
                
                Pguess = -0.8*min(Py,obj.eulerLoad);
                options.TolFun = 1.0e-12;
                options.TolX = abs(1.0e-8*Pguess);
                options.Display = 'off';
                [Pcr_analysis,~,exitflag,out] =  fsolve(...
                    @(P)errorCriticalLoadWithTau(obj,P,Py,tauType),...
                    Pguess,options);
                if exitflag <= 0
                    error('fsolve could not find solution');
                end
                
                if Pcr_design > Pcr_analysis
                    P  = Pcr_design;
                    M2 = 0;
                else
                    P  = Pcr_analysis;
                    id = interactionDiagram2d(designM,designP);
                    M2 = id.findXgivenY(P,'positive');
                    assert(~isempty(M2),'Could not find M2');
                end
                
            else 
                Pguess = -0.8*min(Py,obj.eulerLoad);
                options.TolFun = 1.0e-12;
                options.TolX = abs(1.0e-8*Pguess);
                options.Display = 'off';
                [P,~,exitflag,out] =  fsolve(...
                    @(P)errorPeakLoad(obj,designM,designP,P,M1,Py,tauType),...
                    Pguess,options);
                if exitflag <= 0
                    error('fsolve could not find solution');
                end
                tau = AISC_tau(-P/Py,tauType);
                obj2 = obj.get_copy_with_reduced_stiffness(tau);
                M2 = obj2.maxSecondOrderMoment(P,M1);
            end
        end
        function [P,M2] = determinePeakLoadWithEccentricity(obj,designM,designP,e,Py,tauType)
            noImperfMoment = obj.delta0 == 0 || ~obj.includeInitialGeometricImperfections;
            if e == 0 && noImperfMoment
                % No eccentricity
                id = interactionDiagram2d(designM,designP);
                Pcr_design = id.findYgivenX(0,'negative');
                assert(~isempty(Pcr_design),'Could not find Pcr_design')
                
                Pguess = -0.49*min(Py,obj.eulerLoad);
                options.TolFun = 1.0e-12;
                options.TolX = abs(1.0e-8*Pguess);
                options.Display = 'off';
                [Pcr_analysis,~,exitflag,out] =  fsolve(...
                    @(P)errorCriticalLoadWithTau(obj,P,Py,tauType),...
                    Pguess,options);
                if exitflag <= 0
                    error('fsolve could not find solution');
                end
                
                if Pcr_design > Pcr_analysis
                    P  = Pcr_design;
                    M2 = 0;
                else
                    P  = Pcr_analysis;
                    id = interactionDiagram2d(designM,designP);
                    M2 = id.findXgivenY(P,'positive');
                    assert(~isempty(M2),'Could not find M2');
                end
            else 
                Pguess = -0.49*min(Py,obj.eulerLoad);
                options.TolFun = 1.0e-12;
                options.TolX = abs(1.0e-8*Pguess);
                options.Display = 'off';
                [P,~,exitflag,out] =  fsolve(...
                    @(P)errorPeakLoad(obj,designM,designP,P,P*e,Py,tauType),...
                    Pguess,options);
                if exitflag <= 0
                    error('fsolve could not find solution');
                end
                tau = AISC_tau(-P/Py,tauType);
                obj2 = obj.get_copy_with_reduced_stiffness(tau);
                M2 = obj2.maxSecondOrderMoment(P,P*e);
            end
        end        
        function [Pmax,M2atPmax] = maximumLoad(obj,section,axis,...
                designStrengthType,...
                notionalLoadObject,...
                effectiveLengthFactorType,...
                columnStiffnessReduction,beamStiffnessReduction,tauType)

            % Determine Design Beam Column Interaction
            sectionKL = section;
            sectionKL.Lx = obj.L;
            sectionKL.Ly = obj.L;
            switch lower(effectiveLengthFactorType)
                case 'k'
                    sectionKL.Kx = obj.K;
                    sectionKL.Ky = obj.K;
                case 'storybasedk'
                    sectionKL.Kx = obj.storyBasedK;
                    sectionKL.Ky = obj.storyBasedK;
                case 'one'
                    sectionKL.Kx = 1.0;
                    sectionKL.Ky = 1.0;
                case 'zero'
                    sectionKL.Kx = 0.0;
                    sectionKL.Ky = 0.0;
                otherwise
                    error('Unknown effectiveLengthFactorType: %s',effectiveLengthFactorType);
            end
            [designP,designM] = ...
                sectionKL.beamColumnInteraction2d(axis,['Factored' designStrengthType],'CompPos');

            Py = section.Pnco;
            obj2 = obj.get_copy_with_reduced_stiffness(columnStiffnessReduction);
                
            % Run axial only analysis to get Pn
            [Pmax,M2atPmax] = obj2.determinePeakLoad(designM,designP,0.0,Py,tauType);        
        end
        function [P1,M1,P2,M2] = designInteraction(obj,section,axis,...
                designStrengthType,numPoints,...
                notionalLoadObject,...
                effectiveLengthFactorType,...
                columnStiffnessReduction,beamStiffnessReduction,tauType)

            % Determine Design Beam Column Interaction
            sectionKL = section;
            sectionKL.Lx = obj.L;
            sectionKL.Ly = obj.L;
            switch lower(effectiveLengthFactorType)
                case 'k'
                    sectionKL.Kx = obj.K;
                    sectionKL.Ky = obj.K;
                case 'storybasedk'
                    sectionKL.Kx = obj.storyBasedK;
                    sectionKL.Ky = obj.storyBasedK;
                case 'one'
                    sectionKL.Kx = 1.0;
                    sectionKL.Ky = 1.0;
                case 'zero'
                    sectionKL.Kx = 0.0;
                    sectionKL.Ky = 0.0;
                otherwise
                    error('Unknown effectiveLengthFactorType: %s',effectiveLengthFactorType);
            end
            [designP,designM] = ...
                sectionKL.beamColumnInteraction2d(axis,designStrengthType,'CompPos');

            % Initilize results
            P1 = zeros(numPoints,1);
            M1 = zeros(numPoints,1);
            P2 = zeros(numPoints,1);
            M2 = zeros(numPoints,1);

            Py = section.Pnco;
            obj2 = obj.get_copy_with_reduced_stiffness(columnStiffnessReduction);
            
            % Run axial only analysis to get Pn
            [p,m2] = obj2.determinePeakLoad(designM,designP,0.0,Py,tauType);
            P1(1) = p;
            M1(1) = 0.0;
            P2(1) = p;
            M2(1) = m2;
            
            % Run nonproportional analyses at different axial loads
            Ps = linspace(P1(1),0,numPoints);
            for i = 2:length(Ps)
                p = Ps(i);
                tau = AISC_tau(p/Py,tauType);
                obj2 = obj.get_copy_with_reduced_stiffness(tau*columnStiffnessReduction);
                id = interactionDiagram2d(designM,designP);
                m2 = id.findXgivenY(p,'Positive');
                m1 = obj2.determineAppliedMoment(p,m2);
                P1(i) = p;
                M1(i) = m1;
                P2(i) = p;
                M2(i) = m2;
            end
        end
        function err = errorEIeff_d(obj,P1,x1,vmax,EIeff)
            obj2 = obj.get_copy_with_reduced_stiffness(EIeff/obj.EI);
            err = vmax-obj2.maxSecondOrderDisplacement(P1,x1);
        end
        function err = errorEIeff_f(obj,P1,x1,vmax,EIeff)
            obj2 = obj.get_copy_with_reduced_stiffness(EIeff/obj.EI);
            err = vmax-obj2.maxSecondOrderMoment(P1,x1);
        end
        function [P2e,M2e] = impliedInteraction(obj,P1_FN,M1_FN,...
            notionalLoadObject,...
            columnStiffnessReduction,beamStiffnessReduction,tauType,Py)
        
            P2e = nan(size(P1_FN));
            M2e = nan(size(M1_FN));
            for iPoint = 1:length(P2e)
                P = P1_FN(iPoint);
                                
                tau = AISC_tau(P/Py,tauType);
                obj2 = obj.get_copy_with_reduced_stiffness(tau*columnStiffnessReduction);
                
                if P > obj2.eulerLoad
                    continue
                end
                
                M = obj2.determineAppliedLoadFromM1(P,M1_FN(iPoint));
                
                P2e(iPoint) = P;
                M2e(iPoint) = obj2.maxSecondOrderMoment(P,M);
            end
        end
        function [drift_1,drift_2,driftRatio] = computeDrifts(obj,P1,M1)
            assert(isequal(size(P1),size(M1)),'The size of P1 and M1 must be the same');
            drift_1 = zeros(size(P1));
            drift_2 = zeros(size(P1));
            driftRatio = ones(size(P1));
        end
        function p = determineLoadForMomentRatio(obj,momentRatio,...
                columnStiffnessReduction,beamStiffnessReduction,tauType,Py_tau)
            
            % Check for quick exit
            tau = AISC_tau(1,tauType);
            if tau ~= 0
                obj2 = obj.get_copy_with_reduced_stiffness(tau*columnStiffnessReduction);
                if obj2.secondOrderMomentRatio(-Py_tau) < momentRatio;
                    p = Py_tau;
                    return
                end
            end
               
            % Estimate the necessary lateral load
            Pest   = 0.95*min(Py_tau,obj.eulerLoad);
            
            % fsolve to determine lateral load
            options = struct;
            options.Display = 'off';
            [p,~,exitflag] =  fsolve(...
                @(P)errorPforMomentRatio(obj,P,momentRatio,columnStiffnessReduction,tauType,Py_tau),...
                Pest,options);
            if exitflag <= 0
                %error('fsolve could not find solution');
                fprintf('fsolve could not find solution\n');
                p = NaN;
            end
        end 
        function obj2 = get_copy_with_reduced_stiffness(obj,column_reduction)
            obj2 = BenchmarkAnalysis2d_Elastic_Sidesway_Inhibited(...
                obj.EI*column_reduction,...
                obj.L,...
                obj.beta,...
                obj.delta0);
            obj2.includeInitialGeometricImperfections = obj.includeInitialGeometricImperfections;
        end        
    end

    methods (Static)
        function t = type()
            t = 'SideswayInhibited';
        end
        function t = type2()
            t = 'I';
        end
    end
end



% function x = errorPeakLoad(obj,designM,designP,P,M1)
% if P < -obj.eulerLoad
%     x = NaN;
% else
%     M2 = obj.maxSecondOrderMoment(P,M1);
%     x = interactionDiagram2d.checkPoints(designM,designP,M2,P)-1;
% end
% end

function x = errorPeakLoad(obj,designM,designP,P,M1,Py,tauType)
tau = AISC_tau(-P/Py,tauType);
obj2 = obj.get_copy_with_reduced_stiffness(tau);
if P < -obj2.eulerLoad
    x = NaN;
else
    M2 = obj2.maxSecondOrderMoment(P,M1);
    id = interactionDiagram2d(designM,designP);
    x  = id.checkPoints(M2,P)-1;
end
end

function x = errorPforMomentRatio(obj,P,momentRatio,columnStiffnessReduction,tauType,Py_tau)

if P < 0
    x = NaN;
else
    tau = AISC_tau(P/Py_tau,tauType);
    obj2 = obj.get_copy_with_reduced_stiffness(tau*columnStiffnessReduction);

    if P >= obj2.eulerLoad
        x = NaN;
    else
        x = obj2.secondOrderMomentRatio(-P)-momentRatio;
    end
end

end

function x = errorCriticalLoadWithTau(obj,P,Py,tauType)
tau = AISC_tau(P/Py,tauType);
obj2 = obj.get_copy_with_reduced_stiffness(tau);
x = -P - obj2.eulerLoad;
end
