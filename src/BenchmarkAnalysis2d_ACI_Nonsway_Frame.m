classdef BenchmarkAnalysis2d_ACI_Nonsway_Frame < BenchmarkAnalysis2d_ACI_Base

    properties
        beta % End moment ratio
    end
    
    methods
        function obj = BenchmarkAnalysis2d_ACI_Nonsway_Frame(section,axis,L,beta)
            if nargin == 1
                data = section;
                obj.section = data.section;
                obj.axis    = data.axis;
                obj.L       = data.L;
                obj.beta    = data.beta;
            else
                obj.section = section;
                obj.axis    = axis;
                obj.L       = L;
                obj.beta    = beta;
            end
            
            num_points = 50;            
            sc = obj.section.strainCompatibilityAciObject;
            switch lower(obj.axis)
                case {'x','z'}
                    [P,M,~,et] = sc.interactionSweep(0,num_points);
                case 'y'
                    [P,~,M,et] = sc.interactionSweep(pi/2,num_points);
                otherwise
                    error('Bad axis');
            end
            phi = obj.section.phi(et);
            
            % Store results
            obj.InteractionDiagram_P   = P;
            obj.InteractionDiagram_M   = M;
            obj.InteractionDiagram_et  = et;
            obj.InteractionDiagram_phi = phi;
            obj.id2d_nom = interactionDiagram2d(M,P);
            obj.id2d_red = interactionDiagram2d(phi.*M,phi.*P);
        end
        function cm = Cm(obj)
            if obj.beta < -1
                m1m2 = 1/obj.beta;
            elseif obj.beta < 1
                m1m2 = -obj.beta;
            else
                m1m2 = -1/obj.beta;
            end
            cm = 0.6-0.4*m1m2;
        end
        function results = run_to_peak_proportional(obj,e)
            if nargin < 2
                e = 0;
            end
            
            % Compression is negative
            num_points = 1000;
            
            % Run analysis
            if strcmpi(obj.EIeff_type,'c') || strcmpi(obj.EIeff_type,'6.6.4.4.4c')
                assert(e == 0,'Only implemented for e = 0')
                
                P  = nan(1,num_points);
                M1 = zeros(1,num_points);
                M2 = linspace(0,1.01*max(obj.InteractionDiagram_M),num_points);
                delta = nan(1,num_points);
                
                for i = 1:num_points
                    opts = optimoptions('fsolve',...
                        'Display','off');
                    if i == 1
                        P_guess = 0;
                    else
                        P_guess = P(i-1);
                    end 
                    [P(i),~,exitflag,output] = fsolve(@(p)error_EI_type_c(obj,p,0,M2(i)),P_guess,opts);

                    if exitflag <= 0
                        fprintf('%s\n',output.message);
                        error('fsolve unable to find a solution')
                    end
                    
                    % Compute Delta
                    Pc = pi^2*obj.EIeff('nonsway',M2(i),P(i))/obj.L^2;
                    if obj.include_stiffness_reduction
                        delta(i) = max(obj.Cm./(1-(-P(i))/(0.75*Pc)),1);
                    else
                        delta(i) = max(obj.Cm./(1-(-P(i))/Pc),1);
                    end
                    
                    % Check stopping points
                    if P(i) > 0.98*min(P)
                        % Stability limit reached
                        break
                    end
                    
                    if obj.include_strength_reduction
                        if obj.id2d_red.checkPoints(M2(i),P(i)) > 1
                            % Cross section strength reached
                            break
                        end
                    else
                        if obj.id2d_nom.checkPoints(M2(i),P(i)) > 1
                            % Cross section strength reached
                            break
                        end
                    end
                end
            else
                Pc  = pi^2*obj.EIeff('nonsway')/obj.L^2;

                if obj.include_stiffness_reduction
                    P  = linspace(0,-min(0.851*obj.section.Po,0.749*Pc),1000);
                    M1 = P*max(obj.minimum_eccentricity,e);
                    delta = max(obj.Cm./(1-(-P)/(0.75*Pc)),1);
                else
                    P  = linspace(0,-min(0.851*obj.section.Po,0.99*Pc),1000);
                    M1 = P*max(obj.minimum_eccentricity,e);
                    delta = max(obj.Cm./(1-(-P)/Pc),1);
                end
                M2 = delta.*M1;
            end
            
            results = struct;
            results.path.P  = P;
            results.path.M1 = M1;
            results.path.M2 = M2;
            results.path.second_order_moment_ratio = delta;
            
            % Determine Peak
            results.limit_point = obj.find_limit_point(results,'P');
        end
        function results = run_to_peak_nonproportional(obj,Pa)
            % Compression is negative
            num_points = 1000;
            
            % Run analysis
            if strcmpi(obj.EIeff_type,'c') || strcmpi(obj.EIeff_type,'6.6.4.4.4c')
                
                P  = Pa*ones(1,num_points);
                M1 = nan(1,num_points);
                delta = nan(1,num_points);
                
                if obj.include_strength_reduction
                    M2_max = obj.id2d_red.findXgivenY(Pa,'pos');
                else
                    M2_max = obj.id2d_nom.findXgivenY(Pa,'pos');
                end
                M2     = linspace(0,M2_max,num_points);
                
                for i = 1:num_points
                    Pc = pi^2*obj.EIeff('nonsway',M2(i),P(i))/obj.L^2;
                    if obj.include_stiffness_reduction
                        if P > -0.75*Pc
                            delta(i) = max(obj.Cm./(1-(-P(i))/(0.75*Pc)),1);
                            M1(i)    = M2(i)/delta(i);
                        end
                    else
                        if P > -Pc
                            delta(i) = max(obj.Cm./(1-(-P(i))/Pc),1);
                            M1(i)    = M2(i)/delta(i);
                        end
                    end
                end
                
            else
                Pc = pi^2*obj.EIeff('nonsway')/obj.L^2;           
                P  = Pa*ones(1,num_points);
                M1 = linspace(0,max(obj.InteractionDiagram_M),num_points);
                if obj.include_stiffness_reduction
                    delta = max(obj.Cm./(1-(-P)/(0.75*Pc)),1);
                else
                    delta = max(obj.Cm./(1-(-P)/Pc),1);
                end
                M2 = delta.*max(M1,(-P)*obj.minimum_eccentricity);
            end
            
            results = struct;
            results.path.P  = P;
            results.path.M1 = M1;
            results.path.M2 = M2;
            results.path.second_order_moment_ratio = delta;
            
            % Determine Peak
            results.limit_point = obj.find_limit_point(results,'M1');
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


function err = error_EI_type_c(obj,P,M1,M2)
if P > 0
    err = Inf;
    return
end
Pc = pi^2*obj.EIeff('nonsway',M2,P)/obj.L^2;
if obj.include_stiffness_reduction
    if P < -0.75*Pc
        err = Inf;
        return
    end
    delta = max(obj.Cm./(1-(-P)/(0.75*Pc)),1);
else
    if P < -Pc
        err = Inf;
        return
    end
    delta = max(obj.Cm./(1-(-P)/Pc),1);
end
M1n = max(M1,(-P)*obj.minimum_eccentricity);
M2n = delta.*M1n;
err = M2n-M2;
end
