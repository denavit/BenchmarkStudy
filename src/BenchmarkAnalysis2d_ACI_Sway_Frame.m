classdef BenchmarkAnalysis2d_ACI_Sway_Frame < BenchmarkAnalysis2d_ACI_Base

    properties
        EcIgb_over_Lb_top % Ratio of gross moment of itertia to length of the beam at the top
        EcIgb_over_Lb_bot % Ratio of gross moment of itertia to length of the beam at the bottom
    end
    
    methods
        function obj = BenchmarkAnalysis2d_ACI_Sway_Frame(section,axis,L,EcIgb_over_Lb)
            if nargin == 1
                data = section;
                obj.section = data.section;
                obj.axis    = data.axis;
                obj.L       = data.L;
                obj.EcIgb_over_Lb_top = data.EcIgb_over_Lb;
                obj.EcIgb_over_Lb_bot = data.EcIgb_over_Lb;
            else
                obj.section = section;
                obj.axis    = axis;
                obj.L       = L;
                obj.EcIgb_over_Lb_top = EcIgb_over_Lb;
                obj.EcIgb_over_Lb_bot = EcIgb_over_Lb;
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
            if obj.EcIgb_over_Lb_top == obj.EcIgb_over_Lb_bot
                m1m2 = 1;
            else
                error('Not yet implemented');
            end
            cm = 0.6-0.4*m1m2;
        end
        function k = K(obj)
            % The stiffness of the column is taken as 0.7EcIg 
            % The stiffness of the beams is taken as 0.35EcIg 
            % per Section R6.6.4.4.3 of ACI 318-19
            Ec = obj.section.Ec;
            Ig = obj.section.Ig(obj.axis);
            G_top = (0.7*Ec*Ig/obj.L)/(0.35*obj.EcIgb_over_Lb_top);
            G_bot = (0.7*Ec*Ig/obj.L)/(0.35*obj.EcIgb_over_Lb_bot);
            k = nomographK_sidesway_uninhibited(G_top,G_bot);
        end
        function results = run_to_peak_proportional(obj,e)
            if nargin < 2
                e = 0;
            end
            
            assert(e == 0,'Only implemented for e = 0')
            
            % Compression is negative
            num_points = 1000;
            
            % Run analysis
            if strcmpi(obj.EIeff_type,'c') || strcmpi(obj.EIeff_type,'6.6.4.4.4c')
                                
                P  = nan(1,num_points);
                Mc = linspace(0,1.01*max(obj.InteractionDiagram_M),num_points);
                second_order_moment_ratio = nan(1,num_points);
                
                for i = 1:num_points
                    opts = optimoptions('fsolve',...
                        'Display','off');
                    if i == 1
                        P_guess = 0;
                    else
                        P_guess = P(i-1);
                    end 
                    [P(i),~,exitflag,output] = fsolve(@(p)error_EI_type_c(obj,p,0,Mc(i)),P_guess,opts);

                    if exitflag <= 0
                        fprintf('%s\n',output.message);
                        error('fsolve unable to find a solution')
                    end
                    
                    % Compute Second Order Moment Ratio
                    EIeff   = obj.EIeff('sway',Mc(i),P(i));
                    Pc      = pi^2*EIeff/(obj.K*obj.L)^2;
                    if obj.include_stiffness_reduction
                        delta   = max(obj.Cm./(1-(-P(i))/(0.75*Pc)),1);
                        delta_s = max(1./(1-(-P(i))/(0.75*Pc)),1);
                    else
                        delta   = max(obj.Cm./(1-(-P(i))/Pc),1);
                        delta_s = max(1./(1-(-P(i))/Pc),1);
                    end
                    second_order_moment_ratio(i) = max(delta,delta_s);
                    
                    % Check stopping points
                    if P(i) > 0.98*min(P)
                        % Stability limit reached
                        break
                    end
                    
                    if obj.include_strength_reduction
                        if obj.id2d_red.checkPoints(Mc(i),P(i)) > 1
                            % Cross section strength reached
                            break
                        end
                    else
                        if obj.id2d_nom.checkPoints(Mc(i),P(i)) > 1
                            % Cross section strength reached
                            break
                        end
                    end
                end
            else
                EIeff   = obj.EIeff('sway');
                Pc      = pi^2*EIeff/(obj.K*obj.L)^2;
                
                if obj.include_stiffness_reduction
                    P = linspace(0,-min(0.851*obj.section.Po,0.749*Pc),num_points);
                    delta = max(obj.Cm./(1-(-P)/(0.75*Pc)),1);
                    delta_s = max(1./(1-(-P)/(0.75*Pc)),1);
                else
                    P = linspace(0,-min(0.851*obj.section.Po,0.99*Pc),num_points);
                    delta = max(obj.Cm./(1-(-P)/Pc),1);
                    delta_s = max(1./(1-(-P)/Pc),1);
                end
                
                % eccentricity is zero, M1s is zero
                M2min = (-P)*obj.minimum_eccentricity;
                Mc = delta.*M2min;    
                second_order_moment_ratio = max(delta,delta_s);
            end
            
            results = struct;
            results.path.P  = P;
            results.path.M1 = zeros(size(P));
            results.path.M2 = Mc;
            results.path.second_order_moment_ratio = second_order_moment_ratio;
            
            % Determine Peak
            results.limit_point = obj.find_limit_point(results,'P');
        end
        function results = run_to_peak_nonproportional(obj,Pa)
            % Compression is negative
            num_points = 1000;
            
            % Run analysis
            if strcmpi(obj.EIeff_type,'c') || strcmpi(obj.EIeff_type,'6.6.4.4.4c')
                
                P       = Pa*ones(1,num_points);
                if obj.include_strength_reduction
                    Mc_max  = obj.id2d_red.findXgivenY(Pa,'pos');
                else
                    Mc_max  = obj.id2d_nom.findXgivenY(Pa,'pos');
                end
                Mc      = linspace(0,Mc_max,num_points);
                
                M2s     = nan(1,num_points);
                
                for i = 1:num_points
                    EIeff   = obj.EIeff('sway',Mc(i),P(i));
                    Pc      = pi^2*EIeff/(obj.K*obj.L)^2;
                    
                    if obj.include_stiffness_reduction
                        if P > -0.75*Pc
                            delta   = max(obj.Cm./(1-(-P(i))/(0.75*Pc)),1);
                            M2      = Mc(i)/delta;

                            delta_s = max(1./(1-(-P(i))/(0.75*Pc)),1);
                            M2s(i)  = M2/delta_s;
                        end
                    else
                        if P > -Pc
                            delta   = max(obj.Cm./(1-(-P(i))/Pc),1);
                            M2      = Mc(i)/delta;

                            delta_s = max(1./(1-(-P(i))/Pc),1);
                            M2s(i)  = M2/delta_s;
                        end
                    end
                end     
            else
                EIeff   = obj.EIeff('sway');
                Pc      = pi^2*EIeff/(obj.K*obj.L)^2;
                
                P       = Pa*ones(1,num_points);
                M2s     = linspace(0,max(obj.InteractionDiagram_M),num_points);
                
                if obj.include_stiffness_reduction
                    delta_s = max(1./(1-(-P)/(0.75*Pc)),1);
                    delta   = max(obj.Cm./(1-(-P)/(0.75*Pc)),1);
                else
                    delta_s = max(1./(1-(-P)/Pc),1);
                    delta   = max(obj.Cm./(1-(-P)/Pc),1);
                end
                M2 = delta_s.*M2s; % M1ns is zero
                Mc = delta.*max(M2,-P*obj.minimum_eccentricity);
            end
            
            results = struct;
            results.path.P  = P;
            results.path.M1 = M2s; % Applied moment
            results.path.M2 = Mc;  % Moment used for design
            results.path.second_order_moment_ratio = Mc./max(M2s,-P*obj.minimum_eccentricity);
            
            % Determine Peak
            results.limit_point = obj.find_limit_point(results,'M1');
        end        
    end

    methods (Static)
        function t = type()
            t = 'SideswayUninhibited';
        end
        function t = type2()
            t = 'U';
        end
    end
end


function err = error_EI_type_c(obj,P,M1,M2)
assert(M1==0,'Only implemented for M1 = 0');
if P > 0
    err = Inf;
    return
end
EIeff   = obj.EIeff('sway',M2,P);
Pc      = pi^2*EIeff/(obj.K*obj.L)^2;

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

M2min   = (-P)*obj.minimum_eccentricity;
Mc      = delta.*M2min;
err     = Mc-M2;
end
