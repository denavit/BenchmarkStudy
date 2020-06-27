classdef BenchmarkAnalysis2d_ACI_Base
        
    properties
        section         % reinforced concrete section
        axis            % bending axis
        L               % Length of the member
        beta_ns = 0     % 
        beta_dns = 0    % sustained axial load ratio
        
        InteractionDiagram_P
        InteractionDiagram_M
        InteractionDiagram_et
        InteractionDiagram_phi
        id2d_nom
        id2d_red
        
        include_strength_reduction = true;
        include_stiffness_reduction = true;
        second_order_moment_ratio_limit = 1.4;
        EIeff_type = 'a';        
    end
     
    methods    
        function [P1,M1,P2,M2,type] = design_interaction(obj,num_points)
            
            % Initilize results
            P1 = nan(num_points,1);
            M1 = nan(num_points,1);
            P2 = nan(num_points,1);
            M2 = nan(num_points,1);            
            type = cell(num_points,1);
            
            % Determine peak load
            results = obj.run_to_peak_proportional(0);
            P1(1) = results.limit_point.P;
            M1(1) = results.limit_point.M1;
            P2(1) = results.limit_point.P;
            M2(1) = results.limit_point.M2;
            type{1} = results.limit_point.type;
            
            % Determine nonproportional analyses at different axial loads
            Ps = [P1(1) linspace(0.999*P1(1),0,num_points-1)];
            for i = 2:length(Ps)
                results = obj.run_to_peak_nonproportional(Ps(i));
                P1(i) = results.limit_point.P;
                M1(i) = results.limit_point.M1;
                P2(i) = results.limit_point.P;
                M2(i) = results.limit_point.M2;
                type{i} = results.limit_point.type;
            end
            
            if nargout < 5
                clear type
            end
        end
        function ei = EIeff(obj,sway,Mu,Pu)
            Ec = obj.section.Ec;
            Ig = obj.section.Ig(obj.axis);
                    
            switch obj.EIeff_type
                case {'a','6.6.4.4.4a'}
                    ei = 0.4*Ec*Ig;
                case {'b','6.6.4.4.4b'}
                    Es  = obj.section.Es;
                    Isr = obj.section.Isr(obj.axis); 
                    ei  = 0.2*Ec*Ig+Es*Isr;
                case {'c','6.6.4.4.4c'}
                    % Compression is negative
                    assert(Pu<=0,'Function only valid for compressive loads: P = %g',Pu)
            
                    Ast = obj.section.Asr;
                    Ag  = obj.section.Ag;
                    h   = obj.section.depth(obj.axis);
                    Po  = obj.section.Po;
                    if Mu == 0 && Pu == 0
                        I = (0.80+25*Ast/Ag)*Ig;
                    else
                        I = (0.80+25*Ast/Ag)*(1-Mu/((-Pu)*h)-0.5*(-Pu)/Po)*Ig;
                    end
                    if I > 0.875*Ig
                        I = 0.875*Ig;
                    end
                    if I < 0.35*Ig
                        I = 0.35*Ig;
                    end
                    ei = Ec*I;
                otherwise
                    error('Unknown type: %s',obj.EIeff_type)
            end
            
            switch lower(sway)
                case 'nonsway'
                    ei = ei/(1+obj.beta_dns);
                case 'sway'
                    ei = ei/(1+obj.beta_ns);
                otherwise
                    error('Unknown sway type: %s',obj.sway)
            end
        end        
        function emin = minimum_eccentricity(obj)
            switch obj.section.units
                case 'US'
                    e1 = 0.6; % inches
                otherwise
                    error('Unknown unit system: %s',obj.section.units)
            end
            emin = e1 + 0.03*obj.section.depth(obj.axis);
        end        
        function limit_point = find_limit_point(obj,results,applied_load_type)
            
            % Stability Limit
            switch applied_load_type
                case 'P'
                    [P,ind] = min(results.path.P);
                    limit_stability.P    = P;
                    limit_stability.M1   = results.path.M1(ind);
                    limit_stability.M2   = results.path.M2(ind);
                    limit_stability.time = ind;
                    limit_stability.type = 'Stability Limit';
                case 'M1'
                    [M1,ind] = max(results.path.M1);
                    limit_stability.P    = results.path.P(ind);
                    limit_stability.M1   = M1;
                    limit_stability.M2   = results.path.M2(ind);
                    limit_stability.time = ind;
                    limit_stability.type = 'Stability Limit';
                otherwise
                    error('Unknown applied_load_type: %s',applied_load_type);
            end
                    
            % Cross Section Strength Limit
            limit_cross_section_strength = struct;
            if obj.include_strength_reduction
                [M2,P,ind,x] = obj.id2d_red.findIntersection(results.path.M2,results.path.P);
            else
                [M2,P,ind,x] = obj.id2d_nom.findIntersection(results.path.M2,results.path.P);
            end
            if isempty(ind)
                limit_cross_section_strength.P    = nan;
                limit_cross_section_strength.M1   = nan;
                limit_cross_section_strength.M2   = nan;
                limit_cross_section_strength.time = Inf;
                limit_cross_section_strength.type = '';
            else
                limit_cross_section_strength.P    = P;
                limit_cross_section_strength.M1   = interpolate_vector(results.path.M1,ind,x);
                limit_cross_section_strength.M2   = M2;
                limit_cross_section_strength.time = ind+x;
                limit_cross_section_strength.type = 'Cross Section Strength Limit';
            end
            
            % Second Order Ratio Limit
            limit_second_order_ratio = struct;
            if isempty(obj.second_order_moment_ratio_limit)
                limit_second_order_ratio.P    = nan;
                limit_second_order_ratio.M1   = nan;
                limit_second_order_ratio.M2   = nan;
                limit_second_order_ratio.time = Inf;
                limit_second_order_ratio.type = '';
            else
                [ind,x] = find_limit_point_in_vector(results.path.second_order_moment_ratio,obj.second_order_moment_ratio_limit);
                if isempty(ind)
                    limit_second_order_ratio.P    = nan;
                    limit_second_order_ratio.M1   = nan;
                    limit_second_order_ratio.M2   = nan;
                    limit_second_order_ratio.time = Inf;
                    limit_second_order_ratio.type = '';
                else
                    limit_second_order_ratio.P    = interpolate_vector(results.path.P,ind,x);
                    limit_second_order_ratio.M1   = interpolate_vector(results.path.M1,ind,x);
                    limit_second_order_ratio.M2   = interpolate_vector(results.path.M2,ind,x);
                    limit_second_order_ratio.time = ind+x;
                    limit_second_order_ratio.type = 'Second Order Ratio Limit';
                end
            end
            
            % Find Controlling Limit
            limit_point_times = [
                limit_stability.time
                limit_cross_section_strength.time
                limit_second_order_ratio.time      ];
            controlling_limit_point_time = min(limit_point_times);
            
            if limit_stability.time == controlling_limit_point_time
                limit_point = limit_stability;
            elseif limit_cross_section_strength.time == controlling_limit_point_time
                limit_point = limit_cross_section_strength;
            elseif limit_second_order_ratio.time == controlling_limit_point_time
                limit_point = limit_second_order_ratio;
            else
                error('error');
            end
        end
    end
end
 
