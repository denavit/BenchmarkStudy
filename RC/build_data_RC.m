function data = build_data_RC()

% Concrete Strength
fc          = [4 6 8 10 12];
Ec          = 57*sqrt(1000*fc);

% Steel Strength
fy          = [60  80];
fu          = [90 105];

% Cross Sections
[~,section_data] = csvread2('section_data_RC.csv');
num_sections = length(section_data.shape);

% Slenderness
L_over_H    = [5 10 15 20 25 30 35 40];

% Boundary Conditions
num_boundary_conditions = 9;
Psi         = [0 1 2 4 8];
beta        = [-0.5 0.0 0.5 1.0];

% Constants
axis = 'x';
fyt = 60;
Es  = 29000;

% Build Sections
data(length(fc)*length(fy)*num_sections*length(L_over_H)*num_boundary_conditions) = struct;

iData = 1;
for i1 = 1:length(fc)
    for i2 = 1:length(fy)
        for i3 = 1:num_sections
            for i4 = 1:length(L_over_H)
                for i5 = 1:num_boundary_conditions
                    data(iData).section_type            = 'RC';
                    data(iData).units                   = 'US';
                    data(iData).axis                    = axis;
                    data(iData).section_name            = sprintf('Section %i',i3);

                    switch section_data.shape{i3}
                        case 'square' 
                            % Define gross cross section
                            H = section_data.B(i3);
                            B = section_data.B(i3);
                            conc_cross_section = Rectangle_Shape(H,B);

                            longitudinal_bar_size   = reinf_bar_lookup(section_data.longitudinal_bar_size{i3});
                            transverse_bar_size     = reinf_bar_lookup(section_data.transverse_bar_size{i3});
                            
                            % Define longitudinal reinforcement
                            Bc = B - 2*section_data.cover(i3) - 2*transverse_bar_size.diameter - longitudinal_bar_size.diameter;
                            Hc = H - 2*section_data.cover(i3) - 2*transverse_bar_size.diameter - longitudinal_bar_size.diameter;
                            num_bars = sscanf(section_data.longitudinal_config{i3},'%ix-%iy');
                            nbx = num_bars(1);
                            nby = num_bars(2);
                            reinforcement = reinf_rect(Bc,Hc,0,0,nbx,nby,longitudinal_bar_size.area,longitudinal_bar_size.diameter);

                            % Define cross section
                            data(iData).section = RC(fc(i1),fy(i2),conc_cross_section,reinforcement,'US');
                            data(iData).section.transverse_reinf_type = section_data.transverse_reinf_type{i3};
                            data(iData).section.Ec = Ec(i1);

                            switch section_data.transverse_config{i3}
                                case 'A'
                                    nLegX = 2;
                                    nLegY = 2;
                                case 'B'
                                    nLegX = 3;
                                    nLegY = 3;
                                case 'C'
                                    nLegX = 4;
                                    nLegY = 4;
                                case 'D'
                                    nLegX = 2;
                                    nLegY = 3;
                                case {'E','F'}
                                    nLegX = 2;
                                    nLegY = 4;
                                case 'D*'
                                    nLegX = 3;
                                    nLegY = 2;
                                case {'E*','F*'}
                                    nLegX = 4;
                                    nLegY = 2;
                                otherwise
                                    error('Unknown transverse config: %s',section_data.transverse_config{i3});
                            end
                            
                            % Define other data
                            data(iData).shape               = 'square';
                            data(iData).B                   = B;
                            data(iData).H                   = H;
                            data(iData).cover               = section_data.cover(i3);
                            data(iData).fc                  = fc(i1);
                            data(iData).Ec                  = Ec(i1);
                            data(iData).L                   = L_over_H(i4)*H;
                            data(iData).nBarX               = nbx;
                            data(iData).nBarY               = nby;
                            data(iData).Es                  = Es;
                            data(iData).fy                  = fy(i2);
                            data(iData).fu                  = fu(i2);
                            data(iData).db                  = longitudinal_bar_size.diameter;
                            data(iData).Ab                  = longitudinal_bar_size.area;
                            data(iData).fyt                 = fyt;
                            data(iData).dbt                 = transverse_bar_size.diameter;
                            data(iData).Abt                 = transverse_bar_size.area; 
                            data(iData).nLegX               = nLegX;
                            data(iData).nLegY               = nLegY;
                            data(iData).s                   = str2double(section_data.s{i3});
                            data(iData).longitudinal_config     = section_data.longitudinal_config{i3};
                            data(iData).transverse_config       = section_data.transverse_config{i3};
                            data(iData).longitudinal_bar_size   = section_data.longitudinal_bar_size{i3};
                            data(iData).transverse_bar_size     = section_data.transverse_bar_size{i3};                   
                    
                        case 'circle'
                            % Define gross cross section
                            D = section_data.B(i3);
                            conc_cross_section = Circle_Shape(D);

                            longitudinal_bar_size   = reinf_bar_lookup(section_data.longitudinal_bar_size{i3});
                                                        
                            % transverse reinf. 
                            switch section_data.transverse_reinf_type{i3}
                                case 'ties'
                                    transverse_bar_size = reinf_bar_lookup(section_data.transverse_bar_size{i3});
                                    s = str2double(section_data.s{i3});
                                case 'spiral'
                                    Dc  = D - 2*section_data.cover(i3);
                                    Ag  = pi/4*D^2;
                                    Ach = pi/4*Dc^2;
                                    
                                    % Try #3 bars
                                    transverse_bar_size = reinf_bar_lookup('#3');
                                    rho_s_min =0.45*(Ag/Ach-1)*fc(i1)/fyt;
                                    s = min([3 4*transverse_bar_size.area/rho_s_min/Dc]);
                                    
                                    if s < 1.25
                                        % Try #4 bars
                                        transverse_bar_size = reinf_bar_lookup('#4');
                                        rho_s_min =0.45*(Ag/Ach-1)*fc(i1)/fyt;
                                        s = min([3 4*transverse_bar_size.area/rho_s_min/Dc]);
                                        
                                        if s < 1.25
                                            error('Spacing too small with #4 bars');
                                        end
                                    end
                                otherwise
                                    error('Unknown transverse reinf type: %s',section_data.transverse_reinf_type{i3});
                            end
                            
                            % Define longitudinal reinforcement
                            rc = 0.5*D - section_data.cover(i3) - transverse_bar_size.diameter - 0.5*longitudinal_bar_size.diameter;
                            num_bars = sscanf(section_data.longitudinal_config{i3},'%i bars');
                            reinforcement = reinf_circ(0,0,rc,num_bars,longitudinal_bar_size.area,longitudinal_bar_size.diameter);

                            % Define cross section
                            data(iData).section = RC(fc(i1),fy(i2),conc_cross_section,reinforcement,'US');
                            data(iData).section.transverse_reinf_type = section_data.transverse_reinf_type{i3};
                            data(iData).section.Ec = Ec(i1);

                            % Define other data
                            data(iData).shape               = 'circle';
                            data(iData).D                   = D;                        
                            data(iData).L                   = L_over_H(i4)*D;
                            data(iData).cover               = section_data.cover(i3);
                            data(iData).fc                  = fc(i1);
                            data(iData).Ec                  = Ec(i1);
                            data(iData).num_bars            = num_bars;
                            data(iData).Es                  = Es;
                            data(iData).fy                  = fy(i2);
                            data(iData).fu                  = fu(i2);
                            data(iData).db                  = longitudinal_bar_size.diameter;
                            data(iData).Ab                  = longitudinal_bar_size.area;
                            data(iData).fyt                 = fyt;
                            data(iData).dbt                 = transverse_bar_size.diameter;
                            data(iData).Abt                 = transverse_bar_size.area;
                            data(iData).s                   = s;
                            data(iData).transverse_reinf_type = section_data.transverse_reinf_type{i3};
                            
                        otherwise
                            error('Unknown shape: %s',section_data.shape{i3})
                    end
                    
                    if i5 <= 5
                        data(iData).frame_type = 'Sidesway_Uninhibited';
                        data(iData).Psi = Psi(i5);
                        EcIgc_over_L = data(iData).Ec * data(iData).section.Ig(data(iData).axis) / data(iData).L;
                        data(iData).EcIgb_over_Lb = 2*EcIgc_over_L/Psi(i5);
                        data(iData).kqtop = 6*0.4*data(iData).EcIgb_over_Lb; % Without stiffness reduction
                        data(iData).kqbot = 6*0.4*data(iData).EcIgb_over_Lb; % Without stiffness reduction
                        data(iData).gamma = 0.0;
                        data(iData).delta0 = data(iData).L/1000;
                        data(iData).Delta0 = data(iData).L/500;
                    else
                        data(iData).frame_type = 'Sidesway_Inhibited';
                        data(iData).beta = beta(i5-5);
                        data(iData).delta0 = data(iData).L/1000;
                    end
                    
                    iData = iData+1;
                end
            end
        end
    end
end        
        

end