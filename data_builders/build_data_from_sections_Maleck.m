function data = build_data_from_sections_Maleck(sections,include_UR)

% Frame set based on:
% Maleck, A. E. (2001). "Second-Order Inelastic and Modified Elastic Analysis 
% and Design Evaluation of Planar Steel Frames." Ph.D. Dissertation, Georgia 
% Institute of Technology, Atlanta, Georgia.

if nargin < 2
    include_UR = true;
end

if include_UR
    num_frames = 23;
else
    num_frames = 18;
end

% Define Data
num_sections = length(sections);
data(num_frames*num_sections) = struct;

i = 1;
for iSection = 1:num_sections
    
    axis = sections(iSection).axis;
    r    = sections(iSection).section.r(axis);
    EI   = sections(iSection).section.E*sections(iSection).section.I(axis);

    % Copy data from sections structure
    for j = i:(i+num_frames-1)
        f = fields(sections(iSection));
        for k = 1:length(f)
            data(j).(f{k}) = sections(iSection).(f{k});
        end
    end
    
    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 20*r;
    data(i).kqtop       = Inf;
    data(i).kqbot       = 0;
    data(i).gamma       = 3;
    data(i).frame_name  = 'UP_20_G0_a3'; 
    data(i).frame_index = 1;
    i = i+1;
    
    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 20*r;
    data(i).kqtop       = (6*EI)/(1*data(i).L);
    data(i).kqbot       = 0;
    data(i).gamma       = 3;
    data(i).frame_name  = 'UP_20_G1_a3'; 
    data(i).frame_index = 2;
    i = i+1;
    
    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 40*r;
    data(i).kqtop       = Inf;
    data(i).kqbot       = 0;
    data(i).gamma       = 2;
    data(i).frame_name  = 'UP_40_G0_a2'; 
    data(i).frame_index = 3;
    i = i+1;
    
    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 40*r;
    data(i).kqtop       = (6*EI)/(1*data(i).L);
    data(i).kqbot       = 0;
    data(i).gamma       = 2;
    data(i).frame_name  = 'UP_40_G1_a2'; 
    data(i).frame_index = 4;
    i = i+1;

    if include_UR
        data(i).frame_type  = 'Sidesway_Uninhibited';
        data(i).L           = 40*r;
        data(i).kqtop       = Inf;
        data(i).kqbot       = Inf;
        data(i).gamma       = 2;
        data(i).frame_name  = 'UR_40_G0_a2'; 
        data(i).frame_index = 5;
        i = i+1;

        data(i).frame_type  = 'Sidesway_Uninhibited';
        data(i).L           = 40*r;
        data(i).kqtop       = Inf;
        data(i).kqbot       = Inf;
        data(i).gamma       = 3;
        data(i).frame_name  = 'UR_40_G0_a3'; 
        data(i).frame_index = 6;
        i = i+1;

        data(i).frame_type  = 'Sidesway_Uninhibited';
        data(i).L           = 80*r;
        data(i).kqtop       = Inf;
        data(i).kqbot       = Inf;
        data(i).gamma       = 1;
        data(i).frame_name  = 'UR_80_G0_a1'; 
        data(i).frame_index = 7;
        i = i+1;

        data(i).frame_type  = 'Sidesway_Uninhibited';
        data(i).L           = 80*r;
        data(i).kqtop       = Inf;
        data(i).kqbot       = Inf;
        data(i).gamma       = 2;
        data(i).frame_name  = 'UR_80_G0_a2'; 
        data(i).frame_index = 8;
        i = i+1;

        data(i).frame_type  = 'Sidesway_Uninhibited';
        data(i).L           = 80*r;
        data(i).kqtop       = Inf;
        data(i).kqbot       = Inf;
        data(i).gamma       = 3;
        data(i).frame_name  = 'UR_80_G0_a3'; 
        data(i).frame_index = 9;
        i = i+1;
    end
    
    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 20*r;
    data(i).kqtop       = Inf;
    data(i).kqbot       = 0;
    data(i).gamma       = 0;
    data(i).frame_name  = 'SP_20_G0'; 
    data(i).frame_index = 10;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 40*r;
    data(i).kqtop       = Inf;
    data(i).kqbot       = 0;
    data(i).gamma       = 0;
    data(i).frame_name  = 'SP_40_G0'; 
    data(i).frame_index = 11;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 40*r;
    data(i).kqtop       = (6*EI)/(3*data(i).L);
    data(i).kqbot       = 0;
    data(i).gamma       = 0;
    data(i).frame_name  = 'SP_40_G3'; 
    data(i).frame_index = 12;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 60*r;
    data(i).kqtop       = Inf;
    data(i).kqbot       = 0;
    data(i).gamma       = 0;
    data(i).frame_name  = 'SP_60_G0'; 
    data(i).frame_index = 13;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 80*r;
    data(i).kqtop       = Inf;
    data(i).kqbot       = 0;
    data(i).gamma       = 0;
    data(i).frame_name  = 'SP_80_G0'; 
    data(i).frame_index = 14;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 80*r;
    data(i).kqtop       = (6*EI)/(3*data(i).L);
    data(i).kqbot       = 0;
    data(i).gamma       = 0;
    data(i).frame_name  = 'SP_80_G3'; 
    data(i).frame_index = 15;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 40*r;
    data(i).kqtop       = Inf;
    data(i).kqbot       = Inf;
    data(i).gamma       = 0;
    data(i).frame_name  = 'SR_40_G0';
    data(i).frame_index = 16;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 40*r;
    data(i).kqtop       = (6*EI)/(3*data(i).L);
    data(i).kqbot       = (6*EI)/(3*data(i).L);
    data(i).gamma       = 0;
    data(i).frame_name  = 'SR_40_G3';
    data(i).frame_index = 17;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 80*r;
    data(i).kqtop       = Inf;
    data(i).kqbot       = Inf;
    data(i).gamma       = 0;
    data(i).frame_name  = 'SR_80_G0';
    data(i).frame_index = 18;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Uninhibited';
    data(i).L           = 80*r;
    data(i).kqtop       = (6*EI)/(3*data(i).L);
    data(i).kqbot       = (6*EI)/(3*data(i).L);
    data(i).gamma       = 0;
    data(i).frame_name  = 'SR_80_G3';
    data(i).frame_index = 19;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Inhibited';
    data(i).L           = 80*r;
    data(i).beta        = 1.0;
    data(i).frame_name  = 'SCB_80';
    data(i).frame_index = 20;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Inhibited';
    data(i).L           = 120*r;
    data(i).beta        = 1.0;
    data(i).frame_name  = 'SCB_120';
    data(i).frame_index = 21;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Inhibited';
    data(i).L           = 80*r;
    data(i).beta        = -0.5;
    data(i).frame_name  = 'DCB_80';
    data(i).frame_index = 22;
    i = i+1;

    data(i).frame_type  = 'Sidesway_Inhibited';
    data(i).L           = 120*r;
    data(i).beta        = -0.5;
    data(i).frame_name  = 'DCB_120';
    data(i).frame_index = 23;
    i = i+1;
end

% Define initial geometric imperfections
for i = 1:length(data)
    switch data(i).frame_type
        case 'Sidesway_Uninhibited'
            data(i).Delta0 = data(i).L/500;
            data(i).delta0 = data(i).L/1000;
        case 'Sidesway_Inhibited'
            data(i).Delta0 = 0;
            data(i).delta0 = data(i).L/1000;
        otherwise
            error('Unknwon frame_type: %s',data(i).frame_type)
    end
end

end

