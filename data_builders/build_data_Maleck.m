function data = build_data_Maleck(axis)

% Define Section
shape_data      = steel_shape_lookup('W8x31');
section         = WF(shape_data.d,shape_data.tw,shape_data.bf,shape_data.tf,shape_data.tf,36,'US');
section.Fu      = 58;
section.eu      = 0.20;
section.neglectLocalBuckling = true;
section.neglectLateralTorsionalBuckling = true;

r = section.r(axis);
EI = 29000*section.I(axis);

switch lower(axis)
    case 'strong'
        section_id   = 1;
        section_name = 'SA';
    case 'weak'
        section_id   = 2;
        section_name = 'WA';
    otherwise
        error('Unknown axis: %s',axis)
end

% Define Data
data(23) = struct;
for i = 1:length(data)
    data(i).section         = section;
    data(i).section_type    = 'WF';
    data(i).section_id      = section_id;
    data(i).section_name    = section_name;
    data(i).axis            = axis;
end

data(1).frame_type  = 'Sidesway_Uninhibited';
data(1).L           = 20*r;
data(1).kqtop       = Inf;
data(1).kqbot       = 0;
data(1).gamma       = 3;
data(1).frame_name  = sprintf('UP_%s20_G0_a3',section_name(1)); 

data(2).frame_type  = 'Sidesway_Uninhibited';
data(2).L           = 20*r;
data(2).kqtop       = (6*EI)/(1*data(2).L);
data(2).kqbot       = 0;
data(2).gamma       = 3;
data(2).frame_name  = sprintf('UP_%s20_G1_a3',section_name(1)); 

data(3).frame_type  = 'Sidesway_Uninhibited';
data(3).L           = 40*r;
data(3).kqtop       = Inf;
data(3).kqbot       = 0;
data(3).gamma       = 2;
data(3).frame_name  = sprintf('UP_%s40_G0_a2',section_name(1)); 

data(4).frame_type  = 'Sidesway_Uninhibited';
data(4).L           = 40*r;
data(4).kqtop       = (6*EI)/(1*data(4).L);
data(4).kqbot       = 0;
data(4).gamma       = 2;
data(4).frame_name  = sprintf('UP_%s40_G1_a2',section_name(1)); 

data(5).frame_type  = 'Sidesway_Uninhibited';
data(5).L           = 40*r;
data(5).kqtop       = Inf;
data(5).kqbot       = Inf;
data(5).gamma       = 2;
data(5).frame_name  = sprintf('UR_%s40_G0_a2',section_name(1)); 

data(6).frame_type  = 'Sidesway_Uninhibited';
data(6).L           = 40*r;
data(6).kqtop       = Inf;
data(6).kqbot       = Inf;
data(6).gamma       = 3;
data(6).frame_name  = sprintf('UR_%s40_G0_a3',section_name(1)); 

data(7).frame_type  = 'Sidesway_Uninhibited';
data(7).L           = 80*r;
data(7).kqtop       = Inf;
data(7).kqbot       = Inf;
data(7).gamma       = 1;
data(7).frame_name  = sprintf('UR_%s80_G0_a1',section_name(1)); 

data(8).frame_type  = 'Sidesway_Uninhibited';
data(8).L           = 80*r;
data(8).kqtop       = Inf;
data(8).kqbot       = Inf;
data(8).gamma       = 2;
data(8).frame_name  = sprintf('UR_%s80_G0_a2',section_name(1)); 

data(9).frame_type  = 'Sidesway_Uninhibited';
data(9).L           = 80*r;
data(9).kqtop       = Inf;
data(9).kqbot       = Inf;
data(9).gamma       = 3;
data(9).frame_name  = sprintf('UR_%s80_G0_a3',section_name(1)); 

data(10).frame_type  = 'Sidesway_Uninhibited';
data(10).L           = 20*r;
data(10).kqtop       = Inf;
data(10).kqbot       = 0;
data(10).gamma       = 0;
data(10).frame_name = sprintf('SP_%s20_G0',section_name(1)); 

data(11).frame_type  = 'Sidesway_Uninhibited';
data(11).L           = 40*r;
data(11).kqtop       = Inf;
data(11).kqbot       = 0;
data(11).gamma       = 0;
data(11).frame_name = sprintf('SP_%s40_G0',section_name(1)); 

data(12).frame_type  = 'Sidesway_Uninhibited';
data(12).L           = 40*r;
data(12).kqtop       = (6*EI)/(3*data(12).L);
data(12).kqbot       = 0;
data(12).gamma       = 0;
data(12).frame_name = sprintf('SP_%s40_G3',section_name(1)); 

data(13).frame_type  = 'Sidesway_Uninhibited';
data(13).L           = 60*r;
data(13).kqtop       = Inf;
data(13).kqbot       = 0;
data(13).gamma       = 0;
data(13).frame_name = sprintf('SP_%s60_G0',section_name(1)); 

data(14).frame_type  = 'Sidesway_Uninhibited';
data(14).L           = 80*r;
data(14).kqtop       = Inf;
data(14).kqbot       = 0;
data(14).gamma       = 0;
data(14).frame_name = sprintf('SP_%s80_G0',section_name(1)); 

data(15).frame_type  = 'Sidesway_Uninhibited';
data(15).L           = 80*r;
data(15).kqtop       = (6*EI)/(3*data(15).L);
data(15).kqbot       = 0;
data(15).gamma       = 0;
data(15).frame_name = sprintf('SP_%s80_G3',section_name(1)); 

data(16).frame_type  = 'Sidesway_Uninhibited';
data(16).L           = 40*r;
data(16).kqtop       = Inf;
data(16).kqbot       = Inf;
data(16).gamma       = 0;
data(16).frame_name = sprintf('SR_%s40_G0',section_name(1)); 

data(17).frame_type  = 'Sidesway_Uninhibited';
data(17).L           = 40*r;
data(17).kqtop       = (6*EI)/(3*data(17).L);
data(17).kqbot       = (6*EI)/(3*data(17).L);
data(17).gamma       = 0;
data(17).frame_name = sprintf('SR_%s40_G3',section_name(1)); 

data(18).frame_type  = 'Sidesway_Uninhibited';
data(18).L           = 80*r;
data(18).kqtop       = Inf;
data(18).kqbot       = Inf;
data(18).gamma       = 0;
data(18).frame_name = sprintf('SR_%s80_G0',section_name(1)); 

data(19).frame_type  = 'Sidesway_Uninhibited';
data(19).L           = 80*r;
data(19).kqtop       = (6*EI)/(3*data(19).L);
data(19).kqbot       = (6*EI)/(3*data(19).L);
data(19).gamma       = 0;
data(19).frame_name = sprintf('SR_%s80_G3',section_name(1)); 

data(20).frame_type  = 'Sidesway_Inhibited';
data(20).L           = 80*r;
data(20).beta        = 1.0;
data(20).frame_name = sprintf('SCB_%s80',section_name(1));

data(21).frame_type  = 'Sidesway_Inhibited';
data(21).L           = 120*r;
data(21).beta        = 1.0;
data(21).frame_name = sprintf('SCB_%s120',section_name(1)); 

data(22).frame_type  = 'Sidesway_Inhibited';
data(22).L           = 80*r;
data(22).beta        = -0.5;
data(22).frame_name = sprintf('DCB_%s80',section_name(1)); 

data(23).frame_type  = 'Sidesway_Inhibited';
data(23).L           = 120*r;
data(23).beta        = -0.5;
data(23).frame_name = sprintf('DCB_%s120',section_name(1));

for i = 1:length(data)
    data(i).frame_id = i;
    data(i).Delta0 = data(i).L/500;
    data(i).delta0 = data(i).L/1000;
end

if strcmpi(axis,'weak')
    data = data([1:4 10:end]);
end

end