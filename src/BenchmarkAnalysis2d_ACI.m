function obj = BenchmarkAnalysis2d_ACI(data,EIeff_type)
switch data.frame_type
    case {'Sidesway_Uninhibited','Sway_Frame'}
        obj = BenchmarkAnalysis2d_ACI_Sway_Frame(data);
    case {'Sidesway_Inhibited','Nonsway_Frame'}
        obj = BenchmarkAnalysis2d_ACI_Nonsway_Frame(data);
    otherwise
        error('Unknown frame type: %s',data.frame_type);
end
if nargin > 1
    obj.EIeff_type = EIeff_type;
end
end