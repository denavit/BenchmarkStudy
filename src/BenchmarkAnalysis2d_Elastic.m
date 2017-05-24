function obj = BenchmarkAnalysis2d_Elastic(data,EI)
if nargin > 1
    if isnumeric(EI)
        data.EI = EI;
    elseif ischar(EI)
        data.EI = data.section.EI(data.axis,EI);
    else
        error('Unknown EI type: %s',type(EI));
    end
end
switch data.frame_type
    case 'Sidesway_Uninhibited'
        obj = BenchmarkAnalysis2d_Elastic_Sidesway_Uninhibited(data);
    case 'Sidesway_Inhibited'
        obj = BenchmarkAnalysis2d_Elastic_Sidesway_Inhibited(data);
    otherwise
        error('Unknown frame type: %s',data.frame_type);
end
end