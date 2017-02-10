function obj = BenchmarkAnalysis2d_Elastic(data)
switch data.frame_type
    case 'Sidesway_Uninhibited'
        obj = BenchmarkAnalysis2d_Elastic_Sidesway_Uninhibited(data);
    case 'Sidesway_Inhibited'
        obj = BenchmarkAnalysis2d_Elastic_Sidesway_Inhibited(data);
    otherwise
        error('Unknown frame type: %s',data.frame_type);
end
end