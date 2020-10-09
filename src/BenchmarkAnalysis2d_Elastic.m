function obj = BenchmarkAnalysis2d_Elastic(data,EI)
if nargin > 1
    if isnumeric(EI)
        data.EI = EI;
    elseif ischar(EI)
        if strcmpi(EI,'RC_Study_With_Stiffness_Reduction')
            data.EI = 0.7*data.section.Ec*data.section.Ig(data.axis);
            if strcmpi(data.frame_type,'Sidesway_Uninhibited')
                data.kqtop = 6*0.35*data.EcIgb_over_Lb;
                data.kqbot = 6*0.35*data.EcIgb_over_Lb;
            end
        elseif strcmpi(EI,'RC_Study_Without_Stiffness_Reduction')
            data.EI = 0.8*data.section.Ec*data.section.Ig(data.axis);
            if strcmpi(data.frame_type,'Sidesway_Uninhibited')
                data.kqtop = 6*0.4*data.EcIgb_over_Lb;
                data.kqbot = 6*0.4*data.EcIgb_over_Lb;
            end
        else
            data.EI = data.section.EI(data.axis,EI);
        end
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