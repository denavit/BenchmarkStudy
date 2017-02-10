function data = build_data_from_sections(sections)


frames = build_frames;
numSections = length(sections);
numFrames   = length(frames);
numData     = numSections*numFrames;

data(numData) = struct;
for iSection = 1:numSections
    axis = sections(iSection).axis;
    EIg  = sections(iSection).section.EI(axis,'Gross');
    Pnog = sections(iSection).section.getSectionData('GrossSectionCompressionStrength');
    
    for iFrame = 1:numFrames
        L = frames(iFrame).lambdaoe1g*pi*sqrt(EIg/Pnog);
        
        iData = iFrame + (iSection-1)*numFrames;
        
        section_fields = fields(sections(iSection));
        for i = 1:length(section_fields)
            data(iData).(section_fields{i}) = sections(iSection).(section_fields{i});
        end
        data(iData).section_id  = iSection;
        data(iData).frame_type  = frames(iFrame).frame_type;
        data(iData).frame_name  = frames(iFrame).frame_name;
        data(iData).frame_id    = iFrame;
        data(iData).L           = L;
        data(iData).Delta0      = L/500;
        data(iData).delta0      = L/1000;
        
        switch data(iData).frame_type
            case 'Sidesway_Uninhibited'
                data(iData).kqtop = (6*EIg)/(frames(iFrame).Ggtop*L);
                data(iData).kqbot = (6*EIg)/(frames(iFrame).Ggbot*L);
                data(iData).gamma = frames(iFrame).gamma;
            case 'Sidesway_Inhibited'
                data(iData).beta  = frames(iFrame).beta;
            otherwise
                error('Unknown frame type: %s',data(iData).frame_type)
        end
    end
    
end

end


function frames = build_frames
iFrame = 1;

%% Unbraced Frames
% Parameters
lambdaoe1g = [0.22 0.45 0.67 0.90];
gamma = [0 1 2 3];
endRestraint0 = [
    0 0
    3 3
    0 Inf
    3 Inf ];
endRestraint = [
    0 0
    1 1
    0 Inf
    1 Inf ];
% Define Frames
for i = 1:length(lambdaoe1g)
    for j = 1:length(gamma)
        for k = 1:length(endRestraint)
            frames(iFrame).frame_type = 'Sidesway_Uninhibited';
            frames(iFrame).lambdaoe1g = lambdaoe1g(i);
            frames(iFrame).gamma      = gamma(j);
            if ( frames(iFrame).gamma == 0 )
                frames(iFrame).Ggtop = endRestraint0(k,1);
                frames(iFrame).Ggbot = endRestraint0(k,2);
            else
                frames(iFrame).Ggtop = endRestraint(k,1);
                frames(iFrame).Ggbot = endRestraint(k,2);
            end
            frames(iFrame).frame_name = sprintf('U%s-%i-g%i',...
                upper(listLetter(k)),...
                100*frames(iFrame).lambdaoe1g,...
                frames(iFrame).gamma);
            iFrame = iFrame+1;
        end
    end
end

%% Braced Frames
% Parameters
lambdaoe1g = [0.45 0.90 1.35 1.90];
beta = [-0.5 0.0 0.5 1.0];
% Define Frames
for i = 1:length(lambdaoe1g)
    for j = 1:length(beta)
        frames(iFrame).frame_type = 'Sidesway_Inhibited';        
        frames(iFrame).lambdaoe1g = lambdaoe1g(i);
        frames(iFrame).beta       = beta(j);
        frames(iFrame).frame_name = sprintf('I-%i-b%i',...
            100*frames(iFrame).lambdaoe1g,j);
        iFrame = iFrame+1;
    end
end

end