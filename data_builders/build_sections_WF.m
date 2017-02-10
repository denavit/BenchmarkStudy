function sections = build_sections_WF(axis)

iSection = 1;

Fy = [36 50 65];
Fu = [58 65 80];
shapeNames = {'W8x31','W14x43','W14x120','W14x311','W14x730',...
    'W18x311','W21x62','W40x264','W40x392','W40x593'};

for i_Fy = 1:length(Fy)
    for i_shape = 1:length(shapeNames)
        shp = steel_shape_lookup(shapeNames{i_shape});
        
        sections(iSection).axis = axis;
        sections(iSection).section = WF(...
            shp.d,shp.tw,shp.bf,shp.tf,shp.kdes,Fy(i_Fy),'US');
        sections(iSection).section.Fu = Fu(i_Fy);
        sections(iSection).section.eu = 0.20;
        sections(iSection).section.neglectLocalBuckling = true;
        sections(iSection).section.neglectLateralTorsionalBuckling = true;
        sections(iSection).section_type = 'WF';
        sections(iSection).section_name = ...
            sprintf('WF%s-%s-%i',lower(axis(1)),upper(listLetter(i_shape)),Fy(i_Fy));
        iSection = iSection+1;
    end
end

end