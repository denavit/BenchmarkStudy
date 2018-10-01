function sections = build_sections_SRC(axis)

iSection = 1;

H = 28;
B = 28;
Fy = 50;
Fu = 65;
Fyr = 60;
Fur = 90;
shapeNames = {'W14x311','W14x233','W12x120','W8x31'};
reinfConfig = {'6x-6y','4x-4y','2x-2y'};
db = {'#11','#10','#8'};
dbTies = '#3';
s = 12;
fc = [4 8 16];
cover = 1.5;

sections(length(shapeNames)*length(reinfConfig)*length(fc)) = struct;
for j = 1:length(shapeNames)
    for k = 1:length(reinfConfig)
        for l = 1:length(fc)
            shp = steel_shape_lookup(shapeNames{j});
            
            sections(iSection).axis = axis;
            sections(iSection).section = SRC(...
                shp.d,shp.tw,shp.bf,shp.tf,Fy,H,B,fc(l),db{k},reinfConfig{k},Fyr,...
                dbTies,s,Fyr,cover,'US');
            sections(iSection).section.Fu = Fu;
            sections(iSection).section.Fulr = Fur;
            sections(iSection).section.Ec = 57*sqrt(1000*fc(l));
            sections(iSection).section_type = 'SRC';
            sections(iSection).section_name = sprintf('SRC%s-%s-%i',...
                lower(axis(1)),upper([listLetter(j) listLetter(k)]),fc(l));
            iSection = iSection+1;
        end
    end
    
end

end