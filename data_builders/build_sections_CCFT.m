function sections = build_sections_CCFT()

iSection = 1;
fc = [4 8 16];

Fy = 42;
Fu = 58;
chosenSections = [
    7.000 0.500
    10.000 0.500
    12.750 0.375
    16.000 0.250
    24.000 0.125];
D = chosenSections(:,1);
t = 0.93*chosenSections(:,2);

sections(length(D)*length(fc)) = struct;
for i = 1:length(D)
    for j = 1:length(fc)
        sections(iSection).axis = 'strong';
        sections(iSection).section = CCFT(D(i),t(i),Fy,fc(j),'US');
        sections(iSection).section.Fu = Fu;
        sections(iSection).section.neglectLocalBuckling = true;
        sections(iSection).section.Ec = 57*sqrt(1000*fc(j));
        sections(iSection).section_type = 'CCFT';
        sections(iSection).section_name = ...
            sprintf('CCFT-%s-%i',upper(listLetter(i)),fc(j));
        iSection = iSection+1;
    end
end

end