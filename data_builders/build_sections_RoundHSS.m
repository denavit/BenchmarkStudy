function sections = build_sections_RoundHSS()

iSection = 1;

Fy = [35 42 50 65];
Fu = [60 58 65 80];
sectionLetters = 'OABCDE';
chosenSections = [
     5.000 0.500
     7.000 0.500
    10.000 0.500
    12.750 0.375
    16.000 0.250
    24.000 0.125];
D = chosenSections(:,1);
t = 0.93*chosenSections(:,2);

for i_shape = 1:length(D)
    for i_Fy = 1:length(Fy)
        sections(iSection).axis = 'strong';
        sections(iSection).section = ...
            RoundHSS(D(i_shape),t(i_shape),Fy(i_Fy),'US');
        sections(iSection).section.Fu = Fu(i_Fy);
        sections(iSection).section.neglectLocalBuckling = true;
        sections(iSection).section_type = 'RoundHSS';
        sections(iSection).section_name = ...
            sprintf('RoundHSS-%s-%i',sectionLetters(i_shape),Fy(i_Fy));
        iSection = iSection+1;
    end
end

end