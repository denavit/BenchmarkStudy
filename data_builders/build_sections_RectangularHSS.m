function sections = build_sections_RectangularHSS()

iSection = 1;

Fy = [36 46 50 65];
Fu = [58 58 65 80];
sectionLetters = 'OABCDE';
chosenSections = [
    6.0 6.0 1/2
    9.0 9.0 1/2
    8.0 8.0 1/4
    9.0 9.0 1/8
    14.0 14.0 1/8];
D = chosenSections(:,1);
B = chosenSections(:,2);
t = 0.93*chosenSections(:,3);

for i_shape = 1:length(D)
    for i_Fy = 1:length(Fy)
        sections(iSection).axis = 'strong';
        sections(iSection).section = ...
            RectangularHSS(D(i_shape),B(i_shape),t(i_shape),Fy(i_Fy),'US');
        sections(iSection).section.Fu = Fu(i_Fy);
        sections(iSection).section.neglectLocalBuckling = true;
        sections(iSection).section_type = 'RectangularHSS';
        sections(iSection).section_name = ...
            sprintf('RectHSS-%s-%i',sectionLetters(i_shape),Fy(i_Fy));
        iSection = iSection+1;
    end
end

end