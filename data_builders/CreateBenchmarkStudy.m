clear all; close all; clc;

study_name = 'xxxx';
study_path = fullfile(pathOf.BenchmarkStudyResults,study_name);

create_study    = false;
build_data      = false;
build_options   = false;

%% Create Study
if create_study
    study = BenchmarkStudy('create',study_path,study_name);
end

%% Build Data
if build_data
    study = BenchmarkStudy('open',study_path);
    
    switch study_name
        case 'CCFT'
            data = build_data_from_sections(build_sections_CCFT());
        case 'RCFT'
            data = build_data_from_sections(build_sections_RCFT());
        case 'SRCs'
            data = build_data_from_sections(build_sections_SRC('strong'));
        case 'SRCw'
            data = build_data_from_sections(build_sections_SRC('weak'));
        otherwise
            error('Unknown study name: %s',study_name)
    end

    save(study.path_of_data(),'data')
end

%% Build Options
if build_options
    study = BenchmarkStudy('open',study_path);

    % Define options structures
    fiber_section_definition_options = struct;
    analysis_options = struct;

    fiber_section_definition_options.nf1 = 30;
    fiber_section_definition_options.includePackageDefinition = true;
    analysis_options.absoluteStrainLimit = 0.05;

    switch study_name
        case {'CCFT','RCFT'}
            fiber_section_definition_options.SteelMaterialType = 'AbdelRahman';
            fiber_section_definition_options.ConcreteMaterialType = 'ProposedForDesign';
        case {'SRCs','SRCw'}
            fiber_section_definition_options.SteelMaterialType = 'ElasticSmallStiffness';
            fiber_section_definition_options.ConcreteMaterialType = 'ProposedForDesign';
            fiber_section_definition_options.ReinforcementMaterialType = 'ElasticSmallStiffness';
        otherwise
            error('Unknown study name: %s',study_name)
    end

    study.set_fiber_section_definition_options('Strength',fiber_section_definition_options)
    study.set_analysis_options('Strength',analysis_options)
end
    