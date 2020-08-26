clear all; close all; clc;

study_name = 'RC';
study_path = fullfile(pathOf.BenchmarkStudyResults,study_name);

create_study    = false;
build_data      = true;
build_options   = true;

%% Create Study
if create_study
    study = BenchmarkStudy('create',study_path,study_name);
end

%% Build Data
if build_data
    study = BenchmarkStudy('open',study_path);
    data = build_data_RC();
    save(study.path_of_data(),'data')
end

%% Build Options
if build_options
    study = BenchmarkStudy('open',study_path);

    % Define fiber_section_definition_options
    fiber_section_definition_options = struct;
    fiber_section_definition_options.nf1 = 30;
    fiber_section_definition_options.includePackageDefinition = true;
    fiber_section_definition_options.conc_material = 'Concrete04';
    fiber_section_definition_options.steel_material = 'ElasticPP';
    study.set_fiber_section_definition_options('Strength',fiber_section_definition_options)
    
    % Define analysis_options
    analysis_options = struct;    
    analysis_options.absoluteStrainLimit = 0.05;
    analysis_options.store_extra_data = false;
    study.set_analysis_options('Strength',analysis_options)
end
    