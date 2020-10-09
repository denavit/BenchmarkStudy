function runBenchmarkStudy(study_path,action,option,delete_previous,save_figures)

if nargin == 0
    runBenchmarkStudyGUI
    return
end

if strcmp(action,'--------------------------')
    return
end

tic

%% Open Study
if exist(study_path,'dir') ~= 7
    study_path = fullfile(pathOf.BenchmarkStudyResults,study_path);
end
study = BenchmarkStudy('open',study_path);

study.scratch_path = fullfile(pathOf.scratch,study.study_name);
if exist(study.scratch_path,'dir') ~= 7
    mkdir(study.scratch_path)
end

%% Run Action
switch action 
    case 'List Results Files'
        study.list_results_files();
        
    case 'Analysis Interaction'
        tag = 'Analysis Interaction';
        
        fiber_section_definition_options = study.get_fiber_section_definition_options('Strength');
        analysis_options = study.get_analysis_options('Strength');
       
        if delete_previous && ~isempty(study.check_results_tag(tag))
            study.remove_results_file(tag);
        end
        runBenchmarkStudy_AnalysisInteraction(study,tag,fiber_section_definition_options,analysis_options);      

    case 'Figure - Analysis Interaction'
        tag = 'Analysis Interaction';
        runBenchmarkStudy_Figures_Interaction(study,tag,{'FN'},save_figures);        
        
    case 'Analysis Interaction (Section)'
        tag = 'Analysis Interaction (Section)';
        
        fiber_section_definition_options = study.get_fiber_section_definition_options('Strength');
        analysis_options = study.get_analysis_options('Strength');
        
        if delete_previous && ~isempty(study.check_results_tag(tag))
            study.remove_results_file(tag);
        end
        runBenchmarkStudy_AnalysisInteraction_Section(study,tag,fiber_section_definition_options,analysis_options);
               
    case 'Figure - Analysis Interaction (Section)'
        tag = 'Analysis Interaction (Section)';
        runBenchmarkStudy_Figures_Interaction(study,tag,{'FN'},save_figures);     
        
    case 'Design Interaction'
        tag    = sprintf('Design Interaction (%s)',option);
        tag_FN = 'Analysis Interaction';
        
        if delete_previous && ~isempty(study.check_results_tag(tag))
            study.remove_results_file(tag);
        end
        
        runBenchmarkStudy_DesignInteraction(study,tag,tag_FN,option);
        
    case 'Figure - Design Interaction'        
        tag    = sprintf('Design Interaction (%s)',option);
        tag_FN = 'Analysis Interaction';
        runBenchmarkStudy_Figures_Interaction(study,{tag,tag_FN},{'Design','FN'},save_figures);   
        
    case 'ACI Interaction'
        
        switch option 
            case 'a'
                EIeff_type = 'a';
                include_strength_reduction  = false;
                include_stiffness_reduction = false;
                second_order_moment_ratio_limit = 1.4;
                
            case 'b'
                EIeff_type = 'b';
                include_strength_reduction  = false;
                include_stiffness_reduction = false;
                second_order_moment_ratio_limit = 1.4;
                
            case 'c'
                EIeff_type = 'c';
                include_strength_reduction  = false;
                include_stiffness_reduction = false;
                second_order_moment_ratio_limit = 1.4;
                
            case 'a - no limit'
                EIeff_type = 'a';
                include_strength_reduction  = false;
                include_stiffness_reduction = false;
                second_order_moment_ratio_limit = Inf;
                
            case 'b - no limit'
                EIeff_type = 'b';
                include_strength_reduction  = false;
                include_stiffness_reduction = false;
                second_order_moment_ratio_limit = Inf;
                
            case 'c - no limit'
                EIeff_type = 'c';
                include_strength_reduction  = false;
                include_stiffness_reduction = false;
                second_order_moment_ratio_limit = Inf;
                
            otherwise
                error('Unknown option: %s',option);
        end
             
        tag = sprintf('ACI Interaction (%s)',option);
        if delete_previous && ~isempty(study.check_results_tag(tag))
            study.remove_results_file(tag);
        end
        runBenchmarkStudy_ACI(study,tag,EIeff_type,include_strength_reduction,include_stiffness_reduction,second_order_moment_ratio_limit)
       
    case 'ACI Interaction (Section)'
        
        tag = sprintf('ACI Interaction (Section) (%s)',option);
        if delete_previous && ~isempty(study.check_results_tag(tag))
            study.remove_results_file(tag);
        end
        runBenchmarkStudy_ACI_Section(study,tag,option)
        
    case 'Figure - ACI Interaction'
        
        section_tag = 'ACI Interaction (Section) (nominal)';
        if ~isempty(study.check_results_tag(section_tag))
            tag  = {section_tag};
            name = {'section'};
        else 
            tag  = cell(0,1);
            name = cell(0,1);
        end
        
        switch option 
            case {'a','b','c','a - no limit','b - no limit','c - no limit'}
                tag  = vertcat(tag,{sprintf('ACI Interaction (%s)',option)});
                name = vertcat(name,{option});
            case 'a & a - no limit'
                tag  = vertcat(tag,{'ACI Interaction (a)';'ACI Interaction (a - no limit)'});
                name = vertcat(name,{'a';'a - no limit'});
            case 'b & b - no limit'
                tag  = vertcat(tag,{'ACI Interaction (b)';'ACI Interaction (b - no limit)'});
                name = vertcat(name,{'b';'b - no limit'});
            case 'c & c - no limit'
                tag  = vertcat(tag,{'ACI Interaction (c)';'ACI Interaction (c - no limit)'});
                name = vertcat(name,{'c';'c - no limit'});
            otherwise
                error('Unknown option: %s',option);
        end
        
        runBenchmarkStudy_Figures_Interaction(study,tag,name,save_figures);

    case 'Figure - ACI Interaction (Section)'
        tag = sprintf('ACI Interaction (%s)',option);
        runBenchmarkStudy_Figures_Interaction(study,tag,{'ACI'},save_figures);        
        
    otherwise
        error('Unknown action: %s',action)
end
   
t = toc;
if t >= 60
    fprintf('\nElapsed time is %s\n',datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'));
end

end
