classdef BenchmarkStudy < handle
    
    properties
        study_name
        study_path
        results_filenames
        results_tags
        scratch_path = ''        
        numPoints = 9;
        
        fiber_section_definition_options_name = cell(0,1);
        fiber_section_definition_options_S = cell(0,1);
        analysis_options_name = cell(0,1);
        analysis_options_S = cell(0,1);
    end
    
    methods 
        
        %% Constructor
        function obj = BenchmarkStudy(action,study_path,study_name)            
            switch action
                case 'create'                       
                    assert(exist(study_path,'dir')~=7,...
                        'Cannot create study, path already exists');
                    
                    % Set Properties
                    obj.study_name          = study_name;
                    obj.study_path          = study_path;
                    obj.results_filenames   = cell(0,1);
                    obj.results_tags        = cell(0,1);
                    
                    % Create Directories
                    mkdir(obj.study_path);
                    
                    % Save study to mat file
                    obj.save   
                    
                case 'open'
                    study_file = fullfile(study_path,'study.mat');
                    
                    assert(exist(study_file,'file')==2,...
                        'Cannot open study, study.mat does not exit')
                    
                    load(study_file)
                    obj.study_path = study_path;
                    
                otherwise
                    error('Invalid action');
            end
        end

        %% Save function
        function save(obj)
            save(obj.path_of_study,'obj');
        end        
        
        %% Functions that modify the study
        function path = addResultsFile(obj,tag)
            % Check tag
            assert(isempty(obj.check_results_tag(tag)),...
                'tag should be unique');
            
            % Create file path
            filename = sprintf('%s.mat',datestr(now,'yyyy-mm-dd-HH-MM-SS'));
            
            % Update study properties
            obj.results_filenames = vertcat(obj.results_filenames,filename);
            obj.results_tags  = vertcat(obj.results_tags,{tag});
            
            % Save study to mat file
            obj.save;
            
            % Create file path
            path = fullfile(obj.study_path,filename);
        end
        function removeResultsFile(obj,tag)            
            % Check tag
            ind = obj.check_results_tag(tag);
            
            % Return if tag does not exist
            if isempty(ind)
                return
            end
            
            % Delete file
            delete(fullfile(obj.study_path,obj.results_filenames{ind}));
            
            % Update study properties
            obj.results_filenames = vertcat(...
                obj.results_filenames(1:ind-1,1),obj.results_filenames(ind+1:end,1));
            obj.results_tags = vertcat(...
                obj.results_tags(1:ind-1,1),obj.results_tags(ind+1:end,1));
            
            % Save study to mat file
            obj.save;
        end
        
        %% Fuctions that provide information on the study
        function ind = check_results_tag(obj,tag)
            ind = find(strcmp(obj.results_tags,tag)==1);           
            assert(numel(ind)<2,'Tag is not unique');
        end
        function path = path_of_study(obj)
            path = fullfile(obj.study_path,'study.mat');
        end
        function path = path_of_data(obj)
            path = fullfile(obj.study_path,'data.mat');
        end
        function path = path_of_results(obj,tag)
            ind = obj.check_results_tag(tag);
            assert(~isempty(ind),'Tag does not exist');
            path = fullfile(obj.study_path,obj.results_filenames{ind});
        end
        function list_results_files(obj)
            numFiles = length(obj.results_filenames);
            
            if numFiles == 0
                fprintf('No results files to list\n');
            else
                fprintf('Results Files:\n');
                for i = 1:numFiles
                    fprintf('  Tag: %s, Path: %s\n',...
                        obj.results_tags{i},...
                        obj.results_filenames{i});
                end
            end
        end
        
        % 
        function set_fiber_section_definition_options(obj,name,S)
            match = strcmp(obj.fiber_section_definition_options_name,name);
            if sum(match) == 0
                obj.fiber_section_definition_options_name = vertcat(...
                    obj.fiber_section_definition_options_name,name);
                obj.fiber_section_definition_options_S = vertcat(...
                    obj.fiber_section_definition_options_S,S);
            elseif sum(match) == 1
                obj.fiber_section_definition_options_S{match == 1} = S;
            else
                error('name exists in study twice or more')
            end
            obj.save;
        end
        function S = get_fiber_section_definition_options(obj,name)
            match = strcmp(obj.fiber_section_definition_options_name,name);
            if sum(match) == 0
                error('fiber section definition options "%s" does not exist',name);
            elseif sum(match) == 1
                S = obj.fiber_section_definition_options_S{match == 1};
            else
                error('name exists in study twice or more')
            end
        end
        function set_analysis_options(obj,name,S)
            match = strcmp(obj.analysis_options_name,name);
            if sum(match) == 0
                obj.analysis_options_name = vertcat(obj.analysis_options_name,name);
                obj.analysis_options_S = vertcat(obj.analysis_options_S,S);
            elseif sum(match) == 1
                obj.analysis_options_S{match == 1} = S;
            else
                error('name exists in study twice or more')
            end
            obj.save;            
        end
        function S = get_analysis_options(obj,name)
            match = strcmp(obj.analysis_options_name,name);
            if sum(match) == 0
                error('analysis options "%s" does not exist',name);
            elseif sum(match) == 1
                S = obj.analysis_options_S{match == 1};
            else
                error('name exists in study twice or more')
            end            
        end
    end
end

