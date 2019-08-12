function runBenchmarkStudy_Figures_Interaction(study,tags,names,saveFigures)

selectedData = [];
if iscell(study)
    selectedData = study{2};
    study = study{1};
end

if ischar(tags)
    tags = {tags};
end

if ischar(names)
    names = {names};
end

if ischar(saveFigures)
    folderName = saveFigures;
    saveFigures = true;
end

%% Load Data
% Load Sections
load(study.path_of_data);
numData = length(data);
if isempty(selectedData)
    selectedData = 1:numData;
end

% Load Results
numResults = length(tags);
results_combined = cell(numResults,1);
for i = 1:numResults
    load(study.path_of_results(tags{i}));
    results_combined{i} = results;
    clear results
end

myColors = lines(numResults);


%% Interaction Diagrams
fs = figureStyle('Display');

h1 = NaN;
for iData = selectedData
    if ishandle(h1)
        close(h1);
    end

    h1 = fs.figure(6,4);
    ha = fs.axes([0.12 0.20 0.85 0.77]);

    legendNames   = [];
    legendHandles = cell(1,0);

    for i = 1:numResults
        h = plotInteraction(results_combined{i}(iData),myColors(i,:));
        legendHandles = horzcat(legendHandles,h);
        legendNames   = horzcat(legendNames,names(i));
    end

    % Set Axis Limtis
    axis_limits('margin',[0.35 0.20])
    axis_limits('quadrant',1)
    
    % Lables
    section_name = strrep(data(iData).section_name,'_','\_');
    if isfield(data(iData),'frame_name')
        frame_name   = strrep(data(iData).frame_name,'_','\_');
    else
        frame_name = '';
    end
    
    extraInfo = sprintf('Case %i',iData');
    %extraInfo = sprintf('Case %i - Section %i: %s - Frame %i: %s',...
    %    iData,...
    %    data(iData).section_id,section_name,...
    %    data(iData).frame_id,frame_name);
    xlabel({'Bending Moment (M)',extraInfo})
    ylabel('Axial Compression (P)')
    legend(legendHandles,legendNames,'Location','NE');

    if saveFigures
        if exist(folderName,'dir') ~= 7
            mkdir(folderName);
        end
        figureName1 = sprintf('figure%02i',iData);
        figureTools.print(h1,fullfile(folderName,[figureName1 '.fig']),'fig')
        %figureTools.print(h1,fullfile(folderName,[figureName1 '.eps']),'eps')
        figureTools.print(h1,fullfile(folderName,[figureName1 '.png']),'png',100)
        close(h1);

        if iSection == selectedSections(end) && iFrame == selectedFrames(end)
            options = struct;
            options.openWebpage = false;
            figureWebpage(folderName,'Benchmark Study Interaction',options)
        end
    else
        if iData ~= selectedData(end)
            pause
        end
    end
end

end


function h = plotInteraction(results,color)

switch results.study_type
    case 'AnalysisInteraction'
        h = plot(results.M1,-results.P1,'-','LineWidth',2,'Color',color);
            plot(results.M2,-results.P2,'--','LineWidth',2,'Color',color);
        
        for i = 1:length(results.P1)
            marker_type = get_marker_type(results.limit_type{i});
            plot(results.M1(i),-results.P1(i),marker_type,'LineWidth',2,'Color',color,'MarkerSize',8);
            plot(results.M2(i),-results.P2(i),marker_type,'LineWidth',2,'Color',color,'MarkerSize',8);
        end

    case 'AnalysisInteraction_Section'
        h = plot(results.M1,-results.P1,'-','LineWidth',2,'Color',color);
        
        for i = 1:length(results.P1)
            marker_type = get_marker_type(results.limit_type{i});
            plot(results.M1(i),-results.P1(i),marker_type,'LineWidth',2,'Color',color,'MarkerSize',8);
        end        
        
    case 'DesignInteraction'
        h = plot(results.M1,-results.P1,'-','LineWidth',2,'Color',color);
            plot(results.M2,-results.P2,'--','LineWidth',2,'Color',color);
            plot(results.M2e,-results.P2e,'-.','LineWidth',2,'Color',color);        
        
    case 'ACI Interaction'
        h = plot(results.M1,-results.P1,'-','LineWidth',2,'Color',color);
            plot(results.M2,-results.P2,'--','LineWidth',2,'Color',color);
            
    case 'ACI Interaction (Section)'
        h = plot(results.M1,-results.P1,'-','LineWidth',2,'Color',color);
            
    otherwise
        error('Unknown study type: %s',results.study_type)
end

end


function marker_type = get_marker_type(limit_type)

if isempty(limit_type)
    fprintf('Unknown limit type: empty\n');
    marker_type = 'o';
	return
end
    
switch limit_type
    case 'Stability Limit (Reached)'
        marker_type = 's';
    case 'Stability Limit (Tolerance Reached)'
        marker_type = 'd';
    case 'Strain Limit (Reached)'
        marker_type = '^';
    case 'Concrete Compressive Strain Limit (Reached)'
        marker_type = 'v';
    case 'Compressive Force Limit (Reached)'
        marker_type = '<';
    otherwise
        fprintf('Unknown limit type: %s\n',limit_type);
        marker_type = 'o';
end

end
