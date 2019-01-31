clear all; close all; clc;

%% Input
type = 'RoundHSS';
%type = 'RectangularHSS';

E  = 29000;
L_over_r = 20:20:200;
units = 'US';
switch type
    case 'RoundHSS'
        Fy = 42;
        shapes = {'HSS5.000x0.500','HSS7.000x0.500','HSS10.000x0.500','HSS12.750x0.375','HSS16.000x0.250'};
    case 'RectangularHSS'
        Fy = 46;
        shapes = {'HSS6x6x1/2','HSS8x8x1/2','HSS9x9x1/2','HSS8x8x1/4','HSS9x9x1/8'};
    otherwise
        error('Unknown section type: %s',type)
end


%% Run Analyses 
numLengths = length(L_over_r);
numShapes  = length(shapes);

P_analysis = nan(numShapes,numLengths);
Py = nan(numShapes,1);

for iShape = 1:numShapes
    shape_data = steel_shape_lookup(shapes{iShape});
    Py(iShape) = shape_data.A*Fy;
    
    switch type
        case 'RoundHSS'
            section = RoundHSS(shape_data.OD,shape_data.tdes,Fy,units);
            r = shape_data.rx;
        case 'RectangularHSS'
            section = RectangularHSS(shape_data.Ht,shape_data.B,shape_data.tdes,Fy,units);
            r = shape_data.rx;
        otherwise
            error('Unknown section type: %s',type)
    end
    
    for iLength = 1:numLengths
        data = struct;
        data.L          = L_over_r(iLength)*r;
        data.axis       = 'Strong';
        data.section    = section;
        data.type       = type;
        data.frame_type = 'Sidesway_Inhibited';
        data.beta       = 1;
        data.delta0     = data.L/1000;
        
        fsDefOptions = struct;
        fsDefOptions.nf1 = 30;
        fsDefOptions.SteelMaterialType = 'AbdelRahman_LowHardening';
        fsDefOptions.StrengthReduction = 0.9;
        fsDefOptions.StiffnessReduction = 0.9;
        fsDefOptions.StrainHardeningRatio = 1/1000;
        fsDefOptions.includePackageDefinition = true;

        [sectionDefinition,numMat] = FiberSectionDefinition(data.section,data.axis,1,1,fsDefOptions);

        analysisOptions = struct;
        analysisOptions.absoluteStrainLimit = 0.05;

        ba = BenchmarkAnalysis2d_OpenSees(data,sectionDefinition,analysisOptions);
        
        iResults = ba.runAnalysis('LimitPoint_Proportional',0,[],1);
        
        if ~iResults.limitPoint.good
            iResults = ba.runAnalysis('LimitPoint_Proportional',0,[],2);
        end
        
        if ~iResults.limitPoint.good
            error('Limit Point Not Obtained: Axial Only');
        end
        
        P_analysis(iShape,iLength) = -iResults.limitPoint.P1;
    end
end

%% Plots
fs = figureStyle('Display');
myColors = lines(numShapes);

% Plot 1
hf = fs.figure(8,6);
ha = fs.axes;

legend_h = nan(1,numShapes+1);
legend_n = cell(1,numShapes+1);

L_over_r_plot = linspace(0,200,200);
P_AISC = AISC_column_curve(Fy/(pi^2*E)*L_over_r_plot.^2);
legend_h(1) = plot(L_over_r_plot,0.9*P_AISC,'k-','LineWidth',2);
legend_n{1} = '\phiP_n';

for iShape = 1:numShapes
    legend_h(iShape+1) = plot(L_over_r,P_analysis(iShape,:)/Py(iShape),'o-','LineWidth',2,'Color',myColors(iShape,:));
    legend_n{iShape+1} = shapes{iShape};
end

xlim([0 200])
ylim([0 1.1])
xlabel('Member Slenderness (L/r)')
ylabel('Normalized Axial Compression (P/A_gF_y)')
legend(legend_h,legend_n,'Location','NE')

% Plot 2
hf = fs.figure(8,6);
ha = fs.axes;

legend_h = nan(1,numShapes+1);
legend_n = cell(1,numShapes+1);

P_AISC = AISC_column_curve(Fy/(pi^2*E)*L_over_r.^2);
legend_h(1) = plot([0 200],[1 1],'k-','LineWidth',2);
legend_n{1} = '\phiP_n';

for iShape = 1:numShapes
    legend_h(iShape+1) = plot(L_over_r,P_analysis(iShape,:)./(0.9*P_AISC*Py(iShape)),'o-','LineWidth',2,'Color',myColors(iShape,:));
    legend_n{iShape+1} = shapes{iShape};
end

xlim([0 200])
ylim([0.9 1.1])
xlabel('Member Slenderness (L/r)')
ylabel('Normalized Axial Compression (P/\phiP_n)')
legend(legend_h,legend_n,'Location','NW')
