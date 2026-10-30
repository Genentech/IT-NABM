%% Calibrate the NOP model, run sensitivity analysis (PRCC/Sobol), and plot
% End-to-end calibration driver for the Neuro Oligo Platform (NOP) model:
% builds the model, loads the cleaned data and parameters, adds variants and
% doses, optionally runs a parameter-estimation calibration
% (args.switches.calibration), optionally computes confidence intervals and
% PRCC/Sobol global sensitivity analyses, saves the PSA results, and plots the
% group fits.
% Models run on the default hour time unit and the data is all transformed to
% have hour time units.

clc
close all
close(findall(0, 'type', 'figure')); % closes fit figures
clear all
warning('off','all')
% restoredefaultpath;

file_path = fileparts(fileparts(mfilename('fullpath'))); % project root (one up from /calibration)

% repo location
addpath(genpath(file_path)); 

% Read a task log for scenarios to run and simulate
% a sheet with Locations 
% a sheet with Tasks

args = loadArgs([file_path,'/args_NOP.xlsx']);
args.projLoc = file_path;

% args.taskLog = loadTaskLog('/Users/jafarnem/Documents/MATLAB/C1s/taskLog_C1s.xlsx');

% switches for specific simulations
args.switches.test = 1;
args.switches.calibration = 0; % 1: run particleswarm calibration, 0: skip

%% Build or load the model
% model = buildModel_avidity_v1(locations,switches);
% modelSbproj = sbioloadproject(args.modelLoc); % load the sbioproj file
% model = modelSbproj.m1; % Main model for the effort
model = buildModel_NOP(args); % Build model
% model = buildModel_NOP_delay(args);
% model = buildModel_NOP_two_indirect(args);
% model = buildModel_NOP_delay_2(args); % Build model
% model = buildModel_NOP_delay_CSF(args); % Build model
% model = buildModel_NOP_aggregate(args);

% TASK: Update the model with CNS PD
% species specific parameters and flagged translation

%% Load the data
% tall data for simulation/data overlay
data = loadClean_dataset_NOP(args);

% read the wide dataset that is similar to the grouped data used for model calibration
dataWide = readtable(args.dataWideLoc);

%% Load parameters
args.params = readtable(args.paramLoc);

%% Add Variant
model = addVariant_NOP(model,data,args);
% model = addVariant_C1s_baseline(model,data,args);

% takes name of the "Name" and "Value" columns to generate the variant
model = add_variantParamTable(model, args,'Name','Value'           ,Name='varFit');
% model = add_variantParamTable(model, args,'Name','GantTronti_calib',Name='GantTronti_calib');
% model = add_variantParamTable(model, args,'Name','Denali_calib'    ,Name='Denali_calib');
% model = add_variantParamTable(model, args,'Name','Genentech_calib' ,Name='Genentech_calib');

%% Add doses
model = addDose_NOP(model,data,args);

%% load the tasks
args = loadTasks(model,data,args);
args = addDose2taskTable(model,data,args);

%% Run the model for the tasks in log and plot the data if available


%% Calibrate the model
% Groups to be used for calibrated in args.groupFit
% MAPT ASO Cyno EastonA - [250805001]
% MAPT ASO Mouse EastonA - [250801001:250801004];
 % Biogen Slide - [250808001:250808004];
 % Internal BIIB080 - [250811001:250811004]
 % EastonA Mouse - [260602001]
 % Kordasiewicz Mouse - [260603001:260603003]
 % Monine - [250812001:250812004]
 % Mummery - [250813001:250813004]
 % BP41660PET - [250814001]
 % Devos - [250822001:250822003]
 % Edwards - [250829001:250829006]
 % Tominersen - [250901001]
 % BrownKM2021 - Cyno - [250902001:250902005] [Slide 30, Slide 29, Slide 31, Slide32]
 % BrownKM2023 - Cyno - [260526001]

 % BrownKM2021 - Human - [250904001:250904003]

 % NorrisDA - Mouse - [250905001]
 % WO_2011_032045_A1_Mouse - [260525001]

 % TabriziSJ - Cyno - [250909001:250909007]
 % HTT - Human - [250910001:250910010]

 % Jacifusen - Mouse - [250911001]
 % Jacifusen - Cyno - [250912001:250912002]
 % Jacifusen - Human - [250913001:250913004]

args.groupFit = [260525001];

args.groupFit = groupOrder(args,data);


if args.switches.calibration
% Specify the calibration settings for the algorithm
calibrationSettings.Algorithm    = ['particleswarm']; % 'particleswarm','patternsearch','ga','lsqnonlin','scattersearch'
responseMapTable = readtable(args.responseMapLoc);
include = find(responseMapTable.Include); 
calibrationSettings.ErrorModel   = repmat({'constant'},1,numel(include)); % 'constant','proportional','combined','exponential'
calibrationSettings.Weights      = []; % Can be function handle weightsFunction = @(y) 1./y.^2;
calibrationSettings.UseParallel  = true;
calibrationSettings.Pooled       = 'auto'; % MUST use 'auto' for categoryVariable based estimates. if specify false or true it overrides the categoryVariable
% calibrationSettings.Pooled       = true;

calibrationSettings.options.MaxIterations      = 'auto'; % 'auto', 500, 1000
calibrationSettings.options.MaxStallIterations = 50*2;

args.calibrationSettings = calibrationSettings;

% Run the calibration
% Follow these steps:
% 1. Select the relevant groups and make sure tasks are updated
% 2. Update the parameter table with the correct paramters of interest and p0, LB/UB etc
% 3. Make sure responseMap is updated with correct species mapping for fitting
% 4. Check to have the desired algorithm settings (e.g. Pooled, ErrorModel, etc) and setup
args = runCalibration(model, args);
% plot(args.calibration.results); % QC plot of the fit and data from simbio

% Update varFit to include the new calibration (TASK: update to work with category-specific fitting)
model = updateVarFit(model, args);

% reload the tasks to ensure that the updated variants are being used
loadTasks(model,data,args);
args = addDose2taskTable(model,data,args);

%% Capture the calibrated parameter set
% close the parameters file
system('osascript -e ''tell application "Microsoft Excel" to close workbook "parameters.xlsx" saving yes''')
% TASK: make it so that it overwrites the Values as well
saveCalibration(args,Name='Genentech_test',Overwrite=1); % Fit210151, GantTronti, ChugaiSC, 
% re-open the parameter excel file to ensure it is not unintentionally overwritten
system(append('open ',args.paramLoc));
end % args.switches.calibration

%% Confidence interval
args.switches.CI = 0;
if args.switches.CI
    % Parameter confidence interval
    ciParam = sbioparameterci(args.calibration.results);
    ci2table(ciParam)
    plot(ciParam)
    
    % In matlab documentation there are exmaples on comparing pooled and unpooled fits
    
    % Prediction confidence interval
    ciPred = sbiopredictionci(args.calibration.results);
    plot(ciPred)
end

%% PSA

% Run PRCC PSA
args.switches.PRCC = 1;
if args.switches.PRCC

    args.paramsPSA = readtable(args.paramPSALoc);
    args.PRCC.Group = [251301001]; % MAPT ASO human translation 60 mg/kg Q12W
    args.PRCC.outputNames = {'V_bsf.mRNA_perc','V_bsf.Prot_perc','V_bsf.Oligo_tot_homog_nM'};
    args.PRCC.outputNamesDisplay = {'mRNA %','Protein %','Brain Oligo (nM)'};

    args = runPRCC(model, args, outputTime=[24*840]); % final timepoint = 4 weeks (hours)
    

    % ------ plotting --------------------
    % Heatmap (Task: add a column to the parameter table for short descriptive name for the variables to be used as labels)
    f=figure; tiledlayout(2,1,TileSpacing="compact");
    nexttile; h = heatmap(args.PRCC.inputNamesDisplay, args.PRCC.outputNamesDisplay, args.PRCC.rho'); h.ColorLimits = [-1,1]; h.Title='PRCC'; h.Colormap=turbo;
    nexttile; h = heatmap(args.PRCC.inputNamesDisplay, args.PRCC.outputNamesDisplay, args.PRCC.pval'); h.Title='p-value'; h.ColorScaling="log"; h.Colormap=flipud(parula);h.ColorLimits = [-3,1];
    exportgraphics(f,fullfile(args.figLoc,['PRCC_heatmap',char(args.figFormat)]),'ContentType','vector','BackgroundColor','none');

    
    % ------  Horizontal Bar chart --------------------
    f=figure; 
    f.Position(4) = max(f.Position(4), 28*numel(args.PRCC.inputNames)*numel(args.PRCC.outputNames));
    h = barh(args.PRCC.rho); ax=gca; ax.XLabel.String ='PRCC'; l=legend; l.String=args.PRCC.outputNamesDisplay; ax.XLim=[-1,1]; ax.YTickLabel = args.PRCC.inputNamesDisplay;
    for i=1:numel(args.PRCC.outputNames)
        for j=1:numel(args.PRCC.inputNames)
            text(1.02,h(i).XEndPoints(j),num2str(args.PRCC.pval(j,i),'%.2e'));
        end
    end
    text(1.02,numel(args.PRCC.inputNames)+0.5,"p-value");
    exportgraphics(f,fullfile(args.figLoc,['PRCC_HorizontalBar',char(args.figFormat)]),'ContentType','vector','BackgroundColor','none');


    % ------  Tiled Bar chart for even more parameters --------------------
    % figure; tiledlayout
    % for i=1:numel(args.PRCC.outputNames)
    %     nexttile; h = barh(args.PRCC.rho(:,i)); h.XLabel="PRCC"; h.Title=args.PRCC.outputNames{i};h.yticklabels = args.PRCC.inputNames;
    % end


    % ------ Figure out a Sharknado alternative -------------------- 
    % plotSharknado(args.PRCC.results );


    % ------ Parameter distributions -------------------- 
    f=figure; [~,ax] = plotmatrix(args.PRCC.input);
    iterations = size(ax,1);
    for i = 1:iterations
        ax(i,1).YLabel.String = args.PRCC.inputNamesDisplay(i);
        ax(1,i).XLabel.String = args.PRCC.inputNamesDisplay(i);
        ax(1,i).XAxisLocation = "top";
        % ax(iterations,i).XLabel.String = args.PRCC.inputNames(i);
    end
    exportgraphics(f,fullfile(args.figLoc,['PRCC_ParamDist',char(args.figFormat)]),'ContentType','vector','BackgroundColor','none');

end

% Run Sobol PSA
args.switches.Sobol = 1;
if args.switches.Sobol

    args.paramsPSA = readtable(args.paramPSALoc);
    args.sobolGroup = [251301001]; % MAPT ASO human translation 60 mg/kg Q12W

    sobolOutputs = {'V_bsf.mRNA_perc','V_bsf.Prot_perc','V_bsf.Oligo_tot_homog_nM'};

    args = runSobol(model, args, sobolOutputs = sobolOutputs, outputTime=[0:1:24*168]);
    
    % plotting
    f=plotData(args.Sobol.results );
    exportgraphics(f,fullfile(args.figLoc,['Sobol_plotData',char(args.figFormat)]),'ContentType','vector','BackgroundColor','none');

    f=plot(args.Sobol.results );
    exportgraphics(f,fullfile(args.figLoc,['Sobol_plot',char(args.figFormat)]),'ContentType','vector','BackgroundColor','none');

    f=bar(args.Sobol.results );
    exportgraphics(f,fullfile(args.figLoc,['Sobol_HorizontalBar',char(args.figFormat)]),'ContentType','vector','BackgroundColor','none');

end

%% Save PSA results for standalone plotting
% Persist only the analysis results (not the whole workspace) so that
% runPlot_PSA_NOP() can regenerate figures without re-running the analysis.
if args.switches.PRCC || args.switches.Sobol
    psaResultsLoc = fullfile(file_path,'calibration','PSA_results_NOP.mat');
    PSA = struct();
    if args.switches.PRCC,  PSA.PRCC  = args.PRCC;  end
    if args.switches.Sobol, PSA.Sobol = args.Sobol; end
    PSA.figLoc    = args.figLoc;
    PSA.figSize   = args.figSize;
    PSA.figFormat = args.figFormat;
    PSA.fontSize  = args.fontSize;
    save(psaResultsLoc,'PSA','-v7.3');

    % Plot and save the PSA figures (PRCC and/or Sobol)
    runPlot_PSA_NOP(args);
end

%% GoF plots and tables
% AUC and Cmax comparisons


% Obs vs Pred



%% Run and plot the model diagnostics with a specified variant & dose
% It can be used in different ways
% 1. a group (prebuilt in taskTable) and run the model and plot (DONE)
% 2. TASK: a doseCellArray and a varaintCellArray and run the model and plot

% runPlot_diagnostics_C1s(model,data,args,type='group',group=210015102, stopTimeDay=30) % pC1s
% runPlot_diagnostics_C1s(model,data,args,type='group',group=210015103, stopTimeDay=30) % C1sRec

%% Run the model and plot the results per study 

% runPlot_C_21_0151(model,data,args)
% runPlot_C1s_C_all(model,data,args)
% runPlot_Grimm2023(model,data,args)
% runPlot_Moshkforoush2023(model,data,args)
% runPlot_C_gDpooled(model,data,args)
% args.taskTable.VariantFit(36,1) = args.taskTable.VariantFit(35,1);
% model.Parameters(92,1).Value = 1;
% args.taskTable.VariantOthers1(42,1).Content{14,1}{1,4} = 1E-4;
% args.taskTable.VariantOthers1(42,1).Content{15,1}{1,4} = 1E-3; %12.86;
% model.Compartments(4,1).Value = 12.86;
% args.taskTable.VariantOthers2(42,1).Content{10,1}{1,4} = 10;
% runPlot_C1s_msur(model,data,args)

% save('Calibration_results_RAY121_4.mat','-v7.3');

runPlot_NOP(model,data,args) % For group fit plotting





