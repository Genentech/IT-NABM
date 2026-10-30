%% Simulate the NOP model and run the selected scenario(s)
% Top-level driver for the Neuro Oligo Platform (NOP): builds the model, loads
% the cleaned data and parameters, adds variants and doses, loads the task
% table, and runs the selected runModel_NOP_* scenario(s). The large commented
% lists below are menus of selectable scenarios/calibrations - uncomment the
% one(s) to run.
% Models run on the default hour time unit and the data is all transformed to
% have hour time units.

clc
close all
close(findall(0, 'type', 'figure')); % closes fit figures
clear all
warning('off','all')
% restoredefaultpath;

file_path = pwd;

% repo location
addpath(genpath(file_path)); 

% Read a task log for scenarios to run and simulate
args = loadArgs([file_path,'/args_NOP.xlsx']);
args.projLoc = file_path;


% switches for specific simulations
args.switches.test = 1;

%% Build or load the model
% model = buildModel_avidity_v1(locations,switches);
% modelSbproj = sbioloadproject(args.modelLoc); % load the sbioproj file
% model = modelSbproj.m1; % Main model for the effort
model = buildModel_NOP(args); % Build model
% model = buildModel_NOP_delay(args);


%% Load the data
% tall data for simulation/data overlay
data = loadClean_dataset_NOP(args);

% read the wide dataset that is similar to the grouped data used for model calibration
dataWide = readtable(args.dataWideLoc);

%% Load parameters
args.params = readtable(args.paramLoc);

%% Add Variant
model = addVariant_NOP(model,data,args);

% takes name of the "Name" and "Value" columns to generate the variant
model = add_variantParamTable(model, args,'Name','Value'           ,Name='varFit');

%% Add doses
model = addDose_NOP(model,data,args);

%% load the tasks
args = loadTasks(model,data,args);
args = addDose2taskTable(model,data,args);

args.UseParallel = 1;

%% Run the model for the tasks in log and plot the data if available


%% Calibrate the model
% Groups to be used for calibrated in args.groupFit
% MAPT ASO Cyno EastonA - [250805001]
% MAPT ASO Mouse EastonA - [250801001:250801004];
 % Biogen - [250808001:250808004];
 % Monine - [250812001:250812003]


% args.groupFit = [250811001:250811004];

 % runModel_NOP_v1(model,data,args);

% runModel_NOP_intMouse(model,data,args);

% runModel_NOP_Monine_cyno(model,data,args);


% runModel_NOP_Mummery_human(model,data,args);

%% HTT
% runModel_NOP_HTT_cyno(model,data,args)

% runModel_NOP_HTT_human(model,data,args)


%% Jacifusen
% runModel_NOP_Jacifusen_mouse(model,data,args)

% runModel_NOP_Jacifusen_cyno(model,data,args)

% runModel_NOP_Jacifusen_human(model,data,args)

%% Caalibration plot
runModel_NOP_cynoPK(model,data,args)

%% Final - IC50 alone from species, rest from human (translation)
%% BIIB080

% runModel_NOP_BIIB080(model,data,args)

%% APP siRNA
% runModel_NOP_APP_siRNA(model,data,args)

%% Tominersen
% runModel_NOP_Tominersen(model,data,args)

%% Jacifusen
% runModel_NOP_Jacifusen(model,data,args)


%% Final - All PD from species, rest from human (translation)
%% BIIB080

% runModel_NOP_BIIB080_2(model,data,args)

%% APP siRNA
% runModel_NOP_APP_siRNA_2(model,data,args)

%% Tominersen
% runModel_NOP_Tominersen_2(model,data,args)

%% Jacifusen
% runModel_NOP_Jacifusen_2(model,data,args)

%% PK calibrations
% runModel_NOP_BIIB080_Monine(model,data,args)

% runModel_NOP_Jacifusen_HumanPK(model,data,args)
% runModel_NOP_HTT_HumanPK(model,data,args)


%% PD calibrations
% runModel_NOP_PD(model,data,args)