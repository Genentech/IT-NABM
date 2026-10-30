%% Run the Neuro Oligo Platform (NOP) model and generate plots
% Alternative top-level driver: builds the NOP model, loads the (pre-built)
% data and parameters, adds variants and doses, loads the task table, and runs
% the diagnostic/scenario plots for selected groups.
% Models run on the default hour time unit and the data is all transformed to
% have hour time units.

%% NOTICE: Go to the parent folder where this repo and tools repos is saved

clc
close all
close(findall(0, 'type', 'figure')); % closes fit figures too
clear all
warning('on','all'); % for debugging
restoredefaultpath;

% Location of the repo
addpath(genpath('./tools')); 

% Read a task log for scenarios to run and simulate
% a sheet with Locations 
% a sheet with Tasks

args = loadArgs('./NeuroOligoPlatform/args_NOP.xlsx');

% switches for specific simulations
% args.switches.test = 1;

%% Build or load the model
model = buildModel_NOP(args);

% TASK: Update the model with CNS PD
% speices specific parameters and flagged translation

%% Load the data
% % building the databse takes a few more seconds but updates from the spider
% % and source google sheets
% data = loadClean_dataset_NOP(args);

% read the pre-built data table
data = readtable(args.dataTallLoc);

% read the wide dataset that is similar to the grouped data used for model calibration
% dataWide = readtable(args.dataWideLoc);

% Data QC example
% dataTemp = getData(data, Group=210015103, Analyte='PK', Matrix={'CSF'}, Time==0)

%% Load parameters
args.params = readtable(args.paramLoc);

%% Add Variants
model = addVariant_NOP(model,data,args);
% model = addVariant_2TC_bl(model,data,args);

% takes name of the "Name" and "Value" columns to generate the variant
model = add_variantParamTable(model, args,'Name','Value'           ,Name='varFit'); % adds values column

% TASK: Make it so that all the param columns are added as variants to the model
model = add_variantParamTable(model, args,'Name','VarBaseMouse'    ,Name='VarBaseMouse'); % 
model = add_variantParamTable(model, args,'Name','VarBaseMonkey'   ,Name='VarBaseMonkey'); % 
model = add_variantParamTable(model, args,'Name','VarBaseHuman'    ,Name='VarBaseHuman'); % 

% model = commit_fit2varFit(model, args,varOut='varFit',varIn='Fit210151'); % adds/overwrites VarIn on varOut

%% Add doses
model = addDose_NOP(model,data,args);

%% load the tasks
args = loadTasks(model,data,args); % the varaints are layered on from left to right (right overwrites)
args = addDose2taskTable(model,data,args);

%% Run and plot the model diagnostics with a specified variant & dose
% It can be used in different ways
% 1. a group (prebuilt in taskTable) and run the model and plot 
% 2. a doseCellArray and a varaintCellArray and run the model and plot

% Diagnosis
% runPlot_Diag_NOP(model,data,args,Type='doseParamGroup',group=101,Format='Oligo',Route='ICV',stopTimeDay=30*6,DoseAmountMpk=10/3.8,DoseIntervalDay=100,DoseInfusionDurMin=1,DoseRepeatCount=0) % 
runPlot_Diag_NOP(model,data,args,Type='doseParamGroup',group=102,Format='Oligo',Route='IT' ,stopTimeDay=30*6,DoseAmountMpk=10/3.8,DoseIntervalDay=100,DoseInfusionDurMin=1,DoseRepeatCount=0) % 
% runPlot_Diag_NOP(model,data,args,Type='doseParamGroup',group=103,Format='OligoBrS',Route='IV',stopTimeDay=30*6,DoseAmountMpk=20,DoseIntervalDay=100,DoseInfusionDurMin=1,DoseRepeatCount=0) % gD
% runPlot_Diag_NOP(model,data,args,Type='doseParamGroup',group=104,Format='OligoBrS',Route='IV',stopTimeDay=30*6,DoseAmountMpk=20,DoseIntervalDay=100,DoseInfusionDurMin=1,DoseRepeatCount=0) % 

%% Run the model and plot the results per study
% for each new simulation the data needs to be added, dose and variants
% need to be added to addDose and addVariant and task table and a runplot
% needs to be created

% Cyno scenarios
% runPlot_C_scenarios(model,data,args)

% Run a predefined scenario simulation from tasks
% runPlot_Scenario_NOP(model,data,args,Type='doseParamGroup',group=101,Format='Oligo',Route='ICV',stopTimeDay=30*6,DoseAmountMpk=100/3.8,DoseIntervalDay=100,DoseInfusionDurMin=1,DoseRepeatCount=0) % gD
