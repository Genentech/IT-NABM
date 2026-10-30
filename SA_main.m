%% Sensitivity analysis driver (one parameter at a time)
% Builds the NOP model, loads data/parameters, adds variants and doses, then
% runs a series of one-at-a-time sensitivity sweeps (brain oligo half-life,
% brain protein half-life, IC50, cell uptake, mRNA half-life) by setting
% args.ParaVar / args.ParaVarValues and calling runPlot_SA_scenarios for each.
% Models run on the default hour time unit and the data is all transformed to
% have hour time units.

clc
close all
close(findall(0, 'type', 'figure')); % closes fit figures
clear all
warning('off','all')
restoredefaultpath;

file_path = pwd;

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

%% Build or load the model
model = buildModel_NOP(args); % Build model

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


%% Add doses
model = addDose_NOP(model,data,args);

%% load the tasks
args = loadTasks(model,data,args);
args = addDose2taskTable(model,data,args);


% Define the parameter for analysis and the range (restricted to 10 values)
args.ParaVar = "t_hl_Oligo_b"; 
args.ParaVarValues = [16,26,36,46,56]; 



% Brain oligo half life
args.ParaVar = "t_hl_Oligo_b";
args.ParaVarValues = [16,26,36,46,56]; 

runPlot_SA_scenarios(model,data,args)
set(gcf,"Position",[100,100,400,400])

% Brain protein half life
args.ParaVar = "t_hl_Prot_b"; 
args.ParaVarValues = [57, 86, 114, 143, 171]; 

runPlot_SA_scenarios(model,data,args)
set(gcf,"Position",[100,100,400,400])


% IC50
args.ParaVar = "IC50_mRNA_Oligo_b"; 
args.ParaVarValues = [1,10,162,500,1000]; 

runPlot_SA_scenarios(model,data,args)
set(gcf,"Position",[100,100,400,400])

% Cell uptake
args.ParaVar = "t_hl_cellUptake"; 
args.ParaVarValues = [0.04,0.2,1,5,25]; 

runPlot_SA_scenarios(model,data,args)
set(gcf,"Position",[100,100,400,400])

% mRNA half life
args.ParaVar = "t_hl_mRNA"; 
args.ParaVarValues = [1,12,24,48,120]; 

runPlot_SA_scenarios(model,data,args)
set(gcf,"Position",[100,100,400,400])