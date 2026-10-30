function args = runCalibration(model, args, options)
% runCalibration  Set up and run a SimBiology parameter-estimation fit.
%   args = runCalibration(model,args) assembles a fitproblem from the grouped
%   data, the estimated-parameter table (args.params), the per-group variants and
%   doses (args.fitTable), the response map and the calibration settings, runs
%   the fit, and stores the results (and %CV) back into args.calibration.
%
%   Inputs:
%     model   - SimBiology model to calibrate.
%     args    - struct with data/task/param locations and calibrationSettings.
%     options - name-value options (Name).
%   Outputs:
%     args    - settings struct with the calibration results added.

arguments
    model SimBiology.Model
    args struct
    options.Name string = ''
end

% ------Configset-------------------------------------------------------------------
% Configset Adjustment?

% ------Data-------------------------------------------------------------------
% Define a description of the data.
gDObj = load(args.dataGrpLoc);
groupedDataObj = gDObj.groupedDataObj;
group = args.groupFit;
% % ids = args.groupids;

% Set up category variables
opts = detectImportOptions(args.taskLoc);
opts.VariableTypes = ["double",repmat("char",1,numel(opts.VariableTypes)-1)];
T = readtable(args.taskLoc,opts);

for ii=1:length(group)
    groupedDataObj.Cat1(ismember(groupedDataObj.Group,group(ii))) = str2double(T.Cat1{find(T.Group==group(ii))});
    groupedDataObj.Cat2(ismember(groupedDataObj.Group,group(ii))) = str2double(T.Cat2{find(T.Group==group(ii))});

end

% % groupedDataObj.NewCat = zeros(3188,1);
% % groupedDataObj.NewCat2 = zeros(3188,1);
% % groupedDataObj.NewCat(2536:2582,:) = 1;
% % groupedDataObj.NewCat(2880:3002,:) = 2;
% % groupedDataObj.NewCat2(2536:2582,:) = 1;
% % groupedDataObj.NewCat2(2880:3002,:) = 1;

groupedDataObj=groupedDataObj(ismember(groupedDataObj.Group,group)&groupedDataObj.Time>0,:);
% % groupedDataObj=groupedDataObj(ismember(groupedDataObj.Group,group)&ismember(groupedDataObj.ID,ids)&groupedDataObj.Time>0,:);

% ------Parameters-------------------------------------------------------------------
% Define objects being estimated and their initial estimates.
params = args.params;
% replace 'linear' with '' to satisfy the fit function
params.Transformation(find(contains(lower(params.Transformation),'linear')))={''};
include = find(params.Include); 
for i=1:numel(include)
    estimatedInfoObj(i) = estimatedInfo(params.Name(include(i)));
    estimatedInfoObj(i).InitialValue = params.Value(include(i)) ;
    estimatedInfoObj(i).Bounds = [params.LowerBound(include(i)), params.UpperBound(include(i))];
    estimatedInfoObj(i).Transform = char(params.Transformation(include(i)));
    estimatedInfoObj(i).CategoryVariableName = char(params.CategoryVar(include(i)));
    % % estimatedInfoObj(i).CategoryVariableName = 'NewCat';
    % CategoryVariableName can be used to make the parameter category-specific. It can be used for group spicific estimation, male/female, different molecules for pooled data across multiple studies, etc
end

% % estimatedInfoObj(numel(include)).CategoryVariableName = 'NewCat2';

% ------Tasks-------------------------------------------------------------------
group = args.groupFit;
% get tasks with the same order as group
[Lia,Locb] = ismember(args.fitTable.Group,group);
tasks = args.fitTable(Lia,:);
[~,p] = sort(Locb(Lia));
tasks = tasks(p,:);

% tasks = args.fitTable(find(ismember(args.fitTable.Group,group)),:);

% Table of variants for fit.
variantsForFit = tasks.VariantArrayCell; % or VariantArray?

% Build table of doses.
dosesForFit = tasks.DoseArray;

% ------ResponseMap-------------------------------------------------------------------
% Define response information.
responseMapTable = readtable(args.responseMapLoc);
include = find(responseMapTable.Include); 
for i=1:numel(include)
    responseMap{i} = [char(responseMapTable.Model(include(i))),' = ', char(responseMapTable.Data(include(i)))];
end

% ------Algorithm Setting-------------------------------------------------------------------
% TASK: make this setting a function where you can choose the algo and the
% settings chosen are include in that function or table

% Define Algorithm options.
options                    = struct;
options.MaxIterations      = args.calibrationSettings.options.MaxIterations; % 'auto', 500, 1000
options.FunctionTolerance  = 1e-08;
options.MaxStallIterations = args.calibrationSettings.options.MaxStallIterations; % 50
options.MaxTime            = Inf;
options.NumInitialPoints   = 'auto';
options.NumTrialPoints     = 'auto';
options.XTolerance         = 1e-06;
options.LocalSolver        = 'fminsearch';%'trust-region-reflective';

% Build local solver options
localSolverOptions         = optimset;
localSolverOptions.TolX    = 1e-08;
localSolverOptions.TolFun  = 1e-08;
localSolverOptions.MaxIter = 400;

% Set local solver options to scattersearch
options.LocalOptions = localSolverOptions;

% Define fit problem.
f              = fitproblem('FitFunction', 'sbiofit');
f.Data         = groupedDataObj;
f.Model        = model;
f.Estimated    = estimatedInfoObj;
f.ResponseMap  = responseMap;
f.ErrorModel   = args.calibrationSettings.ErrorModel; % 'constant','proportional','combined','exponential'
f.Variants     = variantsForFit;
f.Doses        = dosesForFit;
f.FunctionName = args.calibrationSettings.Algorithm; % 'particleswarm','patternsearch','ga','lsqnonlin','scattersearch'
f.Options      = options;
f.Weights      = args.calibrationSettings.Weights; % Can be function handle weightsFunction = @(y) 1./y.^2;
f.ProgressPlot = true;
f.UseParallel  = args.calibrationSettings.UseParallel;
f.Pooled       = args.calibrationSettings.Pooled; % false, true, This overrides the category variable options in estimatedInfoObj

% ------Run Calibration-------------------------------------------------------------------
% Estimate parameter values.
tic % start the timer

[results, simdataI] = f.fit;

elapsedTime = toc;  % stop the timer and calculate elapsed time
fprintf('The calibration time is %.2f seconds.\n', elapsedTime);

% ------Results-------------------------------------------------------------------
% Configure the names of the resulting simulation data.
for i = 1:length(simdataI)
    simdataI(i).Name = 'SimData Individual';
end

% add %CV to estimates
if ~isempty(results.ParameterEstimates)
    args.calibration.ParamEstimates = results.ParameterEstimates;
    args.calibration.ParamEstimates.CV_percent = results.ParameterEstimates.StandardError./results.ParameterEstimates.Estimate *100;
end

% Assign output arguments.
args.calibration.results  = results;
args.calibration.simdataI = simdataI;

end


