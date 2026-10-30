function args = runPRCC(model, args, options)
% runPRCC  Run a partial rank correlation (PRCC) sensitivity analysis.
%   args = runPRCC(model,args) builds an LHS sample of the included PSA
%   parameters (args.paramsPSA), simulates the model across the sampled scenarios
%   for the selected group, and computes Spearman partial rank correlation
%   coefficients between each input parameter and the end-time outputs. Results
%   are stored in args.PRCC.
%
%   Inputs:
%     model   - SimBiology model.
%     args    - struct with the PSA parameter table and PRCC settings.
%     options - name-value options (Name, outputTime).
%   Outputs:
%     args    - settings struct with the PRCC results added.

arguments
    model SimBiology.Model
    args struct
    options.Name string = ''
    options.outputTime double = [0:1:24*7*4]
end

tic % start the timer

% ------Configset-------------------------------------------------------------------

% Configset Adjustment?


% ------Set up the task-------------------------------------------------------------------

% find the task or group that is run for this purpose
group = args.PRCC.Group;
tasks = args.taskTable(find(ismember(args.taskTable.Group,group)),:);

% Add group similar to gQSPsim task setups
% sObj = SimBiology.Scenarios;
% add(sObj,'elementwise','group',tasks.Group);

sObj = SimBiology.Scenarios();
% Configure RandomSeed to a unique value.
sObj.RandomSeed = sum(double('jafarnem'));

% Add Variants as elementwise
varCol = find(contains(lower(tasks.Properties.VariableNames), lower('Variant')));
for i=varCol
    add(sObj,'elementwise','variants',repmat(tasks.(tasks.Properties.VariableNames{i})(:), args.PRCCsampleSize,1));
end

% Add doses as elementwise
doseCol = find(contains(lower(tasks.Properties.VariableNames), lower('Dose')));
for i=doseCol
    add(sObj,'elementwise','dose',repmat(tasks.(tasks.Properties.VariableNames{i})(:), args.PRCCsampleSize,1));
end

% ------Input Parameters-------------------------------------------------------------------

% Define objects being estimated and their initial estimates.
paramPRCC = args.paramsPSA(args.paramsPSA.Include==1,:);

% TASK: add other distributions, add a column on, distName, mu and sigma to parameters and truncate the distribution with LB and UB 
% Generate scenarios.
for i=1:numel(paramPRCC.Include)
    if contains(paramPRCC.Transformation(i),'log')
        tempDistribution = 'Uniform';
    else
        tempDistribution = 'loguniform';
    end
    paramPRCC.pd(i) = makedist(tempDistribution, 'Lower', paramPRCC.LowerBound(i), 'Upper', paramPRCC.UpperBound(i));
end
add(sObj, 'elementwise', paramPRCC.Name', paramPRCC.pd', 'Number', args.PRCCsampleSize, 'SamplingMethod', 'lhs');


% Define sensitivity inputs.
scenariosTbl = generate(sObj);
inputTable = scenariosTbl(:,paramPRCC.Name');
input = table2array(scenariosTbl(:,paramPRCC.Name'));

% ------run the simulations-------------------------------------------------------------------

% Run the model using the simFunc
simFunc = createSimFunction(model,sObj,args.PRCC.outputNames,[],[],'UseParallel',args.UseParallel, 'AutoAccelerate',1);
simDataArray = simFunc(sObj, [], [], [0:1:options.outputTime]);

% ------Output Parameters-------------------------------------------------------------------

% Define sensitivity outputs.
for i=1:numel(args.PRCC.outputNames)
    [~,x] = selectbyname(simDataArray,args.PRCC.outputNames{i});
    x = cell2mat(x');
    output(:,i) = x(end,:)';
end

% ------Calc PRCC-------------------------------------------------------------------

% Calculate Partial Rank Correlation Coefficients
[rho, pval] = partialcorr(input,output,ones(size(output)),type='Spearman',rows='complete');


% ------Results-------------------------------------------------------------------

% Populate the output structure.
args.PRCC.input   = input;
args.PRCC.output  = output;
args.PRCC.rho     = rho;
args.PRCC.pval    = pval;
args.PRCC.inputNames = inputTable.Properties.VariableNames;
args.PRCC.inputNamesDisplay = paramPRCC.DisplayName';
args.PRCC.scenarios  = scenariosTbl;
args.PRCC.simresults = simDataArray;


elapsedTime = toc;  % stop the timer and calculate elapsed time
fprintf('PRCC analysis time is %.2f seconds.\n', elapsedTime);

end

