function args = runSobol(model, args, options)
% runSobol  Run a Sobol variance-based global sensitivity analysis.
%   args = runSobol(model,args) builds a Sobol sample of the included PSA
%   parameters (args.paramsPSA) and calls sbiosobol to compute first- and
%   total-order Sobol indices for the requested outputs over the selected group's
%   variants and doses. Results are stored in args.Sobol.
%
%   Inputs:
%     model   - SimBiology model.
%     args    - struct with the PSA parameter table and Sobol settings.
%     options - name-value options (Name, sobolOutputs, outputTime).
%   Outputs:
%     args    - settings struct with the Sobol results added.

arguments
    model SimBiology.Model
    args struct
    options.Name string = ''
    options.sobolOutputs cell = {'V_bsf.AUC_D_nM'}
    options.outputTime double = [0:1:24*7*4]
end

% ------Configset-------------------------------------------------------------------
% Configset Adjustment?

% ------Input Parameters-------------------------------------------------------------------
% Define objects being estimated and their initial estimates.
paramSobol = args.paramsPSA(args.paramsPSA.Include==1,:);
% paramSobol.Name
% [paramSobol.LowerBound, paramSobol.UpperBound]

% TASK: add other distributions, add a column on, distName, mu and sigma to parameters and truncate the distribution with LB and UB 
% Generate scenarios.
samples1 = SimBiology.Scenarios();
for i=1:numel(paramSobol.Include)
    if contains(paramSobol.Transformation(i),'log')
        tempDistribution = 'Uniform';
    else
        tempDistribution = 'loguniform';
    end
    paramSobol.pd(i) = makedist(tempDistribution, 'Lower', paramSobol.LowerBound(i), 'Upper', paramSobol.UpperBound(i));
end
add(samples1, 'elementwise', paramSobol.Name', paramSobol.pd', 'Number', args.sobolSampleSize, 'SamplingMethod', 'sobol');
% 
% 
% pd1      = makedist('Uniform', 'Lower', 1, 'Upper', 10);
% pd2      = makedist('Loguniform', 'Lower', 0.1, 'Upper', 3);
% pd3      = makedist('Lognormal', 'mu', 0, 'sigma', 1);
% add(samples1, 'elementwise', {'CL_C1sRec', 'Vm_TfR', 'f_TfR_BBB'}, [pd1 pd2 pd3], 'Number', 1000, 'SamplingMethod', 'sobol');

% Configure RandomSeed to a unique value.
samples1.RandomSeed = sum(double('jafarnem'));

% Define sensitivity inputs.
inputs = samples1;

% ------Output Parameters-------------------------------------------------------------------
% Define sensitivity outputs.
outputs = options.sobolOutputs;

% ------Tasks-------------------------------------------------------------------
group = args.sobolGroup;
tasks = args.fitTable(find(ismember(args.fitTable.Group,group)),:);

% Table of variants for sobol.
variantsForSobol = tasks.VariantArrayCell; % or VariantArray?

% Build table of doses for sobol.
dosesForSobol = tasks.DoseArray;

% ------Algorithm Setting-------------------------------------------------------------------
% Define options.
options = {'Variants', variantsForSobol{:}, ...
           'Doses', dosesForSobol{:},...
           'OutputTime', options.outputTime,...
           'ShowWaitBar',true};

% ,...
%            'Bounds', table2array(paramSobol(:,{'LowerBound','UpperBound'}))


% ------Run Sobol-------------------------------------------------------------------
% Estimate parameter values.
tic % start the timer

% Compute first and total order Sobol indices.
sobolData = sbiosobol(model, inputs, outputs, options{:});

elapsedTime = toc;  % stop the timer and calculate elapsed time
fprintf('Sobol analysis time is %.2f seconds.\n', elapsedTime);

% ------Results-------------------------------------------------------------------

% Populate the output structure.
args.Sobol.results    = sobolData;
args.Sobol.simresults = sobolData.SimulationInfo.SimData;

end
