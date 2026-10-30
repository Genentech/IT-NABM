function model = updateVarFit_plot(model, group, args, options)
% updateVarFit_plot  Apply a single group's calibration estimates to 'varFit'.
%   model = updateVarFit_plot(model,group,args) is like updateVarFit but first
%   restricts the calibrated estimates to the given group's category value before
%   updating the varFit variant (used when plotting per-group fits).
%
%   Inputs:
%     model   - SimBiology model holding the varFit variant.
%     group   - group ID whose category-specific estimates to apply.
%     args    - struct with the calibration results.
%     options - name-value options (Name).
%   Outputs:
%     model   - model with varFit updated.

arguments
    model SimBiology.Model
    group double
    args struct
    options.Name string = ''
end

%% add varFitPreCalib as precalibration parameters
varFit = getvariant(model,'varFit');
varFitPreCalib = copyobj(varFit,model,'force');
% only rename if varFitPreCalib dosen't exist
if isempty(getvariant(model,'varFitPreCalib'))
    varFitPreCalib.Name = 'varFitPreCalib';
end

%% Update the varFit with the new calibration
% Calibrated parameters
ParameterEstimates = args.calibration.results.ParameterEstimates;

ParameterEstimates = ParameterEstimates(ismember(ParameterEstimates.CategoryValue,[categorical(group),categorical(1)]),:);

% convert varFit to a table for ease of search
varFitTable = variant2table(varFit);

% go over all the estimated parameters and either update varFit value or
% add the parameter to the varFit
for i=1:numel(ParameterEstimates.Name)
    [logic,index] = ismember(ParameterEstimates.Name{i},varFitTable.Name);
    if logic
        varFit.content{index}{4} = ParameterEstimates.Estimate(i);
    else
        addcontent(varFit, {'parameter', ParameterEstimates.Name{i} , 'Value', ParameterEstimates.Estimate(i)});
    end
end
