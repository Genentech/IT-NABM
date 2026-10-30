function model = updateVarFit(model, args, options)
% updateVarFit  Apply the latest calibration estimates to the 'varFit' variant.
%   model = updateVarFit(model,args) snapshots the current varFit as
%   'varFitPreCalib' (once), then updates varFit with the calibrated parameter
%   estimates (args.calibration.results), adding any parameters not already
%   present.
%
%   Inputs:
%     model   - SimBiology model holding the varFit variant.
%     args    - struct with the calibration results.
%     options - name-value options (Name).
%   Outputs:
%     model   - model with varFit updated.

arguments
    model SimBiology.Model
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
