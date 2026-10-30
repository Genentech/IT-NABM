function saveCalibration(args, options)
% saveCalibration  Write fitted parameter estimates into the parameters table.
%   saveCalibration(args, Name=..., Overwrite=...) takes the latest calibration
%   estimates (args.calibration.results) and writes them as a named column in the
%   parameters table (args.params), then saves the table back to args.paramLoc.
%   An existing column is only replaced when Overwrite is true.
%
%   Inputs:
%     args    - struct with the calibration results and the params table.
%     options - name-value options (Name, Overwrite).

arguments
    args struct
    options.Name string = ''
    options.Overwrite logical = 0
end


%% Update the varFit with the new calibration
% Calibrated parameters
ParameterEstimates = args.calibration.results.ParameterEstimates;

% check if the specified name is unique
args.params.Properties.VariableNames;
if ismember(options.Name,args.params.Properties.VariableNames)&(~options.Overwrite)
    display('The specified name already exists in the parameters table');
    return
end
if ismember(options.Name,args.params.Properties.VariableNames)&(options.Overwrite)
    % remove the old params column
    args.params = args.params(:,find(~ismember(args.params.Properties.VariableNames,options.Name)));
    display(['Overwrote the calibrated parameters in ',char(options.Name),' column']);
end
%% write the new table column with only fitted parameters in the right locations
% go over all the estimated parameters and either update varFit value or
% add the parameter to the varFit
for i=1:numel(ParameterEstimates.Name)
    [logic,index] = ismember(ParameterEstimates.Name{i},args.params.Name);
    if logic
        args.params.(options.Name){index,1} = ParameterEstimates.Estimate(i);
    else
        display(['Fitted parameter ',num2str(i),' does not exist in the paramlist']);
    end
end

writetable(args.params ,args.paramLoc,"WriteMode","overwritesheet");

