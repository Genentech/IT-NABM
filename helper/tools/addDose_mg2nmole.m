function model = addDose_mg2nmole(model,options)
% addDose_mg2nmole  Add a scheduled dose, converting mg/kg amounts to nanomoles.
%   model = addDose_mg2nmole(model, Name=..., Time_hr=..., Dose_mgkg=..., MW=...,
%   Infusion_min=..., SpeciesName=..., Route=...) builds a schedule dose table by
%   converting each mg dose amount to nanomoles using the molecular weight MW and
%   the infusion rate from Infusion_min, then adds it to the model targeting the
%   given species (or a route-default target).
%
%   Inputs:
%     model   - SimBiology model to add the dose to.
%     options - name-value options (Name, Time_hr, Dose_mgkg, Infusion_min, MW,
%               SpeciesName, Route).
%   Outputs:
%     model   - model with the scheduled dose added.

arguments
    model SimBiology.Model
    options.Name string = 'dose'
    options.Time_hr double = 0
    options.Dose_mgkg double = 0.0
    options.Infusion_min double = 0.01
    options.MW double = 150
    options.SpeciesName string = ''
    options.Route string = 'IV'
end

%------ Set up the dose table with correct units --------------
% set up the time and dose amount values correctly
% Dose Time in hour
Time = options.Time_hr(:);

% Dose Amount
if numel(Time)==numel(options.Dose_mgkg)
    Amount = options.Dose_mgkg(:);
elseif numel(options.Dose_mgkg)==1 
    Amount = repmat(options.Dose_mgkg,numel(Time),1);
else
    display('Dose_mgkg must have one or same number of elements as Time_hr');
    return
end
% Amount in mg/(MW in kg/mole) = umole*1000nano/u = nanomle 
Amount = Amount/options.MW*1000;

% Rate from infusion duration
if numel(Time)==numel(options.Infusion_min)
    Rate = Amount./options.Infusion_min(:);
elseif numel(options.Infusion_min)==1 
    Rate = Amount./repmat(options.Infusion_min,numel(Time),1);
else
    display('Infusion_min must have one or same number of elements as Time_hr');
    return
end
% Rate in nanomole/min *60min/hour = nanomole/hour
Rate = Rate*60;

% create the dose table and fix the units
tbl = table(Time,Amount,Rate);

%----------- Select What species is dosed -------------- 
if ~(options.SpeciesName == '')
    TargetName = options.SpeciesName;
elseif strcmpi(options.Route,'IV')
    TargetName = 'body.AOC_plasma_pmole';
elseif strcmpi(options.Route,'SC')
    TargetName = 'body.AOC_sc_pmole';
end

%------------- Create the dose object-------------------
dose  = adddose(model,options.Name,"schedule");
setTable(dose,tbl);
dose.TargetName      = TargetName;
dose.AmountUnits     = 'nanomole';
dose.TimeUnits       = 'hour';
dose.RateUnits       = 'nanomole/hour';




