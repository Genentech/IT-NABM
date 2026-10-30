%% Add a variant based on the model fits

function model = add_variantParamTable(model,args,nameColLabel,valueColLabel,options)
% add_variantParamTable  Build a SimBiology variant from a parameter table.
%   model = add_variantParamTable(model,args,nameColLabel,valueColLabel) creates
%   a variant (named 'varFit' by default, or options.Name) whose content is taken
%   from the args.params table: each non-NaN row sets the parameter named in the
%   nameColLabel column to the value in the valueColLabel column.
%
%   Inputs:
%     model         - SimBiology model to add the variant to.
%     args          - struct with the params table (args.params).
%     nameColLabel  - name of the column holding parameter names.
%     valueColLabel - name of the column holding parameter values.
%     options       - name-value options (Name).
%   Outputs:
%     model         - model with the variant added.

arguments
    model SimBiology.Model
    args struct
    nameColLabel string
    valueColLabel string
    options.Name string = ''
end

% Select the variant name
if (options.Name == '')
    variantName = 'varFit';
else
    variantName = options.Name;
end

% create the variant
variantObj = addvariant(model,variantName);
for i=1:length(args.params.(nameColLabel))
    if ~isnan(args.params.(valueColLabel)(i))
        addcontent(variantObj, {'parameter', args.params.(nameColLabel){i}, 'Value', args.params.(valueColLabel)(i)});
    end
end