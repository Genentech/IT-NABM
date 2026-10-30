function args = loadTasks(model,data,args,options)
% loadTasks  Resolve the task table's variant names into SimBiology variants.
%   args = loadTasks(model,data,args) reads the task spreadsheet (args.taskLoc),
%   keeps the populated Variant* columns, replaces each variant name with the
%   matching SimBiology variant object (empty cells use varNull), builds a
%   combined variant-array column, and stores the result in args.taskTable and
%   args.fitTable.
%
%   Inputs:
%     model   - SimBiology model holding the named variants.
%     data    - dataset table.
%     args    - struct of project settings (args.taskLoc).
%     options - name-value options (Name).
%   Outputs:
%     args    - settings struct with taskTable and fitTable added.

arguments
    model SimBiology.Model
    data table
    args struct
    options.Name string = ''
end

%% Improt all the variant names as string 
opts = detectImportOptions(args.taskLoc);
opts.VariableTypes = ["double",repmat("char",1,numel(opts.VariableTypes)-1)];
T = readtable(args.taskLoc,opts);

% Only keep non-empty columns
varCol = find(contains(lower(T.Properties.VariableNames), lower({'Variant','DisplayName','Color'})));

% find the non-empty columns
for j=1:numel(varCol)
    if all(cellfun('isempty', T{:,varCol(j)} ))
        varCol(j)=nan;
    end
end
% take the non NaN
varCol = varCol(~isnan(varCol));

% reduce the taskTable to only nonempty columns
T = T(:,[1,varCol]);

% What to put in the empty cells
varNull = getvariant(model,'varNull'); % varNull = '';

%% Variants for the fit
for i=1:numel(T.Group)
    if ~isempty(T.VariantFit{i})
        VariantFit(i,1) = getvariant(model,T.VariantFit{i});
    else
        VariantFit(i,1) = varNull;
    end
end
T.VariantFit = VariantFit;

%% Variants for different molecules
for i=1:numel(T.Group)
    if ~isempty(T.VariantMolecule{i})
        VariantMolecule(i,1) = getvariant(model,T.VariantMolecule{i});
    else
        VariantMolecule(i,1) = varNull;
    end 
end
T.VariantMolecule = VariantMolecule;

%% Variants for different species
for i=1:numel(T.Group)
    if ~isempty(T.VariantSpecies{i})
        VariantSpecies(i,1) = getvariant(model,T.VariantSpecies{i});
    else
        VariantSpecies(i,1) = varNull;
    end    
end
T.VariantSpecies = VariantSpecies;

%% Other variants that are comma separated
varCol = find(contains(lower(T.Properties.VariableNames), lower('VariantOthers')));

% loop over VariantOthers
for j=varCol
        for i=1:numel(T.Group)
            if ~isempty(T.(T.Properties.VariableNames{j}){i})
                VariantTemp(i,1) = getvariant(model,char(T{i,j}));
            else
                VariantTemp(i,1) = varNull;
            end    
        end
    
    T.(T.Properties.VariableNames{j}) = VariantTemp;
    clear VariantTemp
end

%% Add a column that contains all the varaints combined - not useful with scenarios as it cannot take cell array as input
varCol = find(contains(T.Properties.VariableNames, 'Variant'));
T2 = mergevars(T,varCol,NewVariableName='VariantArray');
for i=1:numel(T.Group)
    VariantArray{i,1} = [T2.VariantArray(i,:)];
end
T2.Group = T.Group;
T2.VariantArrayCell = VariantArray;

%% Pass the table back to args
args.taskTable = T;
args.fitTable = T2;
