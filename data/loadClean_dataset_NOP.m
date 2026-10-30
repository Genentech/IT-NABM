%% Load the dataset for the C1s program
%% This only needs to be run if the source dataset is changed either on spider with a new query or on GoogleDrive with manual entry

function [data] = loadClean_dataset_NOP(args)
% loadClean_dataset_NOP  Load, clean and reshape the NOP source dataset.
%   data = loadClean_dataset_NOP(args) reads the curated tall dataset and its
%   summary sheet, harmonises units (converts times to hours), removes
%   duplicate rows, builds the Matrix_Analyte_Unit fit key, reshapes the data
%   into a wide grouped-data object for sbiofit, and writes the cleaned
%   tall/wide tables and grouped object back to the data folder. Only needs to
%   be re-run when the source dataset changes.
%
%   Inputs:
%     args - struct of project settings/paths (dataset locations).
%   Outputs:
%     data - cleaned tall data table.

% TASK: Merge summary sheet with the data sheet to include key features in
% one single table for filterand combining the data

%% Read the data
% read the metadata (summary dataset)
tSum = readtable(args.dataCuratedSummaryLoc);

% read the dataset from modeling database on the SharedGoogleDrive
T = readtable(args.dataCuratedLoc);

% read the dataset from Spider format
% t1 = readtable(args.dataSpiderLoc);

%% Clean up the data

data = T(:,1:10);
data = sortrows(data,{'Group','Analyte','Matrix','Assay','ID','Time'});

%% Unit coversion
% convert all the times to hour
iDay = find(ismember(lower(data.Time_Unit),lower("Day")));
data.Time_Unit(iDay) = {'Hour'};
data.Time(iDay) = round(hours(days(data.Time(iDay))),4,"significant");

% convert all ng/mL or ng/g to ug/mL
% iUnit = find(ismember(lower(data.Conc_Unit),lower({'ng_ml','ng/ml','ng_g','ng/g'}))); % assumes 1gr=1mL
% data.Conc_Unit(iUnit) = {'ug/mL'};
% data.Conc(iUnit) = data.Conc(iUnit)/1000;

% iUnit = find(ismember(lower(data.Conc_Unit),lower({'nM'}))&ismember(data.Group,[231000201])); % assumes 1gr=1mL
% data.Conc_Unit(iUnit) = {'ug/mL'};
% data.Conc(iUnit) = data.Conc(iUnit)*150/1000;

%% Remove duplicates, Add normalized values
% remove duplicate rows if it was resulted from joins
[~,index]=unique(data(:,{'Group','ID','Time','Analyte','Matrix','Conc','Assay'}),'rows');
display([num2str(numel(data.Group)-numel(index)),' duplicate rows were found in the data and removed']);
data = data(index,:);

% for ease of derivatives convert all serum to plasma
data.Matrix(find(strcmpi(data.Matrix,'Serum'))) = {'Plasma'}; 


%% Add the fit-ready column
% Clean up the units
unitsClean = regexprep(data.Conc_Unit, '/', '');
unitsClean = regexprep(unitsClean, '_', '');

% Create a column that  combines Matrix and Analyte
data.Matrix_Analyte_Unit = strcat(data.Matrix,'_',data.Analyte,'_',unitsClean);
data.Matrix_Analyte_Unit = regexprep(data.Matrix_Analyte_Unit, 'Serum_PK_ugmL',               'Plasma_PK_ugmL');

%% Censoring and Exclusions only for Fit Data
dataCensor = data;
% dataCensor.Conc(find(dataCensor.Group==210015103 ...
%                 & dataCensor.ID==3003 ...
%                 & ismember(dataCensor.Time,[504,672]) ...
%                 & ismember(dataCensor.Analyte,{'PK'}))) = NaN;
% dataCensor.Conc(find(dataCensor.Group==210015103 ...
%                 & ismember(dataCensor.Analyte,{'pC1s','pC1s_norm'}))) = NaN;

%% Create the Wide data tables
% Create a wide dataset for sbiofit
T = dataCensor(:,{'Group','ID','Time','Conc','Matrix_Analyte_Unit'});
meanNoNan = @(x) mean(x,"omitnan");
[dataWide, ~] = unstack(T,'Conc','Matrix_Analyte_Unit','AggregationFunction',meanNoNan);
% reorder columns for ease of QC
VariableOrdersWanted = {'Group','ID','Time','Plasma_PK_ugmL','Brain_PK_ugmL','CSF_PK_ugmL'};
% dataWide = reorderDataWide(dataWide,VariableOrdersWanted);
% baseline pooled grouping variable
dataWide.Pooled = ones(length(dataWide.Group),1);
dataWide.Study = dataWide.Group;
% add log transformed plasma PK column
% dataWide.Plasma_PK_ugmL_log = log(dataWide.Plasma_IgGPK_nM);

% compress baseline values so that there is no negative time (not accepted by simbiology)
dataWide = removeNegBaseline(dataWide);

% sort for ease of quality check
dataWide = sortrows(dataWide,{'Group','ID','Time'});



%% Create the grouped data object
groupedDataObj                                    = groupedData(dataWide);
groupedDataObj.Properties.GroupVariableName       = 'Group';
groupedDataObj.Properties.IndependentVariableName = 'Time';
groupedDataObj = sortrows(groupedDataObj,{'Group','Time'});

% add the units
VariableUnits = {'dimensionless','dimensionless','hour'};
VariableNames = groupedDataObj.Properties.VariableNames;
for i = 4:length(groupedDataObj.Properties.VariableNames)
    newStr = split(VariableNames{i},'_');
    switch newStr{end}
        case 'ugmL'
            unit = 'microgram/milliliter';
        case 'nM'
            unit = 'nanomolarity';
        case 'percent'
            unit = 'nanomolarity'; % It is actually 'dimensionless' but we use nM as surrogate because species must have a conc or amount unit
        case 'ngml'
            unit = 'nanogram/milliliter';
        otherwise
            unit = 'dimensionless';
    end
    VariableUnits{i} = unit;
end
groupedDataObj.Properties.VariableUnits = VariableUnits;

%% Save the final data tables and objects
% write the table to the clean data folder for faster future queries
% Round only numeric variables for export readability/consistency.
dataTallToWrite = data;
dataWideToWrite = dataWide;

numericVarsTall = dataTallToWrite.Properties.VariableNames( ...
    varfun(@isnumeric, dataTallToWrite, "OutputFormat", "uniform"));
for i = 1:numel(numericVarsTall)
    dataTallToWrite.(numericVarsTall{i}) = round(dataTallToWrite.(numericVarsTall{i}),2);
end

numericVarsWide = dataWideToWrite.Properties.VariableNames( ...
    varfun(@isnumeric, dataWideToWrite, "OutputFormat", "uniform"));
for i = 1:numel(numericVarsWide)
    dataWideToWrite.(numericVarsWide{i}) = round(dataWideToWrite.(numericVarsWide{i}),2);
end

writetable(dataTallToWrite,fullfile(args.dataTallLoc),"WriteMode","overwritesheet");
writetable(dataWideToWrite,fullfile(args.dataWideLoc),"WriteMode","overwritesheet");
save(fullfile(args.dataGrpLoc),"groupedDataObj");

display("database was updated");
