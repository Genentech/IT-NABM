% Creates a percentile plot of scan runs
function dataOut = getData(data, options)
% getData  Filter a tall dataset by Group/ID/Analyte/Matrix/Time.
%   dataOut = getData(data, Group=..., ID=..., Analyte=..., Matrix=..., Time=...)
%   returns the rows of DATA matching every supplied (non-empty) filter; analyte
%   and matrix matching are case-insensitive.
%
%   Inputs:
%     data    - tall data table.
%     options - name-value filters (Group, ID, Analyte, Matrix, Time).
%   Outputs:
%     dataOut - filtered table.

% Sample call: dataTemp = getData(data, Group=[1001,1002],ID=[1,2], Analyte='TDB_Cen_ugmL',Time=[2,96])

% TASK:     
    % Update getData to include flitering based on summary first finding
    % groups that tick boxes on study, molecule, dose, species ,etc then
    % run through the tall data tables filters on Group, ID, Analyte,
    % Assay, Time
% TASK: 
    % include assay similar to analyte
   

arguments
    data table  
    options.Group double = []
    options.ID double = []
    options.Analyte string = []
    options.Matrix string = []
    options.Time double = []
end

% data.Analyte = categorical(data.Analyte);

dataOut = data;

if ~isempty(options.Group)
    dataOut = dataOut(find(ismember(dataOut.Group,options.Group)),:);
end
if ~isempty(options.ID)
    dataOut = dataOut(find(ismember(dataOut.ID,options.ID)),:);
end
if ~isempty(options.Analyte)
    dataOut = dataOut(find(ismember(lower(dataOut.Analyte),lower(options.Analyte))),:);
end
if ~isempty(options.Matrix)
    dataOut = dataOut(find(ismember(lower(dataOut.Matrix),lower(options.Matrix))),:);
end
if ~isempty(options.Time)
    dataOut = dataOut(find(ismember(dataOut.Time,options.Time)),:);
end
