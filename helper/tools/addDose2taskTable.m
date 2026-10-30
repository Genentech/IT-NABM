function args = addDose2taskTable(model,data,args,options)
% addDose2taskTable  Attach each group's dose objects to the task/fit tables.
%   args = addDose2taskTable(model,data,args) finds, for every group in the task
%   table, the model doses named doseGroup<Group>-* and stores them as a dose
%   array on the fit table; it also pads the task table with Dose1..DoseN columns
%   (filled with the null dose) and inserts each group's actual doses.
%
%   Inputs:
%     model   - SimBiology model holding the dose objects.
%     data    - dataset table.
%     args    - struct with taskTable and fitTable.
%     options - name-value options (Name).
%   Outputs:
%     args    - settings struct with the dose columns/arrays added.

arguments
    model SimBiology.Model
    data table
    args struct
    options.Name string = ''
end

%% get the doses related to the groups and add them here
T = args.taskTable;
T2 = args.fitTable;

% Get the doses with the same group code
for i=1:length(T.Group)
    doseTemp = sbioselect(model, 'Where','Name', 'regexp', ['doseGroup',num2str(T.Group(i)),'-*']);
    doseArray{i,1} = [doseTemp(:)];
end
T2.DoseArray = doseArray;

% Find the null dose level
doseNull = getdose(model,'doseNull'); % varNull = '';
% Number of doses specified for each group
doseNum = cellfun(@numel,doseArray);
% Loop over the max number of doses per group and create Null doses for all
for j=1:max(cellfun(@numel,doseArray))    
    colName = ['Dose',num2str(j)];
    T.(colName) = repmat(doseNull,numel(doseNum),1);
end
% Go over doses that actually exists and replace the doseNull with actual
% predefined dose
for i=1:length(T.Group)
    for j=1:doseNum(i)
        colName = ['Dose',num2str(j)];
        T.(colName)(i) = doseArray{i}(j);
    end
end

args.taskTable = T;
args.fitTable = T2;



