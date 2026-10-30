%% load args variable with locations and settings

function args = loadArgs(location)
% loadArgs  Read the project settings/locations spreadsheet into a struct.
%   args = loadArgs(location) reads the args spreadsheet at LOCATION and returns
%   a struct whose fields are the entries of the Name column. Values typed as
%   "double" or "str2num" are converted accordingly, any field whose name
%   contains 'Loc' is resolved relative to args.projLoc, and the figure output
%   folder is created if it does not already exist.
%
%   Inputs:
%     location - path to the args .xlsx file.
%   Outputs:
%     args     - settings struct.

opts = detectImportOptions(location);

% Specify range and varTypes
opts.DataRange = "A2";
opts.VariableTypes = ["string", "string", "string", "string"];
% Specify variable properties
opts = setvaropts(opts, ["Name", "Value", "Type", "Description"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Name", "Value", "Type", "Description"], "EmptyFieldRule", "auto");

%% Import the data
T = readtable(location,opts);

%% add arguments 
for i=1:numel(T.Name)
    args.(T.Name(i)) = T.Value(i);
end

%% set the doubles
for i=1:numel(T.Name)
    if strcmpi(T.Type{i},"double")
        args.(T.Name(i)) = double(args.(T.Name(i)));
    end
end

%% set the str2num
for i=1:numel(T.Name)
    if strcmpi(T.Type{i},"str2num")
        args.(T.Name(i)) = str2num(args.(T.Name(i)));
    end
end

%% Add the project location to the paths
% Any variable name that has Loc will be added to the proj path
for i=1:numel(T.Name)
    if contains(T.Name{i},'Loc')&& ~strcmpi(T.Name{i},"projLoc")
        args.(T.Name(i)) = fullfile(args.projLoc,args.(T.Name(i)));
    end
end

%% Create the folder if doesn't exist already
if not(isfolder(args.figLoc))
    mkdir(args.figLoc)
end
% Add the folders to Matlab paths
% addpath(genpath(args.projLoc)); 

