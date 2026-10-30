%% load taskLog variable with locations and settings

function taskLog = loadTaskLog(location)
% loadTaskLog  Read a task-log spreadsheet into a table.
%   taskLog = loadTaskLog(location) reads the task-log file at LOCATION using
%   auto-detected import options and returns it as a table.
%
%   Inputs:
%     location - path to the task-log file.
%   Outputs:
%     taskLog  - table of logged tasks.

% Get the options in case any changes are needed
opts = detectImportOptions(location);
% Import the data
taskLog = readtable(location,opts);

