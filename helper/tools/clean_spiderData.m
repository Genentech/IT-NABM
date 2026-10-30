function t1 = clean_spiderData(t1)
% clean_spiderData  Rename/standardise Spider-export columns to the model format.
%   t1 = clean_spiderData(t1) renames and reformats the columns of a Spider query
%   export (Group, ID, Time, Analyte, Matrix, Conc, units, Assay, etc.) into the
%   modelling-friendly naming and value conventions used downstream.
%
%   Inputs/Outputs:
%     t1 - table of Spider data, returned standardised.

% rename spider variables so that they match the modeling friendly format and names
% Group
t1 = renamevars(t1,{'GROUP'},{'Group'});
% ID
for i=1:numel(t1.USUBJID)
    numbers = str2double(extract(t1.USUBJID(i), digitsPattern));
    t1.ID(i) = numbers(end); 
end
% Time
% t1 = renamevars(t1,{'PCTPTNUM'},{'Time'});
t1 = renamevars(t1,{'PCTPT'},{'Time'});
t1.Time = regexprep(t1.Time, ' day', '');
t1.Time = cellfun(@str2double, t1.Time);

% Time_Unit
t1.Time_Unit(:) = {'Day'};
% Analyte
t1.Analyte = t1.PCTEST;
t1.Analyte(find(contains(t1.Analyte,'IGG'))) = {'PK'};
% Matrix
t1 = renamevars(t1,{'PCSPEC'},{'Matrix'}); 
t1.Matrix = regexprep(t1.Matrix, 'FLUID, CEREBROSPINAL', 'CSF');
t1.Matrix = regexprep(t1.Matrix, 'PLASMA', 'Plasma');
t1.Matrix = regexprep(t1.Matrix, 'SERUM', 'Serum');
t1.Matrix = regexprep(t1.Matrix, ' D29', '');
t1.Matrix = regexprep(t1.Matrix, ' D3', '');
t1.Matrix = regexprep(t1.Matrix, ' ', '_');
% Conc
t1 = renamevars(t1,{'PCORRES'},{'Conc'});
% Conc_SD which is NaN for measurements (used for digitized values)
t1.Conc_SD(:) = NaN;
% Conc_Unit
t1 = renamevars(t1,{'PCORRESU'},{'Conc_Unit'}); %t1.Conc_SD = str2double(t1.Conc_SD);
% Assay or method of measurement
t1 = renamevars(t1,{'PCMETHOD'},{'Assay'});
% make the study_id column categorical for easier filtering
t1.study_id = categorical(t1.study_id);