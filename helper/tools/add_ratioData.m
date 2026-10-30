function dataOut = add_ratioData(data,options)
% add_ratioData  Append ratio (numerator/denominator analyte) rows to a dataset.
%   dataOut = add_ratioData(data, Group=..., Matrix=..., AnalyteNom=...,
%   AnalyteDenom=...) divides, per Group/Matrix/subject ID and matched time, the
%   numerator analyte concentration by the denominator analyte concentration and
%   appends the result with a combined analyte name and the unit set to 'ratio'.

% Sample call: data = add_ratioData(data,Group=[210015101,210015102,210015103], Matrix={'CSF','Plasma'}, AnalyteNom='pC4a', AnalyteDenom='C4');
   

arguments
    data table  
    options.Group double = []
    options.ID double = []
    options.AnalyteNom string = []
    options.AnalyteDenom string = []
    options.Matrix string = []
end

dataOut = data;

for i=1:numel(options.Group)
    for j=1:numel(options.Matrix)
            dataNom   = getData(data, Group=options.Group(i), Matrix=options.Matrix{j}, Analyte=options.AnalyteNom);
            dataDenom = getData(data, Group=options.Group(i), Matrix=options.Matrix{j}, Analyte=options.AnalyteDenom);
            uniqueID = unique(dataNom.ID);
            for p=1:numel(uniqueID)
                dataNomTemp   = dataNom(dataNom.ID==uniqueID(p),:);
                dataDenomTemp = dataDenom(dataDenom.ID==uniqueID(p),:);
                [~, idIndex] = ismember(dataDenomTemp.Time, dataNomTemp.Time);
                dataNomTemp.Conc = dataNomTemp.Conc./dataDenomTemp.Conc(idIndex);
                dataNomTemp.Analyte = strcat(dataNomTemp.Analyte,repmat({'_'},numel(dataNomTemp.Analyte),1),dataDenomTemp.Analyte);
                dataNomTemp.Conc_Unit = repmat({'ratio'},numel(dataNomTemp.Analyte),1);
                dataOut = [dataOut; dataNomTemp];
            end
    end
end

