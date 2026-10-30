function dataOut = add_ratioDataMatrix(data,options)
% add_ratioDataMatrix  Append cross-matrix ratio rows (e.g. Brain/Plasma) to data.
%   dataOut = add_ratioDataMatrix(data, Group=..., Analyte=..., MatrixNom=...,
%   MatrixDenom=...) divides, per Group/Analyte/subject ID and matched time, the
%   numerator-matrix concentration by the denominator-matrix concentration (as a
%   percentage) and appends the result with a combined matrix name and the unit
%   set to 'percent'.

% Sample call: data = add_ratioDataMatrix(data,Group=unique(data.Group), Analyte={'PK'}, MatrixNom='Brain', MatrixDenom='Plasma');
   
arguments
    data table  
    options.Group double = []
    options.ID double = []
    options.MatrixNom string = []
    options.MatrixDenom string = []
    options.Analyte string = []
end

dataOut = data;

for i=1:numel(options.Group)
    for j=1:numel(options.Analyte)
            dataNom   = getData(data, Group=options.Group(i), Matrix=options.MatrixNom,   Analyte=options.Analyte{j});
            dataDenom = getData(data, Group=options.Group(i), Matrix=options.MatrixDenom, Analyte=options.Analyte{j});
            uniqueID = unique(dataNom.ID);
            for p=1:numel(uniqueID)
                dataNomTemp   = dataNom(dataNom.ID==uniqueID(p),:);
                dataDenomTemp = dataDenom(dataDenom.ID==uniqueID(p),:);
                uniqueTime = unique(dataNomTemp.Time);
                for k=1:length(uniqueTime)
                    dataNomTempTime   = dataNomTemp(dataNomTemp.Time==uniqueTime(k),:);
                    dataDenomTempTime = dataDenomTemp(dataDenomTemp.Time==uniqueTime(k),:);
                    [~, idIndex] = ismember(dataDenomTempTime.Time, dataNomTempTime.Time);
                    if ~isempty(idIndex)
                        dataNomTempTime.Conc = dataNomTempTime.Conc./dataDenomTempTime.Conc(logical(idIndex))*100;
                        dataNomTempTime.Matrix = strcat(dataNomTempTime.Matrix,repmat({'2'},numel(dataNomTempTime.Matrix),1),dataDenomTempTime.Matrix(logical(idIndex)));
                        dataNomTempTime.Conc_Unit = repmat({'percent'},numel(dataNomTempTime.Analyte),1);
                        dataOut = [dataOut; dataNomTempTime];
                    end
                end
            end
    end
end

