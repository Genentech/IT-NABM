function dataOut = add_normData(data,options)
% add_normData  Append baseline-normalised (% of baseline) rows to a dataset.
%   dataOut = add_normData(data, Group=..., Matrix=..., Analyte=...) computes,
%   for each selected Group/Matrix/Analyte and subject ID, the baseline mean
%   (Time<=0) and appends a copy of the series normalised to 100% of baseline,
%   with the analyte name suffixed '_norm' and the unit set to 'percent'.

% Sample call: norm_Data(data,Group=[1001,1002,1003], Matrix={'CSF','Plasma'}, Analyte={'C1s','aC1s','C1sComplex'});

   

arguments
    data table  
    options.Group double = []
    options.ID double = []
    options.Analyte string = []
    options.Matrix string = []
end

dataOut = data;

for i=1:numel(options.Group)
    for j=1:numel(options.Matrix)
        for k=1:numel(options.Analyte)
            dataTemp = getData(data, Group=options.Group(i), Matrix=options.Matrix{j}, Analyte=options.Analyte{k});
            uniqueID = unique(dataTemp.ID);
            for p=1:numel(uniqueID)
                dataTempNorm = dataTemp(dataTemp.ID==uniqueID(p),:);
                baseline = mean(dataTempNorm.Conc(find(dataTempNorm.Time<=0)),"omitnan");
                dataTempNorm.Conc = dataTempNorm.Conc/baseline*100;
                dataTempNorm.Analyte = strcat(dataTempNorm.Analyte,repmat({'_norm'},numel(dataTempNorm.Analyte),1));
                dataTempNorm.Conc_Unit = repmat({'percent'},numel(dataTempNorm.Analyte),1);
                dataOut = [dataOut; dataTempNorm];
            end
        end
    end
end
