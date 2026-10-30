function dataStat = grpGeoStats(data,grpVar,dataColName)
% grpGeoStats  Geometric mean and SD of a column, grouped by a variable.
%   dataStat = grpGeoStats(data,grpVar,dataColName) computes, for each level of
%   data.(grpVar), the geometric mean and geometric standard deviation of
%   data.(dataColName) and returns them in a table keyed by the grouping value.
%
%   Inputs:
%     data        - data table.
%     grpVar      - name of the grouping column.
%     dataColName - name of the column to summarise.
%   Outputs:
%     dataStat    - table of geomean_/geostd_ columns plus the grouping value.

% Find the unique values for the grouping variable
uniqueGrpVar = unique(data.(grpVar));
% Initiate output
geoMean = zeros(length(uniqueGrpVar),1);
geoStd = zeros(length(uniqueGrpVar),1);

for i=1:length(uniqueGrpVar)
    dataTemp = data(data.(grpVar)==uniqueGrpVar(i),:);
    geoMean(i,1) = geomean(dataTemp.(dataColName));
    geoStd(i,1) = geostd(dataTemp.(dataColName));
end
dataStat = table(geoMean,geoStd);
dataStat.(grpVar) = uniqueGrpVar;
dataStat = renamevars(dataStat,["geoMean","geoStd"],{['geomean_',dataColName],['geostd_',dataColName]});
