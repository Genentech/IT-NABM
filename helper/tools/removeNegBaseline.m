% remove the negative datapoints 
function dataWide = removeNegBaseline(dataWide)
% removeNegBaseline  Collapse pre-dose (Time<=0) rows into a single Time=0 baseline.
%   dataWide = removeNegBaseline(dataWide) drops all Time<=0 rows and, for each
%   Group/ID, replaces them with one averaged baseline row at Time=0 (SimBiology
%   does not accept negative times).
%
%   Inputs/Outputs:
%     dataWide - wide data table, returned with baseline rows collapsed.

% remove all the Time<=0 data as the baseline
dataWideNew = dataWide(dataWide.Time>0,:);

% add the units
groups = unique(dataWide.Group);
for i=groups'
    dataWideGroup = dataWide(dataWide.Group==i,:);
    IDs = unique(dataWideGroup.ID);
    for j=IDs'
        dataWideGroupID = dataWideGroup(dataWideGroup.ID==j,:);
        dataWideGroupIDNeg = dataWideGroupID(dataWideGroupID.Time<=0,:);
        if ~isempty(dataWideGroupIDNeg)
            dataTemp = mean(dataWideGroupIDNeg,'omitnan');
            dataTemp.Time = 0;
            dataWideNew = [dataWideNew; dataTemp];
        end
    end
end

dataWide = dataWideNew;
