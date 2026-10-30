% reorders the wide data columns 
function dataWide = reorderDataWide(dataWide,VariableOrdersWanted)
% reorderDataWide  Reorder wide-table columns to a preferred ordering.
%   dataWide = reorderDataWide(dataWide,VariableOrdersWanted) moves the columns
%   named in VariableOrdersWanted to the front (in that order) and appends the
%   remaining columns in their original order; warns if the reorder is incomplete.
%
%   Inputs:
%     dataWide             - wide data table.
%     VariableOrdersWanted - cell array of column names to place first.
%   Outputs:
%     dataWide             - reordered table.

% add the units
VariableOrderOld = dataWide.Properties.VariableNames;
newIndex = zeros(1,length(VariableOrdersWanted));
counterNew=1; 
counterHold=1; 
for i = 1:length(VariableOrderOld)
    indexTemp = find(ismember(VariableOrdersWanted,VariableOrderOld{i}));
    if ~isempty(indexTemp)
        newIndex(indexTemp) = i;
        counterNew = counterNew+1;
    else
        holdIndex(counterHold)=i;
        counterHold = counterHold+1;
    end
end
newIndex = [newIndex,holdIndex];
if (length(newIndex)~=length(VariableOrderOld)) || (~isempty(find(newIndex==0)))
    display("Warning: Something went wrong in reordering of the data columns");
end
dataWide = dataWide(:,newIndex);
