function groups_reordered = groupOrder(args,data)
% groupOrder  Order the fitted groups by their first appearance in the data.
%   groups_reordered = groupOrder(args,data) returns the groups listed in
%   args.groupFit, ordered as they first appear in data.Group.
%
%   Inputs:
%     args - struct with the groupFit list.
%     data - dataset table with a Group column.
%   Outputs:
%     groups_reordered - ordered vector of group IDs.
[uniqueGroups, first_indices] = unique(data.Group, 'stable');
groupOrderIndices = find(ismember(uniqueGroups,args.groupFit));
groups_reordered = uniqueGroups(groupOrderIndices);
