function variantTable = variant2table(variant)
% variant2table  Convert a SimBiology variant's content into a table.
%   variantTable = variant2table(variant) unpacks the {Type,Name,Property,Value}
%   entries of variant.content into a table with one row per content element.
%
%   Inputs:
%     variant      - SimBiology variant object.
%   Outputs:
%     variantTable - table of the variant's content.

for i=1:length(variant.content)
    ContentIndex(i,1) = i;
    Type{i,1} = variant.content{i}{1};
    Name{i,1} = variant.content{i}{2};
    Property{i,1} = variant.content{i}{3};
    Value(i,1) = variant.content{i}{4};
end
variantTable = table(ContentIndex,Type,Name,Property,Value);