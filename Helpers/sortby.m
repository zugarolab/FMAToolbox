function output = sortby(data,vector2order)

[~,order] = sort(vector2order,'MissingPlacement','last');
output = data(order,:);