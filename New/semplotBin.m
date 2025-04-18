function [handle,matrix,u] = semplotBin(x,y,nBins,varargin)

if size(y,1)==1, y = y'; end
if size(y,2)>1, y = y(:); end
if size(x,1)==1, x = x'; end
if size(x,2)>1, x = x(:); end
if numel(x) ~= numel(y), error('x and y should be vectors of equal length'); end

nans = isnan(x);
x(nans) = []; y(nans) = [];
if nBins==0
    [binned,~,b] = unique(x);
else
    bins = linspace(min(x),max(x),nBins)';
    b = FindClosest(bins,x);
    binned = bins(b);
end
matrix = nan(size(y,1),max(b));
matrix(sub2ind(size(matrix),(1:size(y,1))',b)) = y;
matrix = sort(matrix);
matrix(sum(~isnan(matrix),2)==0,:) = [];

u = unique(binned);
if length(u)<size(matrix,2)
    u = interp1(unique(b),u,(1:size(matrix,2))');
end

handle = semplot(u,matrix,varargin{:});