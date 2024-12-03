function c = WeightedCorr(weights,x,y)

% Provide a matrix of weights, and this function will check the 
% correlation between the X and Y dimensions of the matrix.
% You can provide the X-values and the Y-values for each column
% and row.

if exist('x','var') && ~isempty('x')
    if isvector(x), x = repmat(x(:)',size(weights,1),1); end
else
    [x,~] = meshgrid(1:size(weights,2),1:size(weights,1));
end
if exist('y','var') && ~isempty('y')
    if isvector(y), y = repmat(y(:),1,size(weights,2)); end
else
    [~,y] = meshgrid(1:size(weights,2),1:size(weights,1));
end

x = x(:);
y = y(:);
w = weights(:);

mX = nansum(w.*x)./nansum(w);
mY = nansum(w.*y)./nansum(w);


covXY = nansum(w.*(x-mX).*(y-mY))./nansum(w(:));
covXX = nansum(w.*(x-mX).*(x-mX))./nansum(w(:));
covYY = nansum(w.*(y-mY).*(y-mY))./nansum(w(:));

c = covXY ./ sqrt(covXX.*covYY);

