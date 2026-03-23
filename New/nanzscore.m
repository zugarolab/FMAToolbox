function [z,mu,sigma] = nanzscore(x,w,dim)
%nanzscore - Compute standardized z-scores, ignoring NaN and Inf values.
%
%  USAGE
%
%    [z,mu,sigma] = nanzscore(x,flag,dim)
%
%    x              data to standardize
%    w              indicator for the standard deviation passed to std(x,w):
%                   - 0: normalize by N - 1 to compute the sample std (default)
%                   - 1: normalize by N to compute the population std
%                   - non-negative vector: weigths to compute weighted std
%    dim            dimension along which to operate (default: first non-
%                   singleton dimension)
%
%  OUTPUT
%
%    z              z-scores, has the same dimension as 'x'
%    mu             mean of 'x' along 'dim'
%    sigma          standard deviation of 'x' along 'dim'
%
%  EXAMPLES
%
%    % z-score a vector, ignoring NaN
%    x = [1 2 NaN 4 5];
%    z = nanzscore(x);               
%
%    % row-wise z-score of a matrix
%    x = magic(4);
%    x(2,3) = NaN;
%    [z,mu,sigma] = nanzscore(x, 0, 2)
%
%  SEE
%
%    See also sem, nansem, semedian.

arguments
    x
    w = []
    dim = []
end

% handle empty input
if isempty(x)
  [z,mu,sigma] = deal([]);
  return
end
x(isinf(x)) = NaN;

% default 'dim' value
if isempty(dim)
    dim = find(size(x) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
end

% z-score
try
    mu = mean(x,dim,'omitnan');
    sigma = std(x,w,dim,'omitnan');
catch ME
    throw(ME)
end
sigma0 = sigma;
sigma0(sigma0<eps) = 1; % avoid division by 0
z = bsxfun(@minus,x,mu);
z = bsxfun(@rdivide,z,sigma0);