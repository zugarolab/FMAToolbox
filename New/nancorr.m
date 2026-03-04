function [c,p] = nancorr(x,y,varargin)

if ~exist('y','var'), y = x(:,2); x = x(:,1); end

[c,p] = corr(x,y,'rows','pairwise',varargin{:});

end