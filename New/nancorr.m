function [c,p] = nancorr(x,y,varargin)

[c,p] = corr(x,y,'rows','pairwise',varargin{:});

end