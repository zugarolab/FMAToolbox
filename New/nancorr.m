function [c,p] = nancorr(x,y,varargin);

bad = any(isnan([x,y]),2);
[c,p] = corr(x(~bad,:),y(~bad,:),varargin{:});
end