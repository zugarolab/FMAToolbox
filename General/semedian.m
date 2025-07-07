function s = semedian(x,dim)

%semedian - Compute standard error of the median, ignoring NaNs.
%
%  USAGE
%
%    s = semedian(x)
%
%    x              vector or matrix over which the error should be computed
%    dim            optionally, dimension along which to operate, default
%                   is first dimension of x different than 1

% Copyright (C) 2013 by Michaël Zugaro & (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
  x (:,:)
  dim (1,1) {mustBeNumeric,mustBeInteger} = 0
end

% assign default value
if dim == 0
    dim = find(size(x) ~= 1,1);
end

n = size(x,dim);
m = median(x,dim,'omitnan');
s = sqrt( sum((x-m).^2,dim,'omitnan') / (n*(n-1)) );