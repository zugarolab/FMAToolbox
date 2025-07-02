function s = nansem(x,dim)

%nansem - Compute standard error of the mean (SEM), ignoring NaNs.
%
%  USAGE
%
%    s = nansem(x)
%
%    x              vector or matrix over which the sem should be computed
%    dim            optionally, dimension along which to operate, default
%                   is first dimension of x different than 1

% Copyright (C) 2008-2022 by Michaël Zugaro, Ralitsa Todorova & (C) 2025 by Pietro Bozzo
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

n = sum(~isnan(x),dim);
s = std(x,0,dim,'omitmissing') ./ sqrt(n);