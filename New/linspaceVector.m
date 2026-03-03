function y = linspaceVector(start,stop,n)

%linspaceVector - Produce linearly spaced vector with multiple [start stop] points.
%
% Shorthand for a loop calling the colon operator for each [start stop] pair:
%
% y = [];
% for i = 1 : length(start) 
%     y = [y;(start(i):n(i):stop(i)).'];
% end
%
%  USAGE
%
%    y = linspaceVector(start,stop,n)
%
%    start, stop    start and stop pairs
%    n              increment for each [start, stop] pair (default is 1)
%    
%  EXAMPLE
%
%    >> linspaceVector([2;5;23],[2;8;23])
%    [2;5;6;7;8;23]

% Copyright (C) 2019 by Ralitsa Todorova & (C) 2026 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    start (:,1)
    stop (:,1)
    n (:,1) = []
end

if isempty(n)
    n = ones(size(start));
elseif isscalar(n)
    n = repelem(n,numel(start),1);
end
if numel(start) ~= numel(stop) || numel(start) ~= numel(n) || numel(stop) ~= numel(n)
    error('linspaceVector:inputSize','Arguments ''start'', ''stop'', and ''n'' must have the same number of elements')
end

y = cell(size(start));
for i = 1 : numel(start)
    y{i} = start(i) : n(i) : stop(i);
end
y = [y{:}].';