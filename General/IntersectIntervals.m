function varargout = IntersectIntervals(a)

%IntersectIntervals - Find intersection between sets of sorted intervals.
%
% Find intersection I of sets of intervals A, B, C, ...; e.g.,
%
%   A : [10,20] [25,35] [45,50],    B : [15,30] [55,70],    then
%   I : [15,20] [25,30]
%
%  USAGE
%
%    [intersection,Ia,Ib,Ic,...] = IntersectIntervals(a,b,c,...)
%
%    a, b, c, ...   any number of lists of (start,stop) pairs, intervals containing
%                   NaNs are ignored
%
%  OUTPUT
%
%    intersection   intersection intervals
%    Ii             for each intersection interval, the index of the original
%                   interval in i which contains it
%                   NOTE: for now, ONLY Ia is meaningful, computation of all other 
%                   indeces needs to be correctly implemented
%
%  SEE
%
%    See also ConsolidateIntervals, SubtractIntervals, ExcludeIntervals,
%    InIntervals, Restrict, FindInInterval, CountInIntervals, PlotIntervals.

% Copyright (C) by 2025 Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments (Repeating)
    a (:,2)
end

% the algorithm exploits the fact that, when both interval lists are flattened, elements of a fall in even / odd
% intervals of b depending on whether they would fall outside / inside original intervals of b
% (this requires them to be sorted and non overlapping)

if nargout > numel(a) + 1
    error('IntersectIntervals:nOuts','Too many indeces requested')
end

% return single interval unchanged
if isscalar(a)
    varargout{1} = a{1}; varargout{2} = (1 : size(a{1},1)).';
    return
end

% 1: more than 2 intervals, recursive call
indeces = {};
if numel(a) > 2
    [out{1:numel(a)}] = IntersectIntervals(a{2:end});
    a = [a(1),out(1)];
    indeces = out(3:end);
end

% 2: intersect 2 interval sets

% consolidate and flatten to use discretize
a = cellfun(@(x) ConsolidateIntervals(x).', a,'UniformOutput',false);
b = a{2}(:);
a = a{1}(:);
a_orig = a;

if isempty(b)
    [varargout{1:3}] = deal([]);
    return
end

% ind(i) is odd iff a(i) falls in an interval of b
ind = discretize(a,b);
% a before first interval of b means 0, a after last interval of b means numel(b)
ind(a<b(1)) = 0;
ind(a>b(end)) = numel(b);

% handle case when a ∋ [10,20] and b ∋ [20,35] : interval [20,20] must not be in result
end_nz = ind ~= 0 & (1 - mod(1:numel(ind),2).'); % end_nz(i) is 1 iff ind(i) is not 0 and i is even
change_ind = a(end_nz) == b(ind(end_nz)); % change_ind(j) is 1 iff corresponding interval must be shortened
find_end_nz = find(end_nz);
ind(find_end_nz(change_ind)) = ind(find_end_nz(change_ind)) - 1;

% if [ind(i),ind(i+1)] are equal and even, it's an interval of a to exclude entirely
keep_ind = mod(ind(1:2:end),2) | (ind(1:2:end) ~= ind(2:2:end));
keep_ind = repelem(keep_ind,2);
ind = ind(keep_ind);
a = a(keep_ind);
if isempty(ind)
    [varargout{1:3}] = deal([]);
    return
end

% odd_ind(i) is index of an odd element of ind, which will be replaced by a(odd_ind(i))
odd_ind = find(mod(ind,2));
% round start to upper odd number and stop to upper even number
start = ind(1:2:end) + 1 - mod(ind(1:2:end),2);
stop = ind(2:2:end) + mod(ind(2:2:end),2);
% expand each [start(i),stop(i)] pair to include intervals in between
new_ind = linspaceVector(start,stop);
% remap odd_ind to point elements of new_ind, which will have to be drawn from a
remapping = [ones(size(start)),stop-start].';
remapping = cumsum(remapping(:));
new_odd_ind = remapping(odd_ind);

% initialize intervals as [b(new_ind(i)),b(new_ind(i+1))]
intersection = b(new_ind);
% when ind(j) was odd, replace corresponding element with a(odd_ind(j))
intersection(new_odd_ind) = a(odd_ind);
% return as interval list
intersection = [intersection(1:2:end),intersection(2:2:end)];

% output
varargout{1} = intersection;
if nargout > 1
    varargout{2} = ceil(discretize(intersection(:,1),a_orig)/2); % Ia
end
if nargout > 2
    varargout{3} = ceil(new_ind(1:2:end)/2); % Ib
end
if nargout > 3
    varargout = [varargout,indeces]; % Ic, Id, ...
end