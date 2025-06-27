function intersection = IntersectIntervals(a,b)

%IntersectIntervals - Find intersection of two sets of sorted intervals.
%
% Find intersection I of two sets of intervals A and B, e.g.,
%
%   A : [10,20] [25,35] [45,50],    B : [15,30] [55,70],    then
%   I : [10,35] [45,50] [55,70]
%
%  USAGE
%
%    intersection = IntersectIntervals(a,b)
%
%    a, b           lists of sorted non-overlapping (start,stop) pairs, i.e., 
%                   first column is in ascending order and intervals of the
%                   list are not overlapping (see ConsolidateIntervals to make
%                   an interval list non overlapping)
%
%  SEE
%
%    See also ConsolidateIntervals, SubtractIntervals, ExcludeIntervals,
%    InIntervals, Restrict, FindInInterval, CountInIntervals, PlotIntervals.

% Copyright (C) 2025 Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    a (:,2)
    b (:,2)
end

% transpose and flatten to use discretize
a = a.'; a = a(:);
b = b.'; b = b(:);

% validate input
if any(a(2:end) < a(1:end-1)) || any(b(2:end) < b(1:end-1))
    error('IntersectIntervals:inputFormat','Inputs must be sorted lists of non-overlapping intervals')
end

% ind(i) is odd iff a(i) falls in an interval of b
ind = discretize(a,b);
% NaN at the beginning means 0 (before first interval of b), at the end means numel(b) (after last interval of b)
nan_ind = isnan(ind);
nan_beginning = find(nan_ind(1:end-1) & ~nan_ind(2:end),1,'first');
ind(1:nan_beginning) = 0;
nan_end = find(~nan_ind(1:end-1) & nan_ind(2:end),1,'first');
ind(nan_end+1:end) = numel(b);

% if [ind(i),ind(i+1)] are even, it's an interval of a to exclude entirely
keep_ind = mod(ind(1:2:end),2) | mod(ind(2:2:end),2);
keep_ind = repelem(keep_ind,2);
ind = ind(keep_ind);

% odd_ind(i) is index of an odd element of ind, which will be replaces by a(odd_ind(i))
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