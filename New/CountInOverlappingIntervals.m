function count = CountInOverlappingIntervals(spikes,intervals)

% If you call CountInIntervals with overlapping intervals, a spike already
% counted in the first interval will not be counted again in the next
% (overlapping) interval.
% To prevent this behavior and count the same spike in as many intervals as
% it appears, this helper function exists. 
%
% Copyright (C) 2026 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


m = cummax(intervals(:,1));
bad = intervals(:,1)<m;

if ~any(bad)
    count = CountInIntervals(spikes,intervals);
    return
else
    count = zeros(size(intervals(:,1)));
    count(~bad) = CountInIntervals(spikes,intervals(~bad,:));
    count(bad) = CountInOverlappingIntervals(spikes,intervals(bad,:));
end

