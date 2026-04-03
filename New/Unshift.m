function timestamps = Unshift(timestamps,intervals)

%Unshift - The opposite of option 'shift' in <a href="matlab:help Restrict">Restrict</a>.
%
%  USAGE
%
%    timestamps = Unshift(timestamps,intervals)
%
%    timestamps     time stamps to unshift, must be provided in the first column
%    intervals      each row is [start, stop] time (s) of an interval,
%                   previously used to shift 'timestamps'
%
%  NOTES
%
%    Operates only on the first column of 'timestamps', as <a href="matlab:help InIntervals">InIntervals</a>.

% Copyright (C) 2018 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

shiftedIntervals(:,2) = cumsum(diff(intervals,[],2));
shiftedIntervals(:,1) = [0; shiftedIntervals(1:end-1,2)];
toAdd = intervals(:,1) - shiftedIntervals(:,1);
[ok,w] = InIntervals(timestamps,shiftedIntervals);
timestamps(~ok) = nan;
timestamps(ok) = timestamps(ok) + toAdd(w(ok));