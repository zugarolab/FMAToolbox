function [consolidated,target] = ConsolidateIntervals(intervals,opt)

%ConsolidateIntervals - Consolidate intervals.
%
% Consolidate overlapping intervals, e.g. replace [10,20] [15,25] with [10,25].
%
%  USAGE
%
%    [consolidated,target] = ConsolidateIntervals(intervals,<options>)
%
%    intervals      list of intervals
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'strict'      intervals with common bounds are consolidated ('off')
%                   or kept separate ('on') (default = 'off')
%     'epsilon'     intervals with close enough bounds (distance lesser than
%                   epsilon) are also consolidated (default = 0)
%    =========================================================================
%
%  OUTPUT
%
%    consolidated   consolidated intervals
%    target         for each original interval, the index of the consolidated
%                   interval to which it belongs (empty intervals yield NaN)
%
%  SEE
%
%    See also SubtractIntervals, ExcludeIntervals, IntersectIntervals, InIntervals,
%    Restrict, FindInInterval, CountInIntervals, PlotIntervals.

% Copyright (C) 2004-2011 by Michaël Zugaro, 2026 by Pietro Bozzo (vectorized algorithm)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    intervals (:,2)
    opt.strict {mustBeGeneralLogical} = false
    opt.epsilon (1,1) {mustBeGreaterThanOrEqual(opt.epsilon,0)} = 0
end

if isempty(intervals)
    [consolidated,target] = deal(intervals);
    return
end

% default values
opt.strict = GeneralLogical(opt.strict);
if ~opt.strict && ~opt.epsilon
    opt.epsilon = 1e-10;
end

% widen intervals
if opt.epsilon
    intervals = intervals + [-1,1]*opt.epsilon;
end

% sort by start time, remove empty intervals
is_empty = intervals(:,1) > intervals(:,2) | any(isnan(intervals),2); % remember empty intervals for target
[~,empty_sort_ind] = sortrows(intervals);
intervals = intervals(~is_empty,:);
[intervals,orig_sort_ind] = sortrows(intervals);
if isempty(intervals)
    consolidated = intervals;
    target = nan(size(is_empty));
    return
end

% flatten and argsort to find overlaps
intervals = intervals.';
intervals = intervals(:);
[~,ind] = sort(intervals);

% remove all ind which are followed by at least one smaller element
m = min(ind(end-1:end));
is_ok = true(size(ind));
is_ok(end-1) = ind(end-1) < ind(end);
for i = 2 : numel(ind)-1
    is_ok(end-i) = ind(end-i) < m;
    m = min(ind(end-i),m);
end
ind = ind(is_ok);

% remove consecutive even elements (they all represent interval ends)
is_even = ~mod(ind,2);
is_not_consecutive = [~is_even(1:end-1) | ~is_even(2:end); true];
ind = ind(is_not_consecutive);

% rebuild intervals
consolidated = reshape(intervals(ind),2,[]).';

% re-shorten intervals
if opt.epsilon
    consolidated = consolidated + [1,-1]*opt.epsilon;
end

% target
if nargout > 1
    is_ok(is_ok) = is_not_consecutive;
    is_ok = is_ok(1:2:end);
    target = cumsum(is_ok);
    % resort to match original order
    [~,inv_sort_ind] = sort(orig_sort_ind);
    sorted_target = nan(size(is_empty)); % empty intervals belong to NaN
    sorted_target(~is_empty(empty_sort_ind)) = target(inv_sort_ind);
    target = sorted_target;
end