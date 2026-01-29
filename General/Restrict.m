function [samples, originalIndex, intervalID] = Restrict(samples,intervals,opt)

%Restrict - Keep only samples that fall in a given list of time intervals.
%
% Keep only samples (positions, spikes, LFP, etc.) that fall in a given list of
% time intervals.
%
% The remaining epochs can optionally be 'shifted' next to each other in time,
% removing the time gaps between them (which result from discarded samples),
% in which case they are also shifted globally to start at t = 0.
%
%  USAGE
%
%    [samples, originalIndex, intervalID] = Restrict(samples,intervals,<options>)
%
%    samples        <a href="matlab:help samples">samples</a> to restrict.
%                   They have to be provided in the first column.
%    intervals      list of sorted (start,stop) pairs, where intervals(:,1)
%                   is in ascending order
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'shift'       shift remaining epochs together in time (default = 'off')
%     'matrix'      if 'off' (default), when samples is a row vector it is treated
%                   as a column
%     'verbose'     if 'on' (default), display warning for empty inputs
%    =========================================================================
%
%  NOTE
%
%    For more advanced time restriction of samples, use <a href="matlab:help InIntervals">InIntervals</a>.
%
%  SEE
%
%    See also ConsolidateIntervals, SubtractIntervals, ExcludeIntervals,
%    IntersectIntervals, InIntervals, FindInInterval, CountInIntervals, PlotIntervals.

% Copyright (C) 2004-2011 by Michaël Zugaro, & (C) 2025 by Pietro Bozzo (vectorized algorithm)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    samples
    intervals (:,2)
    opt.shift {mustBeGeneralLogical} = false
    opt.matrix {mustBeGeneralLogical} = false
    opt.verbose {mustBeGeneralLogical} = true
end

% handle empty input
opt.verbose = GeneralLogical(opt.verbose);
if isempty(intervals)
    samples = []; originalIndex = []; intervalID = [];
    opt.verbose && fprintf(1,'Restriction over empty intervals.\n');
    return
end
if isempty(samples)
    samples = []; originalIndex = []; intervalID = [];
    opt.verbose && fprintf(1,'No samples to restrict.\n');
    return
end

samples = samples(~isnan(samples(:,1)),:);
intervals = double(intervals);
samples = double(samples);
opt.shift = GeneralLogical(opt.shift);
opt.matrix = GeneralLogical(opt.matrix);

transpose = false;
if ~opt.matrix && size(samples,1) == 1
    samples = samples(:);
    transpose = true;
end

% 1. Restrict

% consolidate intervals to use vectorized algorithm
intervals = ConsolidateIntervals(intervals);

% flatten and extend intervals by a small amount to compensate for vectorization, which treats intervals as right-open, i.e., [t1,t2)
deltaT = 10^5 * eps;
flattenedIntervals = intervals.';
interIntervalDistance = intervals(2:end,1) - intervals(1:end-1,2);
flattenedIntervals(2,[interIntervalDistance > deltaT;true]) = flattenedIntervals(2,[interIntervalDistance > deltaT;true]) + deltaT;
flattenedIntervals = flattenedIntervals(:);

% assign and index to each sample time, only odd indeces belong to intervals
ind = discretize(samples(:,1),flattenedIntervals);
status = mod(ind,2) == 1; % status(i) == 1 iff samples(i,1) falls in an interval

% get index of interval for every valid sample
if nargout > 2 || opt.shift
    intervalID = (ind(status) + 1) / 2;
end

% restrict samples
samples = samples(status,:);

% get index of samples to keep in original samples
if nargout > 1
    originalIndex = find(status);
end

% 2. Shift
if opt.shift && ~isempty(samples)
	% samples in each interval will be shifted next to end of the previous interval
	% let dti be the time difference between interval i and interval i+1, then interval n+1 must be shifted by sum_i from 1 to n of dti
	% 1) compute the cumulative shifts dt1, dt1+dt2, dt1+dt2+dt3, ... 
	cumulativeShift = [0;cumsum(interIntervalDistance)];
	% 2) assign cumulative shifts to samples
	shifts = cumulativeShift(intervalID);
	% 3) shift
	samples(:,1) = samples(:,1) - shifts - intervals(1,1);
end

if transpose, samples = samples'; end