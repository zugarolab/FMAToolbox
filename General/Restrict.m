function [samples, originalIndex, intervalID] = Restrict(samples,intervals,varargin)

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
%    samples         <a href="matlab:help samples">samples</a> to restrict.
%                   They have to be provided in the first column.
%    intervals      list of (start,stop) pairs
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'shift'       shift remaining epochs together in time (default = 'off')
%    =========================================================================
%
%  NOTE
%
%    For more advanced time restriction of samples, use <a href="matlab:help InIntervals">InIntervals</a>.
%
%  SEE
%
%    See also ConsolidateIntervals, SubtractIntervals, ExcludeIntervals,
%    InIntervals, Restrict, FindInInterval, CountInIntervals, PlotIntervals.

% Copyright (C) 2004-2011 by Michaël Zugaro, 2025 Pietro Bozzo (vectorized algorithm)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
shift = 'off';
transpose = false;
try samples(isnan(samples(:,1)),:) = []; end

% Check number of parameters
if nargin < 2 || mod(length(varargin),2) ~= 0
  error('Incorrect number of parameters (type ''help <a href="matlab:help Restrict">Restrict</a>'' for details).');
end

% Check parameters
intervals = double(intervals);
samples = double(samples);
if ~isdmatrix(intervals) || size(intervals,2) ~= 2
  error('Incorrect intervals (type ''help <a href="matlab:help Restrict">Restrict</a>'' for details).');
end

if size(samples,1) == 1
	samples = samples(:);
    transpose = true;
end

% Parse parameter list
for i = 1:2:length(varargin)
	if ~ischar(varargin{i})
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Restrict">Restrict</a>'' for details).']);
	end
	switch(lower(varargin{i}))
		case 'shift'
			shift = varargin{i+1};
            if ~isastring(shift,'on','off')
                error('Incorrect value for property ''shift'' (type ''help <a href="matlab:help Restrict">Restrict</a>'' for details).');
            end
        otherwise
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Restrict">Restrict</a>'' for details).']);
	end
end

if isempty(samples)
    samples = [];
    disp('No samples to restrict.');
    return
end

% Restrict
% prepare flattened interval array for faster vectorized algorithm
flattenedIntervals = intervals.';
% extend intervals by a small amount to compensate for vectorization, which
% treats intervals as right-open, i.e., [t1,t2)
deltaT = 10^5 * eps;
interIntervalDistance = intervals(2:end,1) - intervals(1:end-1,2);
flattenedIntervals(2,[interIntervalDistance > deltaT;true]) = flattenedIntervals(2,[interIntervalDistance > deltaT;true]) + deltaT;
flattenedIntervals = flattenedIntervals(:);
% check if faster algorithm can be used
if all(flattenedIntervals(2:end) >= flattenedIntervals(1:end-1))
    % assign and index to each sample time, only odd indeces belong to intervals
    ind = discretize(samples(:,1),flattenedIntervals);
    % identify samples to keep; status(i) == 1 iff samples(i,1) falls in an interval
    status = mod(ind,2) == 1;
    % get index of interval for every valid sample
    if nargout > 2 || strcmp(shift,'on')
        intervalID = (ind(status) + 1) / 2;
    end
else
    % original algorithm
    [status,interval] = InIntervals(samples,intervals);
    intervalID = interval(status);
end
% restrict samples
samples = samples(status,:);
% get index of samples to keep in original samples
if nargout > 1
    originalIndex = find(status);
end

% Shift?
if strcmp(shift,'on') && ~isempty(samples)
	% Samples in each interval will be shifted next to end of the previous interval
	% Let dti be the time difference between interval i and interval i+1,
    % then interval n+1 must be shifted by sum_i from 1 to n of dti
	% 1) Compute the cumulative shifts dt1, dt1+dt2, dt1+dt2+dt3, ... 
	cumulativeShift = [0;cumsum(interIntervalDistance)];
	% 2) Assign cumulative shifts to samples
	shifts = cumulativeShift(intervalID);
	% 3) Shift
	samples(:,1) = samples(:,1) - shifts - intervals(1,1);
end

if transpose, samples = samples'; end