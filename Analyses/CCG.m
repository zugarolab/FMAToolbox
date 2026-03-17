function [ccg,t,tau,c] = CCG(times,id,opt)

%CCG - Compute multiple cross- and auto-correlograms, or cross-covariances
%
%  USAGE
%
%    [ccg,t,tau,c] = CCG(times,id,<options>)
%
%    times          times of all events in seconds (see NOTEs); can also be a
%                   cell array where each cell contains event times for an ID
%    id             ID for each event (e.g. unit ID) from 1 to n, optional (default:
%                   ones(size(times)), ignored if 'times' is a cell array)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration in s of each xcorrelogram (default = 2)
%     'nBins'       number of time bins (default = duration/binSize)
%     'smooth'      smoothing size in bins (default = 0, no smoothing)
%     'groups'      group number (1 or 2) for each event, to restrict
%                   cross-correlograms to pairs across two groups of IDs
%                   (see EXAMPLE #2 below)
%     'range'       range of IDs corresponding to 2nd and 3rd dimensions of
%                   'ccg', 2 values if 'groups' is empty (defaul: [1,max(id)]),
%                   4 values otherise (defaul: min and max of IDs of the
%                   two groups)
%     'mode'        'ccg' or 'ccv' (default = 'ccg')
%     'alpha'       significance level to determine correlated pairs (only
%                   for 'ccv' mode)
%     'totalTime'   recording duration in s (default = 0 corresponds to
%                   max(times) - min(times), only for 'ccv' mode)
%     'Fs'          this property is deprecated and has no effect
%     'norm'        normalization of the CCG, either:
%                   - 'counts', CCG contains raw event counts (default)
%                   - 'rate', CCG is in Hz (event / s)
%    =========================================================================
%
%  NOTES
%
%    The size of the cross-correlograms can be supplied either as a duration
%    (property 'duration') or as an number of bins (property 'nBins').
%    Parameters 'times', 'id' and 'group' can be obtained using <a href="matlab:help CCGParameters">CCGParameters</a>.
%    As a special case, when computing the correlograms of spike trains, one
%    can use the output of <a href="matlab:help GetSpikes">GetSpikes</a> either directly or in combination with
%    <a href="matlab:help CCGParameters">CCGParameters</a>. See EXAMPLES below.
%
%  OUTPUT
%      ccg          value of cross-correlograms or cross-covariances
%                   dimensions are (nbins,m,n) where m is the number of
%                   reference time series (e.g. reference units) and n the
%                   number of referenced time series (in general m = n,
%                   except when using option 'groups')
%      t            time bins
%      tau          lag times for a each pair (mode 'ccv' only)
%      c            maximum cross-covariance for a each pair (mode 'ccv' only)
%
%  EXAMPLES
%
%    % Auto- and cross-correlograms between all neurons
%    spikes = GetSpikes('output','numbered');
%    [ccg,t] = CCG(spikes(:,1),spikes(:,2));
%
%    % Only tetrode #1 vs tetrode #2 (e.g. mPFC vs HPC neurons)
%    pfc = GetSpikes([1 -1],'output','numbered');
%    hpc = GetSpikes([2 -1],'output','numbered');
%    [s,ids,groups] = CCGParameters(pfc,hpc,2);
%    [ccg,t] = CCG(s,ids,'groups',groups);
%
%    % Between stimulations and MUA spikes
%    spikes = GetSpikes;
%    stimulatios = GetEvents('Stimulation');
%    d = [spikes(:,1) ones(size(spikes,1)) ; stimulations 2*ones(size(stimulations,1))];
%    d = sortrows(d);
%    [ccg,t] = CCG(d(:,1),d(:,2));
%
%    % To compute cross-covariances
%    [ccv,t,tau,C] = CCG(times,ids,'mode','ccv');
%
%  SEE
%
%    See also CCGParameters, ShortTimeCCG.

% Copyright (C) 2012-2013 by Michaël Zugaro, Marie Goutierre,
%           (C) 2017 by Dan Levenstein (norm and cell options)
%           (C) 2026 by Pietro Bozzo (optimization)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    times
    id (:,1) {mustBeInteger,mustBePositive} = []
    opt.binSize {mustBeScalarOrEmpty,mustBePositive} = []
    opt.duration {mustBeScalarOrEmpty,mustBePositive} = []
    opt.nBins {mustBeScalarOrEmpty,mustBeInteger,mustBePositive} = []
    opt.smooth (1,1) {mustBeNonnegative} = 0
    opt.groups (:,1) {mustBeMember(opt.groups,[1,2])} = []
    opt.range (1,:) {mustBeInteger,mustBePositive} = []
    opt.mode (1,1) string {mustBeMember(opt.mode,["ccg","ccv"])} = "ccg"
    opt.alpha (1,1) {mustBePositive} = 0.05
    opt.totalTime (1,1) {mustBeNonnegative} = 0
    opt.Fs = []
    opt.norm (1,1) string {mustBeMember(opt.norm,["counts","rate"])} = "counts"
end

% Validate input
if iscell(times)
    % convert to double array
    id = repelem((1:numel(times)).',cellfun(@numel, times),1);
    times = vertcat(times{:});
end
if ~isnumeric(times) || ~isvector(times)
    error('Parameter ''times'' must be a vector or cell (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if size(times,1) == 1
    times = times.';
end
if isempty(id) || isscalar(id)
    id = ones(size(times));
    nIDs = 1;
else
    if length(times) ~= length(id)
        error('Parameters ''times'' and ''id'' have different lengths (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
    end
    nIDs = max(id);
end
if ~isempty(opt.groups) 
    if length(times) ~= length(opt.groups)
        error('Property ''groups'' must have the same length as ''times'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
    end
    [id1{1:2}] = bounds(id(opt.groups == 1));
    id1 = [id1{:}];
    [id2{1:2}] = bounds(id(opt.groups == 2));
    id2 = [id2{:}];
end
if ~isempty(opt.range)
    opt.range = reshape(opt.range,2,[]).';
    range_min = min(opt.range,[],'all');
    range_max = max(opt.range,[],'all');
    n_range = range_max - range_min + 1;
    if any(diff(opt.range,1,2) <= 0)
        error('Property ''range'' must specfify increasing ranges. Type ''help <a href="matlab:help CCG">CCG</a>'' for details.')
    end
    if isempty(opt.groups)
        if numel(opt.range) ~= 2
            error('Incorrect value for property ''range'' (expected two elements). Type ''help <a href="matlab:help CCG">CCG</a>'' for details.')
        end
    else
        if numel(opt.range) ~= 4
            error('Incorrect value for property ''range'' (expected four elements). Type ''help <a href="matlab:help CCG">CCG</a>'' for details.')
        end
    end
end

% Default values
is_given = ~cellfun(@isempty, {opt.binSize,opt.duration,opt.nBins});
defaults = [0.01,2]; % bin size, duration
switch sum(is_given)
    % case 0: keep default values
    case 1
        % if binSize or duration is given, keep default value for the other
        % if nBins is given, compute duration
        if is_given(3)
            defaults(2) = defaults(1) * opt.nBins;
        end
    case 2
        if ~is_given(1)
            defaults(1) = opt.duration / opt.nBins;
        elseif ~is_given(2)
            defaults(2) = opt.binSize * opt.nBins;
        end
    case 3
        if opt.duration ~= opt.binSize * opt.nBins
            error('Incompatible ''duration'' and ''nBins'' parameters (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
        end
end
if isempty(opt.binSize), opt.binSize = defaults(1); end
if isempty(opt.duration), opt.duration = defaults(2); end
halfBins = round(opt.duration/opt.binSize/2);
opt.nBins = 2*halfBins + 1;
t = (-halfBins:halfBins).' * opt.binSize;
if opt.smooth == 0
    opt.smooth = 0.2; % corresponding to no smoothing
end
tau = [];
c = [];
if length(times) <= 1
    ccg = zeros(opt.nBins,nIDs,nIDs);
    return
end
if opt.totalTime == 0
    opt.totalTime = max(times) - min(times);
end

% Restrict events to range
if isempty(opt.range)
    if isempty(opt.groups)
        opt.range = [1,nIDs];
        n_range = nIDs;
    else
        opt.range = [id1;id2];
        range_min = min(opt.range,[],'all');
        n_range = nIDs - range_min + 1;
        id = id - range_min + 1;
        opt.range = opt.range - range_min + 1;
        nIDs = nIDs - range_min + 1;
    end
else
    [id,idx] = Restrict(id,[range_min,range_max]);
    times = times(idx);
    id = id - range_min + 1; % remap events in [1,range_max]
    opt.range = opt.range - range_min + 1;
    nIDs = min(nIDs,range_max) - range_min + 1;
    if ~isempty(opt.groups)
        opt.groups = opt.groups(idx);
    end
end

% Sort events in time
[times,idx] = sort(times);
id = id(idx);
if ~isempty(opt.groups)
    opt.groups = opt.groups(idx);
end
% add dummy spike to ensure correct number of IDs in CCGEngine output
times(end+1) = times(end) + 2*opt.duration;
id(end+1) = n_range;

% 1. CCG
counts = double(CCGEngine(times,id,opt.binSize,halfBins));
counts = reshape(counts,[opt.nBins,n_range,n_range]); % reshape to (time bins, IDs, IDs)
% rate normalization: counts / n_spikes / dt to convert to Hz
if opt.norm == "rate"
    for i = 1 : nIDs
        numREFspikes = sum(id==i); % number of reference events for group
        counts(:,i,:) = counts(:,i,:) ./ numREFspikes ./ opt.binSize;
    end
end
counts = smoothdata(flipud(counts),1,'gaussian',5*opt.smooth); % factor 5 chosen to match previous behavior of Smooth
if isempty(opt.groups)
    ccg = counts;
    % copy upper triangular part to lower to enforce ccg_ij(t) = ccg_ji(-t);
    for i = 1 : nIDs
       for j = 1 : i-1
           ccg(:,i,j) = ccg(end:-1:1,j,i);
       end
    end
else
    ccg = counts(:,opt.range(1,1):opt.range(1,2),opt.range(2,1):opt.range(2,2));
end

% 2. CCV
if strcmp(opt.mode,'ccv')

	% mean event rate for each ID
	eventRate = zeros(nIDs,1);
	for i = 1 : nIDs
		eventRate(i) = sum(id==i) / opt.totalTime;
	end

	% standardized cross-covariances
	ccv = zeros(size(ccg));
	tau = zeros(size(ccg,2),size(ccg,3));
	c = zeros(size(ccg,2),size(ccg,3));

	nPairs = size(ccg,2) * size(ccg,3);
	disp(['# pairs: ' int2str(nPairs)]);

	threshold = sqrt(2) * erfinv(1-(opt.alpha/length(t)));

	for i = 1:size(ccg,2)
		for j = 1:size(ccg,3)
		
			% Compute and normalize CCVs from CCGs
			if ~isempty(opt.groups)
				rate = eventRate(opt.range(1,1)+i-1)*eventRate(opt.range(2,1)+j-1);
			else
				rate = eventRate(i)*eventRate(j);
			end
			ccv(:,i,j) = sqrt((opt.binSize*opt.totalTime)/rate) * (ccg(:,i,j)/(opt.binSize*opt.totalTime)-rate);

			% Smooth with a 3-bin boxcar
			data = ccv(:,i,j);
			top = flipud(data(1:size(ccg,1),:));
			bottom = flipud(data(end-size(ccg,1)+1:end,:));
			data = [top;data;bottom];
			data = filter([1 1 1],3,data);
			n = size(data,1);
			d = n - size(ccg,1);
			start = d/2+1;
			stop = start + size(ccg,1) - 1;
			ccv(:,i,j) = data(start:stop);

			% Find the peak lag time and value
			[~,maxIndex] = max(ccv(:,i,j));
			tau(i,j) = median(t(maxIndex));
			c(i,j) = median(ccv(maxIndex,i,j));
			% Previous version of the code (discard?)
			% [~,maxIndex] = max(ccv(:,i,j));
			% tau(i,j) = t(maxIndex);
			% c(i,j) = median(ccv(max(1,maxIndex-3):min(end,maxIndex+3),i,j));

			% Keep only significantly correlated pairs
			if ~any(abs(ccv(:,i,j))>threshold)
				tau(i,j) = NaN;
			end
			
		end
	end

	nCorrelatedPairs = sum(~isnan(tau(:)));
	disp(['# significantly correlated pairs: ' int2str(nCorrelatedPairs)]);

	ccg = ccv;
	
end