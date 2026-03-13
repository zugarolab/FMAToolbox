function [ccg,t,tau,c] = CCG(times,id,opt)

%CCG - Compute multiple cross- and auto-correlograms, or cross-covariances
%
%  USAGE
%
%    [ccg,t,tau,c] = CCG(times,id,<options>)
%
%    times          times of all events in seconds (see NOTEs); can also be a
%                   cell array where each cell contains event times for an ID
%    id             ID for each event (e.g. unit ID) from 1 to n (ignored
%                   if 'times' is a cell array)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration in s of each xcorrelogram (default = 2)
%     'nBins'       number of bins (default = duration/binSize)
%     'smooth'      smoothing size in bins (0 = no smoothing, default)
%     'groups'      group number (1 or 2) for each event, used to restrict
%                   cross-correlograms to pairs across two groups of events
%                   (see EXAMPLE #2 below)
%     'range'       range of possible IDs, in case some units have not fired
%                   a single spike within the interval, but they should not
%                   be omitted so that the IDs are preserved (2 values for
%                   normal (ungrouped) CCGs, or 4 values for grouped CCGs)
%     'mode'        'ccg' or 'ccv' (default = 'ccg')
%     'alpha'       significance level to determine correlated pairs (only
%                   for 'ccv' mode)
%     'totalTime'   recording duration in s (default = 0 is
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
    id (:,1) {mustBeInteger,mustBePositive}
    opt.binSize (1,1) {mustBePositive} = 0.01
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
if isscalar(id)
    id = ones(size(times));
    nIDs = 1;
else
    if length(times) ~= length(id)
        error('Parameters ''times'' and ''id'' have different lengths (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
    end
    nIDs = max(id);
end
if ~isempty(opt.groups) && length(times) ~= length(opt.groups)
    error('Property ''groups'' must have the same length as ''times'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if isempty(opt.groups) && ~ismember(numel(opt.range),[0,2])
    error('Incorrect value for property ''range'' (expected two elements). Type ''help <a href="matlab:help CCG">CCG</a>'' for details.')
end
if ~isempty(opt.groups) && ~ismember(numel(opt.range),[0,4])
    error('Incorrect value for property ''range'' (expected two elements). Type ''help <a href="matlab:help CCG">CCG</a>'' for details.')
end

% Default values
if isempty(opt.nBins)
    if isempty(opt.duration)
        opt.duration = 2;
    end
else
    if isempty(opt.duration)
        opt.duration = opt.nBins * opt.binSize;
    elseif opt.duration ~= opt.binSize * opt.nBins
        error('Incompatible ''duration'' and ''nBins'' parameters (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
    end
end
halfBins = round(opt.duration/opt.binSize/2);
opt.nBins = 2*halfBins + 1;
t = (-halfBins:halfBins).' * opt.binSize;
tau = [];
c = [];
if length(times) <= 1
    ccg = zeros(opt.nBins,nIDs,nIDs);
    return
end

% Sort events in time
[times,idx] = sort(times);
id = id(idx);
if ~isempty(opt.groups)
    opt.groups = opt.groups(idx);
end
if opt.totalTime == 0
    opt.totalTime = times(end) - times(1);
end

% CCG
counts = double(CCGEngine(times,id,opt.binSize,halfBins));
counts = reshape(counts,[opt.nBins nIDs nIDs]); % reshape to (time bins, IDs, IDs)

% Rate normalization: counts/numREFspikes/dt to convert to Hz
switch opt.norm
    case 'rate'
        for gg = 1:nIDs
            numREFspikes = sum(id==gg); % number of reference events for group
            counts(:,gg,:) = counts(:,gg,:)./numREFspikes./opt.binSize;
        end
end

% Add empty correlograms for IDs in requested range (CCGEngine output aleady considers processes in 1 : nIDs)
if numel(opt.range) == 2
    nIDsOld = nIDs;
    nIDs = max(opt.range);
    extraIDs = (nIDsOld + 1) : nIDs;
    counts(:,extraIDs,:) = zeros(opt.nBins, numel(extraIDs), nIDsOld);
    counts(:,:,extraIDs) = zeros(opt.nBins, nIDs, numel(extraIDs));
end

% Restrict the results to inter-group CCGs if requested
if ~isempty(opt.groups)
    if numel(opt.range) == 4
        opt.range = reshape(opt.range,[2 2])'; %making sure the values are arranged [min1 max1; min2 max2] just in case a single-row vector was provided.
        group1 = unique([id(opt.groups == 1) opt.range(1,:)]); % all the IDs of the neurons that have spiked, including the range provided by the user
        group2 = unique([id(opt.groups == 2) opt.range(2,:)]);
    else
        group1 = unique([id(opt.groups == 1)]);
        group2 = unique([id(opt.groups == 2)]);
    end
    group1 = min(group1) : max(group1); 
    group2 = min(group2) : max(group2); 
    nGroup1 = numel(group1);
    nGroup2 = numel(group2);
    ccg = zeros(opt.nBins,nGroup1,nGroup2);
    for i = 1:nGroup1
        for j = 1:nGroup2
            ccg(:,i,j) = Smooth(flipud(counts(:,group1(i),group2(j))),opt.smooth);
        end
    end
else
    ccg = zeros(opt.nBins,nIDs,nIDs);
    % Compute corr(A,B) for each unique unordered pair (A,B)
	for g1 = 1:nIDs
		for g2 = g1:nIDs
			ccg(:,g1,g2) = Smooth(flipud(counts(:,g1,g2)),opt.smooth);
		end
	end
	% corr(B,A) and corr(B,A) symmetric
	for g1 = 1:nIDs
		for g2 = 1:g1-1
			ccg(:,g1,g2) = flipud(squeeze(ccg(:,g2,g1)));
		end
	end
end


if strcmp(opt.mode,'ccv')

	% Determine mean event rate for each ID
	eventRate = zeros(nIDs,1);
	for i = 1:nIDs
		eventRate(i) = sum(id==i)/opt.totalTime;
	end

	% Determine standardized cross-covariances
	ccv = zeros(size(ccg));
	tau = zeros(size(ccg,2),size(ccg,3));
	c = zeros(size(ccg,2),size(ccg,3));

	nPairs = size(ccg,2)*size(ccg,3);
	disp(['# pairs: ' int2str(nPairs)]);

	threshold = sqrt(2)*erfinv(1-(opt.alpha/length(t)));

	for i = 1:size(ccg,2)
		for j = 1:size(ccg,3)
		
			% Compute and normalize CCVs from CCGs
			if ~isempty(opt.groups)
				rate = eventRate(group1(i))*eventRate(group2(j));
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