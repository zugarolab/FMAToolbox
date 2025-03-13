function [ccg,t,tau,c] = CCG(times,id,varargin)

%CCG - Compute multiple cross- and auto-correlograms, or cross-covariances
%
%  USAGE
%
%    [ccg,t,tau,c] = CCG(times,id,<options>)
%
%    times          times of all events in seconds (see NOTE below)
%                   (alternate) - can be {nCells} array of [nSpikes] 
%                   containing the spiketimes for each cell 
%    id             ID for each event (e.g. unit ID) from 1 to n. Ignored
%                   if "times" is a cell (alternate above)
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
%     'alpha'       significance level to determine correlated pairs
%     'totalTime'   total recording duration in s (if different from the
%                   default = max(times) - min(times))
%     'Fs'          sampling rate (default = 20000 Hz)
%     'norm'        normalization of the CCG, 'counts' or 'rate'
%                   'counts' gives raw event/spike count,
%                   'rate' returns CCG in units of spks/second (default: counts)
%    =========================================================================
%
%  NOTES
%
%    The size of the cross-correlograms can be supplied either as a duration
%    (property 'duration') or as an number of bins (property 'nBins').
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
%
%  NOTE
%
%    Parameters 'times', 'id' and 'group' can be obtained using <a href="matlab:help CCGParameters">CCGParameters</a>.
%    As a special case, when computing the correlograms of spike trains, one
%    can use the output of <a href="matlab:help GetSpikes">GetSpikes</a> either directly or in combination with
%    <a href="matlab:help CCGParameters">CCGParameters</a>. See EXAMPLES below.
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
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
d = 2;
duration = [];
binSize = 0.01;
nBins = [];
smooth = 0;
groups = [];
mode = 'ccg';
alpha = 0.05;
IDrange=[];
normtype = 'counts';
Fs = 1/20000;

% Option for spike times to be in {Ncells} array of spiketimes DL2017
if iscell(times)
    clear id
    numcells = length(times);
    for cc = 1:numcells
        id{cc}=cc.*ones(size(times{cc}));
    end
    times = cat(1,times{:}); id = cat(1,id{:});
end

%Sort
[times,idx] = sort(times);
id = id(idx);

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdvector(times),
	error('Parameter ''times'' is not a real-valued vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(id) && ~isdvector(id),
	error('Parameter ''id'' is not a real-valued scalar or vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(id) && length(times) ~= length(id),
	error('Parameters ''times'' and ''id'' have different lengths (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
id = id(:);
times = times(:);
totalTime = max(times)-min(times);

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'binsize',
			binSize = varargin{i+1};
			if ~isdscalar(binSize,'>0'),
				error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'nbins',
			nBins = varargin{i+1};
			if ~isiscalar(nBins,'>0'),
				error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'duration',
			duration = varargin{i+1};
			if ~isdscalar(duration,'>0'),
				error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'smooth',
			smooth = varargin{i+1};
			if ~isdscalar(smooth,'>=0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'groups',
			groups = varargin{i+1};
			if ~isempty(groups) && ~isdvector(groups) && length(times) ~= length(groups)
				error('Incorrect value for property ''groups'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'alpha',
			alpha = varargin{i+1};
			if ~isdscalar(alpha,'>0'),
				error('Incorrect value for property ''alpha'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'mode',
			mode = varargin{i+1};
			if ~isastring(mode,'ccg','ccv'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
            end
        case 'totaltime',
            totalTime = varargin{i+1};
            if ~isdscalar(totalTime,'>0'),
                error('Incorrect value for property ''totaltime'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
            end
        case 'range',
            IDrange = varargin{i+1};
            if ~isdvector(IDrange),
                error('Incorrect value for property ''range'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
            end
            if ~isempty(groups) && numel(IDrange) ~= 4,
                error('Incorrect value for property ''range'' (expected four elements). Type ''help <a href="matlab:help CCG">CCG</a>'' for details.')
            end
            if isempty(groups) && numel(IDrange) ~= 2,
                error('Incorrect value for property ''range'' (expected two elements). Type ''help <a href="matlab:help CCG">CCG</a>'' for details.')
            end
        case {'norm','normtype'}
            normtype = varargin{i+1};
            if ~isastring(normtype,'counts','rate'),
                error('Incorrect value for property ''normtype'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
            end
        case 'fs',
            % ignore, but prevent error to allow compatibility with CellExplorer
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
	end
end

% Determine binSize/duration
if isempty(nBins),
	if isempty(duration),
		duration = d;
	end
else
	if isempty(duration),
		duration = nBins*binSize;
	elseif duration ~= binSize*nBins,
		error('Incompatible ''duration'' and ''nBins'' parameters (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
	end
end

tau = [];
c = [];

% Number of IDs, number of bins, etc.
if length(id) == 1,
	id = ones(length(times),1);
	nIDs = 1;
else
 	nIDs = max(id);
end
halfBins = round(duration/binSize/2);
nBins = 2*halfBins+1;
t = (-halfBins:halfBins)'*binSize;

if length(times) <= 1,
	return
end

% Sort events in time and compute CCGs
[times,i] = sort(times);
id = id(i);
if ~isempty(groups),
	groups = groups(i);
end

if length(times) <= 1,
    ccg = zeros(nBins,nIDs,nIDs);
    return
end

try
    counts = double(CCGEngine(times,id,binSize,halfBins));
catch
    problem = true; count = 1;
    while problem
        count = count+1;
        if count>length(times), error('Need to fix CCGEngine'); end
        try
            x = ([count:length(times)]); counts = double(CCGEngine(times(x),id(x),binSize,halfBins));
            problem = false;
        end
    end
end

% Reshape the results
counts = reshape(counts,[nBins nIDs nIDs]);

%Rate normalization: counts/numREFspikes/dt to put in units of spikes/s. DL
switch normtype
    case 'rate'
        for gg = 1:nIDs
            numREFspikes = sum(id==gg);%number of reference events for group
            counts(:,gg,:) = counts(:,gg,:)./numREFspikes./binSize;
        end
end

% Add empty correlograms for IDs within customised range which did not spike
% As CCGEngine automatically creates empty correlograms for IDs within the
% range of the id-s in the data, only correlograms for id-s above the
% maximal id in the data need to be added.
if numel(IDrange) == 2,
    nIDs = max(IDrange);
    extraIDs = (max(id) + 1) : max(IDrange);
    counts(:,extraIDs,:) = zeros(nBins, numel(extraIDs), max(id));
    counts(:,:,extraIDs) = zeros(nBins, max(IDrange), numel(extraIDs));
end

% Restrict the results to inter-group CCGs if requested
if ~isempty(groups),
    if numel(IDrange) == 4,
        IDrange = reshape(IDrange,[2 2])'; %making sure the values are arranged [min1 max1; min2 max2] just in case a single-row vector was provided.
        group1 = unique([id(groups == 1) IDrange(1,:)]); % all the IDs of the neurons that have spiked, including the range provided by the user
        group2 = unique([id(groups == 2) IDrange(2,:)]);
    else
        group1 = unique([id(groups == 1)]);
        group2 = unique([id(groups == 2)]);
    end
    group1 = min(group1) : max(group1); 
    %The ccg will return output for all the IDs within the range regardless
    %of whether there were spikes or not
    %This is in order not to omit any silent IDs and thus keep the original
    %ID number in the ccg output.
    group2 = min(group2) : max(group2); 
	nGroup1 = numel(group1);
	nGroup2 = numel(group2);
	ccg = zeros(nBins,nGroup1,nGroup2);
	for i = 1:nGroup1,
		for j = 1:nGroup2,
			ccg(:,i,j) = Smooth(flipud(counts(:,group1(i),group2(j))),smooth);
		end
	end
else
	ccg = zeros(nBins,nIDs,nIDs);
	% Compute corr(A,B) for each unique unordered pair (A,B)
	for g1 = 1:nIDs,
		for g2 = g1:nIDs,
			ccg(:,g1,g2) = Smooth(flipud(counts(:,g1,g2)),smooth);
		end
	end
	% corr(B,A) and corr(B,A) symmetric
	for g1 = 1:nIDs,
		for g2 = 1:g1-1,
			ccg(:,g1,g2) = flipud(squeeze(ccg(:,g2,g1)));
		end
	end
end


if strcmp(mode,'ccv'),

	% Determine mean event rate for each ID
	eventRate = zeros(nIDs,1);
	for i = 1:nIDs,
		eventRate(i) = sum(id==i)/totalTime;
	end

	% Determine standardized cross-covariances
	ccv = zeros(size(ccg));
	tau = zeros(size(ccg,2),size(ccg,3));
	c = zeros(size(ccg,2),size(ccg,3));

	nPairs = size(ccg,2)*size(ccg,3);
	disp(['# pairs: ' int2str(nPairs)]);

	threshold = sqrt(2)*erfinv(1-(alpha/length(t)));

	for i = 1:size(ccg,2),
		for j = 1:size(ccg,3),
		
			% Compute and normalize CCVs from CCGs
			if ~isempty(groups),
				rate = eventRate(group1(i))*eventRate(group2(j));
			else
				rate = eventRate(i)*eventRate(j);
			end
			ccv(:,i,j) = sqrt((binSize*totalTime)/rate) * (ccg(:,i,j)/(binSize*totalTime)-rate);

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
			if ~any(abs(ccv(:,i,j))>threshold),
				tau(i,j) = NaN;
			end
			
		end
	end

	nCorrelatedPairs = sum(~isnan(tau(:)));
	disp(['# significantly correlated pairs: ' int2str(nCorrelatedPairs)]);

	ccg = ccv;
	
end
