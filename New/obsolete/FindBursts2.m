function bursts = FindBursts(spikes,varargin)

%FindBursts - Find bursts of action potentials in a spike train.
%
%  USAGE
%
%    bursts = FindBursts(spikes,<options>)
%
%    Find spike bursts by smoothing and thresholding the input single unit or
%    MUA spike train. Two thresholds are used: the upper threshold defines a
%    burst, while the lower threshold identifies the beginning and end of the
%    burst. See e.g. Wierzynski et al. (2009).
%
%    spikes         spike train (either single unit or MUA)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for burst beginning and end, in multiples of
%                   the stdev (default = [2 5])
%     'binSize'     bin size (in ms, default = 1)
%     'smooth'      smoothing kernel width (in ms, default = 16.667)
%    =========================================================================
%

% Copyright (C) 2016 by Ralitsa Todorova, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values

binSize = 0.001;
thresholds = [2 5];
smooth = (0.05/3)/binSize;

for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'thresholds',
			thresholds = varargin{i+1};
			if ~isdvector(thresholds,'#2','>0','<'),
				error('Incorrect value for property ''thresholds'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
			end
		case 'smooth',
		   smooth = varargin{i+1}/binSize;
		   if ~isdscalar(smooth,'>0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
		   end
		case 'binsize',
		   binSize = varargin{i+1};
		   if ~isdscalar(binSize,'>0'),
				error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
		   end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).']);
		end
end

lowThreshold = thresholds(1);
highThreshold = thresholds(2);

% Comute z-scored smoothed spike counts
spikes = spikes(:,1);
t = linspace(spikes(1),spikes(end),1+(spikes(end)-spikes(1))/binSize);
if spikes(end) > t(end), t = [t t(end)+binSize]; end
binned = histc(spikes,t);
binned = Smooth(binned,smooth);
z = zscore(binned);

% WAIT FOR RALY'S NEW VERSION, THEN START FROM HERE
bursts(:,[1 3]) = ToIntervals(t,z>lowThreshold);
aboveLowThreshold = ToIntervals(t,z>lowThreshold);
aboveHighThreshold = ToIntervals(t,z>highThreshold);
bursts = bursts(IntervalsIntersect(t(FindInterval(z>lowThreshold)), t(FindInterval(z>highThreshold))),:);

if max(diff(bursts(:,[1 3]),[],2))/binSize > 1000, % matrix method would require too much mamory
    for i=1:length(bursts),
        restricted = Restrict([t binned],bursts(i,[1 3]));
        [~,localMax] = max(restricted(:,2));
        bursts(i,2) = restricted(localMax(:,1));
    end
    return
end

[bursting, burstID] = InIntervals(t,bursts(:,[1 3]));
binned = binned(bursting); burstID = burstID(bursting); t = t(bursting);
ranks = CumSum(ones(size(binned)),[1;diff(burstID)]);
matrix = nan(size(bursts,1), max(ranks));
tmatrix = nan(size(bursts,1), max(ranks));
matrix(sub2ind(size(matrix),burstID,ranks)) = binned;
tmatrix(sub2ind(size(matrix),burstID,ranks)) = t;

[~,peakID] = max(matrix,[],2);
peaks = tmatrix(sub2ind(size(tmatrix),(1:length(peakID))',peakID));
bursts(:,2) = peaks;

end


%% HELPER FUNCTION
function interval = FindInterval(logical)

logical = [0; logical(:); 0];
starts = strfind(logical',[0 1])';
stops = strfind(logical',[1 0])'-1;

interval = [starts stops];
end
