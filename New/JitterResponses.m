function [thresholdHigh, thresholdLow, observed, t] = JitterResponses(spikes, events, durations, varargin)

% [thresholdHigh, thresholdLow, observed, t] = JitterResponses(spikes,  events, durations, varargin)
% This function computes significance thresholds for the spike responses
% within provided intervals using a bootstrap procedure.
% To this end, the function adds a uniform jitter to the prodivided spike train,
% and then computes a historam of the jittered spikes' response within the intervals
% provided.
% t=time from interval start
%
%
%
%    ================================================================================
%     Properties    Values
%    --------------------------------------------------------------------------------
%     'jitter'          absolute value of the jittering interval (default=0.1s)
%     'p'               confidence interval to use (default=0.05)
%     'binSize'         the size of the histogram bins (default=0.01s)
%     'nIterations'     number of jitters performed to estimate thresholds (default=1000)
%

% Default values
p = 0.05;
nIterations = 1000;
binSize = 0.01;
jitter = 0.1;
% durations = [-0.01 0.02];

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help JitterResponses">JitterResponses</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'jitter',
            jitter = varargin{i+1};
            if ~isvector(jitter) || length(jitter) ~= 1,
                error('Incorrect value for property ''jitter'' (type ''help <a href="matlab:help JitterResponses">JitterResponses</a>'' for details).');
            end
        case 'p',
            p = varargin{i+1};
            if ~isvector(p) || length(p) ~= 1,
                error('Incorrect value for property ''threshold'' (type ''help <a href="matlab:help JitterResponses">JitterResponses</a>'' for details).');
            end
        case 'binsize',
            binSize = varargin{i+1};
            if ~isvector(binSize) || length(binSize) ~= 1,
                error('Incorrect value for property ''nIterations'' (type ''help <a href="matlab:help JitterResponses">JitterResponses</a>'' for details).');
            end
        case 'niterations',
            nIterations = varargin{i+1};
            if ~isvector(nIterations) || length(nIterations) ~= 1,
                error('Incorrect value for property ''nIterations'' (type ''help <a href="matlab:help JitterResponses">JitterResponses</a>'' for details).');
            end
        otherwise,new
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help JitterResponses">JitterResponses</a>'' for details).']);
    end
end

%%

intervals = [events+durations(1) events+durations(2)];
nIntervals = size(intervals,1);
spikes = sortrows(spikes);

if size(spikes,2)==2,
    id = spikes(:,2);
    t = spikes(:,1);
else
    t = spikes;
    id = ones(size(spikes));
end

% check if intervals overlap
consolidated = ConsolidateIntervals([intervals(:,1)-jitter intervals(:,2)+jitter],'strict','on');
if size(consolidated,1)==nIntervals, % no overlap
    [in,w] = InIntervals(t,[intervals(:,1)-jitter intervals(:,2)+jitter]);
    t = t(in) - events(w(in));
    id = id(in);
else %if intervals overlap
    t0 = t; id0 = id;
    t = [];
    for i=1:nIntervals,
        in = InIntervals(t0,[intervals(i,1)-jitter intervals(i,2)+jitter]);
        t = [t;t0(in) - events(i,1)];
        id = [id; id0(in)];
    end
end

interval = durations;

[bins intervalID] = SplitIntervals(interval,'pieceSize',binSize);
first = [1;diff(intervalID)];
binID = CumSum(ones(size(first)),first);
nUnits = max(id);
nBins = max(binID);


h = zeros(nUnits,nBins,nIterations); % histograms for all the jitters
n = size(t,1);
for i=1:nIterations,
    jittered = t + jitter * (2*rand(n,1)-1);
    [in,w] = InIntervals(jittered,bins);
    if nUnits==1,
        h(1,:,i) = Accumulate(binID(w(in)),1,nBins);
    else
        h(:,:,i) = Accumulate([id(in) binID(w(in))],1,[nUnits nBins]);
    end
end

[in,w] = InIntervals(t,bins);
if nUnits==1,
    observed(1,:) = Accumulate(binID(w(in)),1,nBins);
else
    observed = Accumulate([id(in) binID(w(in))],1,[nUnits nBins]);
end

thresholdHigh = quantile(max(h,[],3),1-p,2);
thresholdLow = quantile(min(h,[],3),p,2);



t = mean(bins,2);














