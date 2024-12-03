function ripples = BatchFindRipples(session,channels, noisyChannel,excludeInterval)


% This is the algorithm for finding ripple events developed in May 2016 for
% Celine's data. The reasoning is to take the filtered (ripple-band) signal's
% envelope from multiple sources (carefully picked informative channels) together
% to find ripple events.
% To remove noise to prevent detecting noisy artefects as ripples, a noise
% envelope is removed from the signal (mean envelope of the signal channels).
% This noise envelope can be either:
%           - the ripple [100 250] Hz envelope of a noise channel, if provided
%           - if no noise channel is provided, the noise envelope is the envelope
%             of the original signal filtered in a non0ripple band [300 500] Hz
% Ralitsa Todorova 2016

%% Initial parameters

here = max(strfind(session,'/'));
try cd([session(1:here) '/events/']); catch cd(session(1:here));end
threshold = 1;
peakThreshold = 3;
if str2double(session(end-15:end-13)) == 272; peakThreshold = 4; end 
if str2double(session(end-15:end-13)) == 291; peakThreshold = 5; end %noisy
minInterRippleInterval = 0.030;
minRippleDuration = 0.020;
maxRippleDuration = 0.20;
aroundArtefact = 0.05; epsilon = 0.15; %for consolidating noisy periods
SetCurrentSession(session,'verbose','off');

if exist('noisyChannel','var') & ~isempty(noisyChannel)
    for i=1:length(noisyChannel)
    lfp = GetLFP(noisyChannel(i));
    filtered = FilterLFP(Detrend(lfp), 'passband', 'ripples');
    [~,amplitude] = Phase(filtered);
     noise(:,i)  = amplitude(:,2);
    end
    noisy = nanmean(noise,2);
else
    % If there was not a valid noiseChannel, use the envelopes for higher frequencies.
    % Artefacts tend to bleed to all frequencies, so this way they will be subtracted
%     noise = zeros(length(t),numel(channels));
    for i=1:length(channels)
        lfp = GetLFP(channels(i));
        filtered = FilterLFP(Detrend(lfp), 'passband', [300 500]);
        [~,amplitude] = Phase(filtered);
        noise(:,i) = amplitude(:,2);
    end
    noisy = nanmean(noise,2);
end


t = lfp(:,1);

% Get envelopes
envelope = zeros(length(t),numel(channels));

for i=1:numel(channels)
    lfp = GetLFP(channels(i));
    filtered = FilterLFP(Detrend(lfp), 'passband', 'ripples');
    [~,amplitude] = Phase(filtered);
    envelope(:,i) = amplitude(:,2);
end
%%
me = nanmean(envelope,2); %mean envelope

% Get rid of (saturated) noisy periods

badPeriods = [];
if any(abs(lfp(:,2))>=32766)% saturation @ 32766;
    disp('bad periods detected')
    bad = abs(lfp(:,2))>=32766;
    starts = strfind(bad', [0 1])';
    ends = strfind(bad', [1 0])';
    if starts(1)>ends(1); starts = [1; starts]; end
    if starts(end)>ends(end); ends = [ends; length(t)]; end
    badPeriods = [t(starts)-aroundArtefact t(ends)+aroundArtefact];
    % Consolidate bad periods that are separated by less than epsilon
    badPeriods = ConsolidateIntervals(badPeriods, 'epsilon', epsilon);
    disp(badPeriods(1:10,:))
elseif str2double(session(end-15:end-13)) == 291; % Rat 291;
    lfp = Detrend(GetLFP(channels(1))); 
    bad = abs(lfp(:,2))>=500;
    starts = strfind(bad', [0 1])';
    ends = strfind(bad', [1 0])';
    if starts(1)>ends(1); starts = [1; starts]; end
    if starts(end)>ends(end); ends = [ends; length(t)]; end
    badPeriods = [t(starts)-aroundArtefact t(ends)+aroundArtefact];
    % Consolidate bad periods that are separated by less than epsilon
    badPeriods = ConsolidateIntervals(badPeriods, 'epsilon', epsilon);
end

if exist('excludeInterval','var'), badPeriods = ConsolidateIntervals(sortrows([badPeriods;excludeInterval]), 'epsilon', epsilon); end
% if ~isempty(badPeriods),noisy(InIntervals(t,badPeriods)) = nan;end
% this = nanzscore(noisy); tis(isnan(this)) = 0; noiseRipples = helper_function(t,this,threshold,minInterRippleInterval,0,Inf,threshold);
% badPeriods = ConsolidateIntervals(sortrows([badPeriods;noiseRipples(:,[1 3])]), 'epsilon', epsilon);
if ~isempty(badPeriods),me(InIntervals(t,badPeriods)) = nan;end

% Find ripples
signal = me-noisy; % remove noise from envelope
signal(me>noisy) = nanzscore(signal(me>noisy));
% signal(nanzscore(noisy)>threshold) = 0;
signal(~(me>noisy)) = 0;

ripples = helper_function(t,signal,threshold,minInterRippleInterval,minRippleDuration,maxRippleDuration,peakThreshold);

%%Save ripples

% dlmwrite('ripples',r,'precision','%f8');
% SaveRippleEvents('ripples.rip.evt', r, 0,'overwrite', 'on');

function r= helper_function(t,signal,threshold,minInterRippleInterval,minRippleDuration,maxRippleDuration,peakThreshold)

% get events above threshold sd-s
sig = signal>threshold;
starts = strfind([0 sig' 0],[0 1])';
stops = strfind([0 sig' 0],[1 0])' - 1;
r = t([starts stops]);
if starts(1)==1, r(1,:) = []; end
if stops(end)==length(t), r(end,:) = []; end

% get the local minima around the event as the start/stop points
minima = t(FindLocalMinima(signal)); %original signal, no nans
minima = [minima(1:end-1) minima(2:end)];

if r(1)<minima(1),r(1,:) = [];end
if r(end)>minima(end),r(end,:) = [];end

starts = minima(CountInIntervals(r(:,1)+eps,minima)>0,1);
stops = minima(CountInIntervals(r(:,2)-eps,minima)>0,2);
r = [starts stops];

% merge close ones
iri = r(2:end,1) - r(1:end-1,2);
newSize =  r(2:end,2) - r(1:end-1,1);
toMerge = iri<minInterRippleInterval & newSize<maxRippleDuration;

while sum(toMerge)>0,
    % Get the index of the first ripple in a pair to be merged
    rippleStart = strfind([0 toMerge'],[0 1])';
    rippleEnd = rippleStart+1; % next ripple
    % Adjust end (incorporate second ripple into first)
    r(rippleStart,2) = r(rippleEnd,2);
    % Remove now-redundant second ripple
    r(rippleEnd,:) = [];
    iri = r(2:end,1) - r(1:end-1,2);
    newSize =  r(2:end,2) - r(1:end-1,1);
    toMerge = iri<minInterRippleInterval & newSize<maxRippleDuration;
end

% peak threshold check
[in,w] = InIntervals(t,r);
[peakV,peak] = Accumulate(w(in),signal(in),'mode','max');
ok = peakV>peakThreshold;
tIn = t(in);
peak = tIn(peak(ok));
r = [r(ok,1) peak r(ok,end) peakV(ok)];

% duration check
d = r(:,3) - r(:,1);
r(d<minRippleDuration | d>maxRippleDuration,:) = [];
