function [estimations,actual,errors,average,lambda] = ReconstructPosition(positions,spikes,windows,varargin)

% Bayesian reconstruction of positions from spike trains.
%
% Instantaneous positions are reconstructed using a Bayesian algorithm.
% Instantaneous population firing rates can be estimated either over fixed time
% windows, or over fractions of the theta cycle (or of any other brain rhythm).
% Similarly, positions will be reconstructed either over time or phase windows.
% The model is first trained on a subset of the data, then tested on the rest.
%
% USAGE
%
%    [estimations,positions,errors,average] = ReconstructPosition(positions,spikes,windows,<options>)
%
%    positions      linear or two-dimensional positions <a href="matlab:help samples">samples</a>, in [0..1]
%    spikes         list of (t,ID) couples (obtained via e.g. <a href="matlab:help GetSpikes">GetSpikes</a>,
%                   using numbered output). 
%                   Alternatively, provide a structure with spikes.timestamps in cell format
%    windows        exact window intervals
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'training'    time interval over which the model should be trained
%                   (see NOTE below for defaults)
%     'type'        two letters (one for X and one for Y) indicating which
%                   coordinates are linear ('l') and which are circular ('c')
%                   - for 1D data, only one letter is used (default 'll')
%     'nBins'       firing curve or map resolution (default = 200 for 1D and
%                   [200 200] for 2D data)
%     'smooth'      smooth that will be applied to the firing maps (default
%                   = 5 for 1D and 10 for 2D positions)
%     'prior'       taking the map occupancy into account or not (default = 'off')
%                   This is what makes the reconstruction Bayesien, but if animal
%                   occupancy is very skewed, the recustruction would be good even
%                   without informative spikes, which is why the default is off.
%     'interpolate' When estimating errors, compute the actual position by
%                   interpolating the positions data to the decoding window
%                   (default = 'on'). If positions data is not continuous (i.e.
%                   there are discontinuities), set to 'off' to take the position
%                   value closest in time without interpolating.
%     'lambda'      Give a previously computed lambda and avoid training
%                   the model on the place fields. 'training' will be
%                   ignored if this input is provided.
%     'id'          bin number for each window, used to compute the avarage
%                   error (for example, phases for each window, to compute the
%                   average error by phase). This should be a vector with one
%                   element (integer) per window (default - no error computed).
%    =========================================================================
%
%   OUTPUT
%
%     estimations   estimated position across time or phase windows
%     positions     real position across time or phase windows
%     errors        estimation error across time or phase windows
%     average       average estimation error in each phase window
%
% Copyright (C) 2012-2015 by Michaël Zugaro, (C) 2012 by Karim El Kanbi (initial, non-vectorized implementation),
% (C) 2015-2025 by Ralitsa Todorova (log-exp fix, flattening dimensions, options)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
nBins = 200;
training = 1;
type = 'l';
nDimensions = 1;
intervalID = [];
prior = 'off';
interpolate = 'on';
minSpikes = 100;
smooth = 5;

% If spikes are provided in cell format, transform them to matrix
if isstruct(spikes) && isfield(spikes,'times')
    spikes = spikes.times';
    spikes = Group(spikes{:});
end

% Check number of parameters
if nargin < 3
    builtin('error','Incorrect number of parameters (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end

% Check parameter sizes
if ~isempty(positions) && ~isdmatrix(positions)
    builtin('error','Incorrect positions (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end
if ~isdmatrix(spikes,'@2') && ~isdmatrix(spikes,'@3')
    builtin('error','Incorrect spikes (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end
if size(positions,2) >= 3
    nDimensions = 2; type='lll'; nBins = repmat(nBins,1,2);
end
if ~isdmatrix(windows,'@2')
    builtin('error','Incorrect value for property ''windows'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end

if nDimensions==2, smooth = 10; end

% Parse parameters
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        builtin('error',['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
    end
    switch(lower(varargin{i}))
        case 'training'
            training = varargin{i+1};
            if ~isdvector(training,'<') && ~isdmatrix([1 2],'@2','>0')
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
            end
        case 'nbins'
            nBins = varargin{i+1};
            if ~isiscalar(nBins) && ~isdmatrix([100 100],'@2','>0')
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
            end
        case 'smooth'
            smooth = varargin{i+1};
            if ~isdvector(smooth,'>=0') || length(smooth) > 2
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
            end
        case 'prior'
            prior = varargin{i+1};
            if ~isastring(prior,'on','off')
                error(['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
            end
        case 'lambda'
            lambda = varargin{i+1};
            if ~ismatrix(lambda) || any(lambda(:)<0)
                error(['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
            end
        case 'interpolate'
            interpolate = varargin{i+1};
            if ~isastring(interpolate,'on','off')
                error(['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
            end
        case 'id'
            intervalID = varargin{i+1};
            if ~isivector(intervalID,'>0')
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
            end
        case 'minspikes'
            minSpikes = varargin{i+1};
            if ~isiscalar(minSpikes)
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
            end
        case 'type'
            type = lower(varargin{i+1});
            if (nDimensions == 1 && isastring(type(1),'c','l'))
                type = type(1);
            elseif nDimensions == 2 && length(type)>1 && isastring(type(1:2),'cl','cc','ll','lc')
                type = type(1:2);
            else
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
            end
        case 'window'
            builtin('error',['Property ''' varargin{i} ''' no longer supported (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
        case 'mode'
            builtin('error',['Property ''' varargin{i} ''' no longer supported (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
        case 'px'
            builtin('error',['Property ''' varargin{i} ''' no longer supported (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
        otherwise
            builtin('error',['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
    end
end

% Make sure that position timestamps are unique
positions = sortrows(positions);
ignore = diff(positions(:,1))<=0;
positions(ignore,:) = [];



% Convert from legacy format for backward compatibility with previous versions of the code (spikes)
if isdmatrix(spikes,'@3')
    % List units, assign them an ID (number them from 1 to N), and associate these IDs with each spike
    % (IDs will be easier to manipulate than (group,cluster) pairs in subsequent computations)
    [units,~,i] = unique(spikes(:,2:end),'rows');
    nUnits = length(units);
    index = 1:nUnits;
    id = index(i)';
    spikes = [spikes(:,1) id];
    warning('Spikes were provided as Nx3 samples - this is now obsolete (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
else
    id = spikes(:,2);
    nUnits = max(id);
end
%% TRAINING

if exist('lambda','var')
    if size(lambda,1) ~= nBins || size(lambda,2)~=nUnits,
        error('Wrong size for input ''lambda''. Aborting...');
    end
    tooFew = nansum(lambda)'==0;
else
    % Defaults (training)
    if isdscalar(training) %If 'training' was provided as a portion, convert it to an interval
        training = [0 positions(1,1)+training*(positions(end,1)-positions(1,1))];
    end
    
    trainingPositions = Restrict(positions,training,'shift','on');
    trainingSpikes = Restrict(spikes,training,'shift','on');
    
    % Apply minSpikes threshold to ignore units with very few spikes (which
    % would lead to noisy estimates of lambda)
    nSpikes = Accumulate(trainingSpikes(:,2),1,'size',max(id));
    tooFew = nSpikes<minSpikes;
    
    % Compute average firing probability lambda for each unit (i.e. firing maps)
    lambda = nan(prod(nBins),nUnits);
    for i = 1:nUnits
        unit = trainingSpikes(:,2) == i;
        s = trainingSpikes(unit,1);
        if size(trainingPositions,2)<3
            map = Map(trainingPositions,s,'nbins',nBins,'smooth',smooth,'type',[type 'l']);
            % extra letter for 'type' required as input for 'Map' even though this extra 'l' does not refer to anything in the case of a point process (spikes) provided
            map.z = map.z';
        else
            map = Map(trainingPositions,s,'nbins',nBins,'smooth',smooth,'type',type);
        end
        try
        lambda(:,i) = map.z(:); % squeeze space dimensions to 1
        catch
                    keyboard
        end
        lambda(:,i) = lambda(:,i) +eps;
    end
end

if strcmp(prior,'off')
    Px = ones(prod(nBins),1)/prod(nBins); % set a uniform prior
else
    if exist('map','var')
        % Compute occupancy probability P(x) (i.e. normalized occupancy map)
        Px = map.time(:);
    else
        % Defaults (training)
        if isdscalar(training) %If 'training' was provided as a portion, convert it to an interval
            training = [0 positions(1,1)+training*(positions(end,1)-positions(1,1))];
        end
        trainingPositions = Restrict(positions,training,'shift','on');
        map = Map(trainingPositions,spikes(:,1),'nbins',nBins,'smooth',smooth,'type',[type 'l']);
        Px = map.time(:);
        Px = Px ./ sum(Px);
    end
end

% We are squeezing the spatial dimensions into one dimension of binsize nBins:
nBinsPerDim = nBins;
nBins = prod(nBins);

%% TEST

nWindows = size(windows,1);

spikecount = zeros(nUnits,nWindows);
for i = 1:nUnits
    % Population spike count matrix
    spikecount(i,:) = CountInIntervals(spikes(id==i),windows);
end

% In rare cases there may be a neuron that didn't fire at all during training.
% This results in multiplying/dividing by zero, so they should be ignored
silent = sum(lambda)==0 | tooFew';

spikecount(silent,:) = [];
lambda(:,silent) = [];

% Sometimes a unit has such a clear place field that the probability of it
% firing farther away is estimated to be zero. This leads to errors because
% the loglambda would be -Inf, which would be estimated to be NaN even when
% the neuron did not fire (-Inf*0 = NaN). To avoid this error, a minimal
% probability is added to all place maps. 
lambda = lambda+eps;

% preset estimations to be uniform
estimations = ones(nBins,nWindows)/nBins;
% They should remain uniform in windows in which there were no spikes; keep track of such windows:
uniform = sum(spikecount)'==0;
spikecount = spikecount(:,~uniform);
dt = diff(windows(~uniform,:),[],2);

% 'estimation' is virtually P(position|spikecount) in the Bayesian formula
% P(position|spikecount) =  P(spikecount|position)*P(position) / P(spikecount) or
% P(n|x)*P(x) / P(n)
% where
% P(x) is the occumpancy map,
% P(n) = sum over x of P(n|x)*P(x)
% P(n|x) = (dt*lambda).^n./factorial(n).*exp(-dt*lambda) for a Poisson approximation
% For P(n|x), we compute the three terms numirator/denominator)*multiplier separately, where
% numerator = (dt*lambda).^n
% denominator = factorial(n)
% multiplier = exp(-dt*lambda)

% To avoid overflow errors due to large values of (dt*lambda).^n and factorial(n), we compute the exponentials of their logs
% numerator = exp(n*log(dt)+n*log(lambda))
% denominator = exp(logfactorial(n))
nlogdt = sum(bsxfun(@times,spikecount,log(dt)'));
nloglambda = log(lambda)*spikecount;
% numerator = exp(bsxfun(@plus,nlogdt,nloglambda));
% denominator = exp(sum(logfactorial(spikecount)));

% Note that as we are computing these probabilities for the whole population
% at once, we are adding neuron-specific terms within exponentials (equivalent
% to multiplying them outside of the exponential).

% multiplier = exp(-sum(lambda,2)*dt');

% P(n|x) = exp(n.*log(dt*lambda)-logfactorial(n)-dt*lambda);
% P(n|x) = [exp(n.*log(dt*lambda))]*[exp(-logfactorial(n))]*[exp(-dt*lambda)];
% P(n|x) = x1*x2*x3

% Pnx = bsxfun(@rdivide,numerator.*multiplier,denominator);

% Again, to avoid overflow it's easier to move the exponents to outside of the numerator.*multiplier/denominator calculation:
numerator = (bsxfun(@plus,nlogdt,nloglambda));
denominator = (sum(logfactorial(spikecount)));
multiplier = -sum(lambda,2)*dt';

Pnx = exp(bsxfun(@minus,numerator+multiplier,denominator));
PnxPx = bsxfun(@times,Pnx,Px);
Pn = nansum(PnxPx);
estimations(:,~uniform) = bsxfun(@rdivide,PnxPx,Pn);

% reshape to original (non-flattened) spatial dimensions
if nDimensions>1
    estimations = reshape(estimations,[nBinsPerDim size(estimations,2)]);
end

% make sure no nans remain
if any(isnan(estimations(:)))
    nans = double(isnan(estimations)); nans(nans==0) = nan;
    % substitute NaNs with uniform probability
    remainingProbability = 1-nansum(estimations);
    remainingProbability(abs(remainingProbability)<0.0000000001) = 0;
    nans = remainingProbability./nansum(nans).*nans;
    estimations(isnan(estimations)) = nans(isnan(estimations));
end

if nargout==1
    return
end

% Estimation error
windowCenters = mean(windows,2);
if strcmp(interpolate,'on')
    interpolated = nan(nWindows,nDimensions);
    for dim = 1:nDimensions
        if strcmp(type(1),'l')
	  x=positions(:,1);
	  
	  y=positions(:,1+dim);
	  %         [x, index] = unique(x);
	  %          y=y(index);
	  interpolated(:,dim) = interp1(x,y,windowCenters);
        else
	  unwrapped = unwrap(2*pi*positions(:,1+dim)/max(positions(:,1+dim)));
	  x=positions(:,1);
	  y=unwrapped;
	  %         [x, index] = unique(x);
	  %         y=y(index);
	  interpolated(:,dim) = wrap(interp1(x,y,windowCenters),2)/(2*pi);
        end
    end
    positions = interpolated;
else
    positions = positions(FindClosest(positions(:,1),windowCenters),2:end);
end
errors = nan(size(estimations));
actual = positions;

% in a 2D matrix, 'x' (first positions columnn) represents the 2nd dimension, while 'y' (second positions column) represents the 1st dimension. Flip them!
corrected = positions;
if nDimensions==2
    corrected = positions(:,[2 1]);
end

ind = cell(1,nDimensions);
[ind{:}, timebin] = ind2sub(size(estimations),(1:numel(estimations))');
for dim = 1:nDimensions
    shift = Bin(corrected(:,dim),[0 1],nBinsPerDim(dim));
    index = (ind{dim} - shift(timebin)) +round(nBinsPerDim(dim)/2);
    if strcmp(type,'c') % if data is circular, wrap estimation around
        ind{dim} = mod(index-1,nBinsPerDim(dim))+1;
    else
        index(index<1 | index>nBinsPerDim(dim)) = nan;
        ind{dim} = index;
    end
end
index = sub2ind(size(errors),ind{:},timebin);
ok = ~isnan(index);
errors(index(ok)) = estimations(ok);
if any(isnan(errors(:)))
    flat = reshape(errors,[],nWindows);
    nans = double(isnan(flat)); nans(nans==0) = nan;
    % substitute NaNs with uniform probability
    remainingProbability = 1-nansum(flat);
    nans = remainingProbability./nansum(nans).*nans;
    flat(isnan(flat)) = nans(isnan(flat));
    errors = reshape(flat,size(errors));
end

if nargout<4
    return
end

if isempty(intervalID)
    intervalID = ones(nWindows,1);
%     warning('No id-s provided. Average error computed for all bins');
end

average = cell(max(intervalID),1);
for i=1:max(intervalID)
    average{i,1} = eval(['nanmean(errors(' repmat(':,',1,nDimensions) 'intervalID==i),nDimensions+1)']);
end
average = cat(nDimensions+1,average{:});

if any(silent) && nargout>4
    l = lambda;
    lambda = nan(size(lambda,1),length(silent));
    lambda(:,~silent) = l;
end

end
% ------------------------------- Helper functions -------------------------------

function data = logfactorial(data)

% We compute log(n!) as the sum of logs, i.e. log(n!) = sum log(i) for i=1:n
% First determine the largest n in the array
m = max(data(:));
% Create a look-up vector of sum log(i) for each i up to the largest n
sums = [0 cumsum(log(1:m))];
% Look-up the value for each item in the array
data(:) = sums(data+1);
end

function [indices,values] = FindClosest(reference, query, mode)

%FindClosest - field the indices of a reference vector closest to a query vector
% Look up a query vector in a reference vector, and for each value in
% the query vector, return the index of the closest value in the
% reference vector.
%
% USAGE
%
%    [indices,values] = FindClosest(reference, query, mode)
%
%    reference      reference signal which will be used as the table look-up
%    query          query signal. For each value of this signal, the function
%                   will look for the closest value in the reference signal
%    mode           a string ('either','higher', or 'lower') indicating
%                   whether the function should find the closest value in
%                   reference higher than the query ('higher) or lower than
%                   the query ('lower') or the closest value regardless of
%                   the direction of the difference (default = 'either')
%
%
% EXAMPLE
% indices = FindClosests(spikes(:,1),deltas(:,2)) will return the indices of spikes closest to the delta wave peak
% d = deltas(:,2) - spikes(indices) will contain all the minimal distances between delta wave peaks and spikes (one distance per each delta wave)
%
% NOTE
% This is an optimized implementation of computing the following trivial code:
% for i=1:length(reference), indices(i,1) = find(abs(query-reference(i)) == min(abs(query-reference(i))), 1, 'first'); end
% values = reference(indices);
%
% Copyright (C) 2016-2019 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ~exist('mode','var'),mode = 'either'; end
if isempty(reference) || isempty(query), indices = []; values = []; return; end
if length(reference)==1,indices = 1; values=reference; return; end

nans = isnan(reference);
nonnans = find(~nans);
reference(nans) = [];

% Make sure values to be matched are unique
[u,i] = unique(reference);
% Find closest index
if strcmp(mode,'higher')
    indices = ceil(interp1(u,(1:length(u))',query));
elseif strcmp(mode,'lower')
    indices = floor(interp1(u,(1:length(u))',query));
else
    indices = round(interp1(u,(1:length(u))',query));
end

% Closest value to an undershooting is the smallest one in the reference list
indices(indices<1 | query<min(u))=1;
% Closest value to an overshooting query is the largest one in the reference list
indices(indices>length(i) | query>max(u)) = length(i);

indices = i(indices);
values = reference(indices);
indices = nonnans(indices);
end
