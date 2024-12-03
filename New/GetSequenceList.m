function varargout = GetSequenceList(spikes, intervals, varargin)

% [sequenceList, times, qualifiedIntervals, logical_qualifiedIntervals, nUniqueUnits] = GetSequenceList(spikes, intervals, options)
% Gives sequenceList matrix to use in PairOrder
% Threshold is the minimum number of different neurons that should discharge 
% within any interval to qualify as a sequence


% Raly & Céline

threshold = [];
repetitions = 'on';

for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help GetSequenceList">GetSequenceList</a>'' for details).']);
    end
    switch(lower(varargin{i}))
        case 'threshold'
            threshold = varargin{i+1};
%             if (threshold),
%                 error('Incorrect value for property ''threshold'' (type ''help <a href="matlab:help PairOrder">PairOrder</a>'' for details).');
%             end
		case {'repetition','repetitions'}
			repetitions = varargin{i+1};
			if ~isastring(repetitions,'on','off')
				error('Incorrect value for property ''repetition'' (type ''help <a href="matlab:help PairOrder">PairOrder</a>'' for details).');
			end
		case 'keep'
			keep = varargin{i+1};
			repetitions = 'off';
			if ~isastring(keep,'first','mid','last')
				error('Incorrect value for property ''keep'' (type ''help <a href="matlab:help PairOrder">PairOrder</a>'' for details).');
			end
		otherwise
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PairOrder">PairOrder</a>'' for details).']);
    end
end

[ok, intervalIDs] = InIntervals(spikes, intervals);
intervalIDs = intervalIDs(ok);
t = spikes(ok,1);
spikes = spikes(ok,2);
ranks = CumSum(ones(size(spikes)),[1;diff(intervalIDs)]);
sequenceList = nan(size(intervals,1), max(ranks));
sequenceList(sub2ind(size(sequenceList),intervalIDs,ranks)) = spikes;

times = nan(size(intervals,1), max(ranks));
times(sub2ind(size(sequenceList),intervalIDs,ranks)) = t;

if ~isempty(threshold)
    matrix = Accumulate([intervalIDs spikes],1,'size',[size(intervals,1) max(spikes)]);
    nUniqueUnits = sum(logical(matrix),2);
    sequenceList(nUniqueUnits<threshold,:) = [];
    qualifiedIntervals = intervals(nUniqueUnits>=threshold,:);
    logical_qualifiedIntervals = nUniqueUnits>=threshold;
else
    qualifiedIntervals = intervals;
    if nargout>2
        matrix = Accumulate([intervalIDs spikes]);
        nUniqueUnits = sum(logical(matrix),2);
        logical_qualifiedIntervals = nan;
    end
end

if strcmp(repetitions, 'off')
	if exist('keep','var')
    sequenceList = RemoveRepetitions(sequenceList,'mode','matrix','keep',keep);
	else
	sequenceList = RemoveRepetitions(sequenceList,'mode','matrix');
	end
end

for i=1:nargout
    if i==1
        varargout{1} = sequenceList;
    elseif i==2
        varargout{2} = times;
    elseif i==3
        varargout{3} = qualifiedIntervals;
    elseif i==4
        varargout{4} = logical_qualifiedIntervals;
    elseif i==5
        varargout{5} = nUniqueUnits;
    end
end

end
%% % HELPER FUNCTION 

function [clean,discarded,cleanIndices] = RemoveRepetitions(vector, varargin);

% Removes repetitions (e.g. bursting neuron). Default mode: all repetitions
% are removed (mode = 'all') - like 'unique' but without changing the order.
% An alternative mode is 'immediate':
% Also, there is mode 'matrix', in which a matrix may be provided.
% Repetitions within each row will be substituted by NaNs.
% Default: keep first value and remove later repetitions. Alternatively, 
% the last ('last') or the median ('mid') value can be kept with earlier 
% (as well as later for 'mid') repetitions discarded
% EXAMPLE:
% mode = 'immediate'
% [clean,cleanIndices,discarded] = RemoveRepetitions([1 1 3 4 5 5 5 2 1]),
% clean = 1 3 4 5 2 1;
% cleanIndices = 1 3 4 5 8 9;
% discarded = 0 1 0 0 0 1 1 0 0;

mode = 'all';
keep = 'first';

% Options
for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help RemoveRepetitions">RemoveRepetitions</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'mode',
            mode = lower(varargin{i+1});
            if ~isastring(mode,'immediate','all','matrix')
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help RemoveRepetitions">RemoveRepetitions</a>'' for details).');
			end
		case 'keep',
			keep = lower(varargin{i+1});
			if ~isastring(keep,'first','mid','last')
				error('Incorrect value for property ''keep'' (type ''help <a href="matlab:help RemoveRepetitions">RemoveRepetitions</a>'' for details).');
			end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help RemoveRepetitions">RemoveRepetitions</a>'' for details).']);
    end
end

if size(vector,2)>size(vector,1) && ~strcmp(mode, 'matrix')
    turn = true;
    vector = vector(:);
else
    turn = false;
end

if strcmp(mode, 'immediate')
    discarded = [false; diff(vector)==0];
    clean = vector(~discarded);
    cleanIndices = find(~discarded);
elseif strcmp(mode, 'matrix')
    % Each nan is unique, causing problems. Make them zeros to avoid such problems:
    discarded = false(size(vector));
    nans = isnan(vector);
    vector(isnan(vector)) = 0; % should be a safe value
    [count, value] = hist(vector', unique(vector));
    if size(count,1)==1, count = count(:); end % avoid problems in case a single sequence exists
    [elementIndex, sequenceIndex] = find(count>1);
    for i=unique(sequenceIndex)'
        repeatedValues = value(elementIndex(sequenceIndex==i));
        repeatedValues(repeatedValues==0) = []; %remove repeats due to NaNs
        for j=1:length(repeatedValues)
            thismany = sum(vector(i,:)==repeatedValues(j));
			if strcmp(keep,'first')
				discarded(i, find(vector(i,:)==repeatedValues(j),round(thismany-1),'last')) = 1; % letting the first one stay, the last n-1 will be discarded
			elseif strcmp(keep,'last')
				discarded(i, find(vector(i,:)==repeatedValues(j),round(thismany-1),'first')) = 1; % letting the last one stay, the first n-1 will be discarded
			else % case of leaving the middle one in
				indices = find(vector(i,:)==repeatedValues(j));
				indices(ceil(length(indices)/2)) = []; % keep the middle index
				discarded(i,indices) = 1;
			end
        end
    end
    cleanIndices = ~discarded;
    clean = vector;
    clean(discarded) = NaN;
    % Put the NaNs back in:
    clean(nans) = NaN;
else
    [count, value] = hist(vector, unique(vector));
    repeatedValues = value(count > 1);
    cleanIndices = find(~ismember(vector,repeatedValues));
	for i=1:length(repeatedValues)
		if strcmp(keep,'first') || strcmp(keep,'last')
			cleanIndices = [cleanIndices; find(vector==repeatedValues(i),1,keep)];
		else
			indices = find(vector==repeatedValues(i)); index = indices(ceil(length(index)/2));
			cleanIndices = [cleanIndices; index];
		end	
    end
    cleanIndices = sort(cleanIndices);
    discarded = ~Unfind(cleanIndices,length(vector));
    clean = vector(cleanIndices);
end


if turn
    clean = clean';
    discarded = discarded';
end
end