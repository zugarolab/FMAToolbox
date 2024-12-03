function [pieces ids] = SplitIntervals(intervals,varargin)

%SplitIntervals - Splits intervals into even pieces. 
% 
%There are two different ways to call this function: by either providing 'pieceSize' or 'nPieces'.
%
% USAGE
%
%   [pieces ids] = SplitIntervals(intervals,<options>)
%
%     intervals     list of [start stop] intervals
%     <options>     optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'pieceSize'   the required size of the split intervals
%     'nPieces'     alternatively, one can provide the required number of pieces.
%                   Thus the piece size is dictated by the overall duration of 
%                   the intervals, divided by nPieces.
%    ===========================================================================
%
% OUTPUTS
%
%    pieces         resulting intervals (split)
%    ids            for each final piece, the index of the original interval
%                   to which it belongs
%    =========================================================================
% Hint
%   to get overlapping windows, simply use:    sortrows([SplitIntervals(intervals, window); SplitIntervals(intervals+window/2, window)])
%
% EXAMPLE 1 (with pieceSize): 
% SplitIntervals([2 5; 6 10], 'pieceSize', 2) (i.e. split these intervals into 2-second pieces)
% gives intervals = [2 4; 6 8; 8 10] and ids = [1; 2; 2]
% EXAMPLE 2 (with nPieces): 
% SplitIntervals([2 5; 6 10], 'nPieces', 2) (i.e. split each of these intervals into 2)
% gives intervals = [2 3.5; 3.5 5; 6 8; 8 0] and ids = [1; 1; 2; 2]

% Copyright (C) 2016-2023 by Ralitsa Todorova, 2023 by Federica Lareno Faccini
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

pieceSize = 0.02;
nPieces = [];

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'piecesize',
            pieceSize = varargin{i+1};
            if ~isvector(pieceSize) || length(pieceSize) ~= 1,
                error('Incorrect value for property ''pieceSize'' (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).');
            end
        case 'npieces',
            nPieces = varargin{i+1};
            if ~isvector(nPieces) || length(nPieces) ~= 1,
                error('Incorrect value for property ''nPieces'' (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).']);
    end
end

%% SPLIT EACH INTERVAL INTO EQUAL PIECES (nPieces)
if ~isempty(nPieces),
    matrix = nan(size(intervals,1),nPieces+1);
    matrix(:,[1 end]) = intervals;
    for i=2:(nPieces),
        matrix(:,i) = matrix(:,1)+(i-1)*(matrix(:,end)-matrix(:,1))/nPieces;
    end
    pieces = [reshape(matrix(:,1:end-1)',[],1) reshape(matrix(:,2:end)',[],1)];
    ids = ceil((1:size(pieces,1))'/nPieces);
    return
end

%% SPLIT INTERVALS INTO PIECES OF EQUAL ABSULUTE LENGTHS (pieceSize)
d = diff(intervals,[],2);
piecesPerInterval = floor(d./pieceSize);
firstPiecePosition = CumSum([1;piecesPerInterval(1:end-1)]);
% create the pieces by filling in the interval starts in their respective positions
pieces = zeros(sum(piecesPerInterval),1);
% --Modification by Federica Lareno Faccini (August 2023)
% Remove the last interval if it's smaller than pieceSize (piecesPerInterval(end)==0)
firstPiecePosition(firstPiecePosition>sum(piecesPerInterval)) = []; firstPiecePosition = unique(firstPiecePosition);
pieces(firstPiecePosition,1) = intervals(piecesPerInterval>0,1);
% --End of modification
%pieces(firstPiecePosition,1) = intervals(:,1);
reset = pieces>0; reset(1)=1;
time2add = ones(size(pieces))*pieceSize;
time2add = CumSum(time2add,reset) - pieceSize; %reset at new bin starts and make sure baseline is 0
if pieceSize<1, time2add = round(time2add/pieceSize)*pieceSize; end %fix underflow errers: time2add should be divisible by pieceSize
pieces = CumSum(pieces,reset) + time2add; %repeat first pieces for the duration of the interval
pieces(:,2) = pieces(:,1)+pieceSize;
ids = CumSum(reset);

% Block added by Federica Lareno Faccini (August 2023)
% To take into account the index of possible intervals too short to be included
% Without this, there was a shift of index number due to the exclusion of too short intervals
piecesTooShort = find(piecesPerInterval==0);
if ~isempty(piecesTooShort),
    if piecesTooShort(1)==1, % In case the first n intervals have been skipped because too short,
                             % remove them from the list and directly add them to the ids (so that it starts from the right value)
        missingHead = find(diff(piecesTooShort)>1, 1, 'first');
        ids = ids+missingHead; 
        piecesTooShort(1:missingHead)=[];
    end
    if length(piecesTooShort) == size(intervals,1)-1,
        ids = size(intervals,1);
        warning('All intervals except last one have been excluded because shorter than the pieceSize selected!')
        return
    end
    if piecesTooShort(end) == size(intervals,1), piecesTooShort(end)=[]; end% Remove last index, avoid error at arrayfun
end
if ~isempty(piecesTooShort),
    lostIndexes = arrayfun(@(x) find(ids>=x,1,'first'), piecesTooShort);
    for i=1:size(lostIndexes,1), ids(lostIndexes(i):end) = ids(lostIndexes(i):end)+1; end
end
end


