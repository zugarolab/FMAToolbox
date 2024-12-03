function [indices,values] = FindClosest(reference, query, mode)

% This function looks up a query vector in a reference vector, and for each
% value in the query vector, returns the index of the closest value in the
% reference vector. The equivalent non-optimised code is trivial:
% for i=1:length(reference), indices(i,1) = find(abs(query-reference(i)) == min(abs(query-reference(i))), 1, 'first'); end
% EXAMPLE 1 (reasoning):
% FindClosest([0.5; 0.2],[0.3 0.45 0.66]) returns [2;1]
% as the closest value to 0.5 was 0.45 (index 2), and to 0.2 was 0.3 (index 1)
% EXAMPLE 2 (usage):
% indices = FindClosests(spikes1,spikes2) will return the indices of spikes2 closest to spikes1
% spikes1 - spikes2(indices) will contain all the minimal distances between spikes1 and spikes2 (one distance per each spike in spike1)
%
% (C) 2016 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
if ~exist('mode','var'),mode = 'either'; end

if isempty(reference) || isempty(query), indices = []; values = []; return; end
if length(reference)==1,indices = 1; values=reference; return; end

% Make sure values to be matched are unique
[u,i] = unique(reference); nans = isnan(u); u(nans) = []; i(nans) =[];
% Find closest index
if strcmp(mode,'higher'),
    indices = ceil(interp1(u,(1:length(u))',query));
elseif strcmp(mode,'lower'),
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
