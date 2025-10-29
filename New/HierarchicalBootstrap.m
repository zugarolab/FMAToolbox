function [means,p] = HierarchicalBootstrap(data,varargin)

% HierarchicalBootstrap:
%
% The function implements a hierarchical bootstrap as described in
% Saravanan et al. (2021) https://doi.org/10.1101/819334 with the addition
% of a within-subject design bootstrap, where data from the same e.g.
% animal is sampled together for both conditions (rather than independently
% for each condition, which would ignore the paired design).
% The function will give you a p-value testing if the values of condition
% 1 are systematically higher than the values corresponding to condition 2,
% after controlling for all the in-between nesting variables.
% The output "means" is the bootstrapped average of the values of conditions
% 1 and 2 for each of the iterations.
%
%  USAGE
%
%    [means,p] = HierarchicalBootstrap(data)
%
%    data           The expected format is long format, with the first 'nColumns'
%                   columns being the values to be tested, and subsequent
%                   columns indicate nested groupings of the data, from
%                   deepest to highest (see EXAMPLE)
%    nColumns       Number of columns in which data is present (default = 1).
%                   Other columns are assumed to contain grouping values.
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'nIterations' how many iterations to perform, resampling from the data
%                   (default = 1000)
%     'fun'         the function handle of the function computing the average
%                   value per condition (default = @mean)
%     'groups'      Boolean indicating if the last column of "data" refers to
%                   the condition of interest (default = true). When true,
%                   the function treats the last column as the condition,
%                   bootstrapping the mean separately for each condition
%                   while using the preceding columns as grouping
%                   variables. Set to false if all data comes from a single
%                   condition, in which case the function will assumes the last
%                   column is the highest-level grouping variable rather
%                   than a condition.
%     'output'      Set to 'cell' if the desired output is more than a
%                   single number and can therefore not be fit in a matrix
%                   (default = 'double').
%    =========================================================================
%
%  OUTPUT
%
%    means          bootstrapped results for the mean of each condition (one for
%                   each iteration)
%    p              In the special case when two conditions are provided, the
%                   p-value of the test that the values of condition 1 are
%                   higher than the values of condition 2 after controlling
%                   for the other nesting variables. When the number of
%                   conditions provided is not 2, the p-value will instead
%                   compare each of the bootstapped means to zero.
%
% EXAMPLE
%
% data = [fieldSize, cellID, session, animal, condition]
% [means,p] = HierarchicalBootstrap(data);
% "p<0.05" indicates that the fieldSizes of condition 1 are significantly
% higher than the fieldSizes of condition 2, after controlling for
% cellID, session, and animal.
% pValues = mean(means<=0); % is data from the each group significantly higher than zero
%
% Copyright (C) 2023-2025 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if iscell(data)
    % find out how deep the cell nesting goes:
    % Somehow group them into a matrix:
    error('Grouping of nested data in cells not yet implemented. Please provide matrix format');
end

% if any(~ismember(data(:,end),[1 2]))
%     error('The last column should provide the grouping variable (1 or 2)');
% end

nIterations = 1000;
average = @nanmean; % change to nanmedian for testing if the group medians are different
nColumns = 1;
grouped = true;
output = 'double';

% Parse parameters
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        builtin('error',['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
    end
    switch(lower(varargin{i}))
        case {'niterations','nshuffles','n'}
            nIterations = varargin{i+1};
            if ~(isscalar(nIterations) && (nIterations>0))
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
            end
        case {'fun','function','handle','average'}
            average = varargin{i+1};
            if ~isa(average,'function_handle')
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
            end
        case {'ncolumns'}
            nColumns = varargin{i+1};
            if ~(isscalar(nColumns) && (nColumns>0) && nColumns<size(data,2))
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
            end
        case {'output','outputs'}
            output = varargin{i+1};
            if ~isastring(lower(output),'cell','double')
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
            end
        case {'groups','group','grouped'}
            grouped = varargin{i+1};
            if isastring(lower(grouped),'on','off')
                grouped = strcmpi(grouped,'on'); % transform to logical
            end
            if ~(islogical(grouped) && length(grouped)==1)
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
            end
        otherwise
            builtin('error',['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
    end
end

if ~grouped % "data" does not contain a condition column
    data(:,end+1) = 1; % add a "condition" in the last column of data.
end

nGroups = max(data(:,end));

if strcmp(output,'cell')
    means = cell(nIterations, nGroups);
else
    means = nan(nIterations, nGroups);
end

% For each column, check if its indices appear in multiple conditions
% This is so that a paired boostrap can be performed if appropriate.
% Example: in a paired bootstrap, if in a given shuffle the bootstrap
% ignores animal2, then the paired bootstrap will ignore it for all
% conditions. This is undesirable if data are not paired: i.e. if
% observations for animal2 only apply to one of the conditions.
levels = nColumns+1:size(data,2)-1;

% Remove empty groups (sessions, animals):
for k=levels
    [~,~,data(:,k)] = unique(data(:,k));
end

pairedBootstrap = false(size(levels));
for i=1:length(levels)
    level=levels(i);
    percentUnique = mean(min(Accumulate(data(:,[level end])),[],2)==0);
    % If most of the cells appear in more than one condition, do a paired boostrap
    if percentUnique<0.5, pairedBootstrap(i) = true; end
end

% Make a nested of cells for easier sampling:
nested = MakeNestedGroups(data,nColumns,pairedBootstrap);

for k = 1:nIterations
    resampledNestedGroups = BootstrapNestedGroups(nested);
    resampled = cat(1, resampledNestedGroups{:});
    while iscell(resampled)
        resampled = cat(1,resampled{:});
    end
    
    for j=1:nGroups
        if strcmp(output,'cell')
            try
                means{k,j} = average(resampled(resampled(:,end)==j,1:nColumns));
            catch
                means{k,j} = nan;
            end
        else
            try
                means(k,j) = average(resampled(resampled(:,end)==j,1:nColumns));
            catch
                means(k,j) = nan;
            end
        end
    end
end

if ~strcmp(output,'cell') && nargout>1
    if size(means,2)==2
        p = mean(means(:,1)<means(:,2));
    else
        p = mean(means<0);
    end
end


% ------------------------------- Helper functions -------------------------------

function nested = MakeNestedGroups(data,nColumns,pairedBootstrap)

if iscell(data) % if it's already a cell, work through its elements
    nested = cell(size(data));
    for i=1:length(data)
        nested{i} = MakeNestedGroups(data{i},nColumns,pairedBootstrap);
    end
    return
end

if size(data,2)==nColumns+1 % data and the (last) condition column
    nested = data;
    return; % no need for further branching
end

if ~pairedBootstrap(end)
    % skip the empty bins (no need to take animal2 when considering
    % condition1 if there are no observations of animal2 in condition1).
    ids = unique(data(:,end-1));
else
    % we don't skip the empty bins if we want to keep the numbering for this
    % level because both conditions exist (so we want to resample from both
    % conditions together in a paired way; overriding the numbering by skipping
    % empty cells would scramble the pairing)
    ids = 1:max(data(:,end-1));
end
nested = cell(length(ids),1);

for j = 1:length(ids)
    nested{j} = data(data(:,end-1) == ids(j), [1:end-2 end]);
end

if size(data,2)>nColumns+2 % We need to call the function more than once:
    for i=1:length(nested)
        nested{i} = MakeNestedGroups(nested{i},nColumns,pairedBootstrap(1:end-1));
    end
    return
end

function resampledNestedGroups = BootstrapNestedGroups(nested)

if ~iscell(nested)
    % Sample with replacement the values within the cell
    n = size(nested,1);
    idx = ceil(rand(n,1)*n);
    resampledNestedGroups = nested(idx,:);
else
    % Sample the cells with replacement
    n = numel(nested);
    idx = ceil(rand(n,1)*n);
    resampledNestedGroups = nested(idx);
    for i = 1:n
        % Sample with replacement the values within each cell
        resampledNestedGroups{i} = BootstrapNestedGroups(resampledNestedGroups{i});
    end
end

