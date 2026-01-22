function [synchronized,Ie,Is] = Sync(samples,sync,varargin)

%Sync - Make sample timestamps relative to synchronizing events.
%
% Select samples that fall around synchronizing events, and make their
% timestamps relative to the synchronizing events. This can be used to
% build e.g. spike raster plots or successive evoked potentials.
%
%  USAGE
%
%    [synchronized,Ie,Is] = Sync(samples,sync,<options>)
%
%    samples        <a href="matlab:help samples">samples</a> to synchronize
%    sync           timestamps to synchronize on (e.g., brain stimulations)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'durations'   durations before and after synchronizing events for each
%                   trial (in s) (default = [-0.5 0.5])
%     'fast'        if 'off' (default), sort 'samples' and 'sync' before
%                   operating, otherwise they are expected to be sorted
%     'verbose'     display information about ongoing processing
%                   (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    synchronized   resynchronized samples (sorted in time)
%    Ie             index of original synchronizing event for each resynchronized sample
%    Is             indeces of resynchronized samples in original samples 
%
%  EXAMPLE
%
%    [raster,indices] = Sync(spikes,stimuli);     % compute spike raster data
%    figure; PlotSync(raster,indices);            % plot spike raster
%
%  SEE
%
%    See also SyncHist, SyncMap, PlotSync, PETHTransition.

% Copyright (C) 2004-2014 by Michaël Zugaro, 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
durations = [-0.5 0.5];
fast = false;
verbose = false;

% Check number of parameters
if nargin < 2 || mod(length(varargin),2) ~= 0
  error('Incorrect number of parameters (type ''help <a href="matlab:help Sync">Sync</a>'' for details).');
end

% Check parameter sizes
if ~isdvector(sync)
	error('Parameter ''sync'' is not a vector (type ''help <a href="matlab:help Sync">Sync</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
		    error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Sync">Sync</a>'' for details).']);
    end
    switch lower(varargin{i})
        case 'durations'
            durations = varargin{i+1};
            if ~isdvector(durations,'#2','<=')
                error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help Sync">Sync</a>'' for details).');
            end
        case 'fast'
            fast = varargin{i+1};
            if ~isastring(fast,'on','off')
                error('Incorrect value for property ''fast'' (type ''help <a href="matlab:help Sync">Sync</a>'' for details).');
            end
            fast = strcmp(fast,'on');
        case 'verbose'
            verbose = varargin{i+1};
            if ~isastring(verbose,'on','off')
                error('Incorrect value for property ''verbose'' (type ''help <a href="matlab:help Sync">Sync</a>'' for details).');
            end
            verbose = strcmp(verbose,'on');
        otherwise
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Sync">Sync</a>'' for details).']);
    end
end

% Make sync a column vector
if size(sync,2) ~= 1
    sync = sync.';
end

% Make sure samples and sync events are sorted in time
sort_samples = [];
sort_sync = [];
if ~fast
    [samples,sort_samples] = sortrows(samples);
    [sync,sort_sync] = sort(sync);
end

nSync = length(sync);
Is = cell(nSync,1);
previous = 1; % index of first sample found in previous cycle
verbose && fprintf(1,[num2str(nSync) ' synchronizing events to process...\n']);
for i = 1 : nSync
    % Find samples within time window around this synchronizing event
  	j = FindInInterval(samples(:,1),[sync(i)+durations(1) sync(i)+durations(2)],previous);
    if ~isempty(j)
        previous = j(1);
        Is{i} = (j(1) : j(2)).';
        verbose && fprintf(1,[' sync ' int2str(i) ' (t=' num2str(sync(i)) '): ' int2str(length(j)) ' samples\n']);
    end
end

Ie = repelem((1:nSync).',cellfun(@numel,Is),1);
Is = vertcat(Is{:});
synchronized = [samples(Is,1)-sync(Ie),samples(Is,2:end)];

% Make Is, Ie point to elements of original samples, sync
if ~isempty(sort_samples)
  Is = sort_samples(Is);
  Ie = sort_sync(Ie);
end