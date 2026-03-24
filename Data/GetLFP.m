function [lfp,indices] = GetLFP(channels,opt)

%GetLFP - Get local field potentials.
%
%  Load local field potentials from disk (unlike spikes or positions, LFP data
%  is usually too large to keep in memory).
%
%  USAGE
%
%    [lfp,indices] = GetLFP(channels,<options>)
%
%    channels       list of channels to load (use keyword 'all' for all)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'restrict'    list of time intervals to read from the LFP file
%     'select'      select channel by ID ('id', counted from 0 a la NeuroScope)
%                   or by number ('number', counted from 1 a la Matlab)
%                   (default = 'id')
%     'chunkSize'   number of samples to read in a single chunk
%                   (default = 10000)
%    =========================================================================
%
%  OUTPUT
%
%    lfp            list of (time,voltage1,...,voltageN) tuples
%    indices        for each tuple, the index of the interval it falls in
%
%  EXAMPLES
%
%    % channel ID 5 (= # 6), from 0 to 120 seconds
%    lfp = GetLFP(5,'restrict',[0 120]);
%    % same, plus from 240.2 to 265.23 seconds
%    lfp = GetLFP(5,'restrict',[0 120;240.2 265.23]);
%    % multiple channels
%    lfp = GetLFP([1 2 3 4 10 17],'restrict',[0 120]);
%    % channel # 3 (= ID 2), from 0 to 120 seconds
%    lfp = GetLFP(3,'restrict',[0 120],'select','number');

% Copyright (C) 2004-2022 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    channels
    opt.restrict (:,2) = []
    opt.intervals (:,2) = [] % retrocompatibility: replaced by 'restrict'
    opt.select (1,1) string {mustBeMember(opt.select,["id","number"])} = "id"
    opt.chunkSize (1,1) {mustBeInteger,mustBePositive} = 10000
end

% 'restrict' has priority over 'intervals'
if isempty(opt.restrict)
    opt.restrict = opt.intervals;
end
if isempty(opt.restrict)
    opt.restrict  = [0 Inf];
end

global DATA;
if isempty(DATA)
    error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end
filename = [DATA.session.path '/' DATA.session.basename '.lfp'];
if ~exist(filename,'file')
    error(['File ''' filename ''' not found.']);
end

nChannels = DATA.nChannels;
if isText(channels,'scalar',true) && strcmpi(channels,'all')
    channels = (1 : nChannels) - 1;
end
if opt.select == "id"
	channels = channels + 1;
end

nIntervals = size(opt.restrict,1);
lfp = cell(nIntervals);
indices = cell(nIntervals);
for i = 1 : nIntervals
	duration = (opt.restrict(i,2)-opt.restrict(i,1));
	start = opt.restrict(i,1);
	% Load data
	data = LoadBinary(filename,'duration',duration,'frequency',DATA.rates.lfp,'nchannels',nChannels,'start',start,'channels',channels,'chunkSize',opt.chunkSize);
	t = start : (1/DATA.rates.lfp) : (start+(length(data)-1)/DATA.rates.lfp);
	lfp{i} = [t.' data];
	indices{i} = i*ones(numel(t),1);
end
lfp = vertcat(lfp{:});
indices = vertcat(indices{:});