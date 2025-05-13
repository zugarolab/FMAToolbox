function times = LoadSpikeTimes(filename,rate)

%LoadSpikeTimes - Load spike times and cluster IDs from disk.
%
%  USAGE
%
%    times = LoadSpikeTimes(filename,rate)
%
%    filename            spike file name (either .clu or .res)
%    rate                sampling rate
%
%  OUTPUT
%
%    The output is a list of (timestamp,group,cluster) t-uples.
%
%  SEE
%
%    See also GetSpikeTimes, PlotTicks.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

[path,basename,extension] = fileparts(filename);
if isempty(path), path = '.'; end

electrodeGroup = str2num(extension(2:end));
[~,basename] = fileparts(basename);

% Load .res file
filename = [path '/' basename '.res.' int2str(electrodeGroup)];
if ~exist(filename,'file')
	error('LoadSpikeTimes:MissingRes',['File ''' filename ''' not found.']);
end
res = readmatrix(filename,FileType='text',CommentStyle='%');

% Load .clu file
filename = [path '/' basename '.clu.' int2str(electrodeGroup)];
if ~exist(filename,'file')
	error('LoadSpikeTimes:MissingClu',['File ''' filename ''' not found.']);
end
clu = readmatrix(filename,FileType='text',CommentStyle='%');

if size(res,1) ~= size(clu,1)-1
    error('LoadSpikeTimes:ResCluSize','Contents of .res and .clu files have incompatible size.');
end

times = [res/rate electrodeGroup*ones(size(res)) clu(2:end)];