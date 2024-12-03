function f = DBGetFigures(query,varargin)

%DBGetFigures - Get all figures that match given criteria.
%
% Open figures and get related information such as code used to generate them.
%
%  USAGE
%
%    f = DBGetFigures(query,<options>)
%    f = DBGetFigures(figures,<options>)
%
%    query          optional figure list query (WHERE clause; see Example 1)
%    figures        optional file name (or cell array of file names) of a PNG
%                   image from an exported gallery (in the 'images' folder)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'output'       'figures' to get the figures but not the information,
%                    i.e. eid, name, comments, parameters, etc. (default),
%                    'info' for the information but not the figures, 'full'
%                    for both, or 'keys' for eid and name only
%    =========================================================================
%
%  OUTPUT
%
%    f is a structure with the following fields:
%
%    eid            experiment ID (identifier string)
%    name           figure descriptive name (identifier string)
%    fig            figure handle
%    comments       comments
%    parameters     figure parameters
%    code           text of m-files used to generate the figure
%    date           date when the figure was saved
%    user           the user that saved the figure
%
%    Each field is a cell array.
%
%  EXAMPLE 1
%
%    Get all figures from experiment "experiment1", the name of which starts
%    with "raster" (for details, see an SQL manual):
%
%    [f,comments,parameters,code,date,user] = ...
%      DBGetFigures('eid="experiment1" and name like "raster%"');
%
%  EXAMPLE 2
%
%    You have two image files named 'spectrogram.png' and 'correlation.png',
%    previously exported using DBExportGallery, but you cannot remember from which
%    experiment they were generated. Still, you would like to be able to edit the
%    original Matlab figures. This can be done with:
%
%    DBGetFigures({'spectrogram.png','correlation.png'});
%
%  SEE
%
%    See also DBAddFigure, DBGetVariables, DBDisplay.
%

% Copyright (C) 2007-2023 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

% Default values
output = 'figures';

% Optional query provided?
if nargin == 0,
	query = '';
elseif isastring(query,'output'),
	varargin = {query varargin{:}};
	query = '';
end

% Edit query
if ischar(query) && exist(query,'file'),
	query = {query};
end
if iscell(query),
	files = query;
	for i = 1:length(files),
		md5 = CalcMD5(files{i},'File');
		if i == 1,
			query = [' where md5="' md5 '"'];
		else
			query = [query ' OR md5="' md5 '"'];
		end
	end
else
	query = strtrim(query);
	query = regexprep(query,'^where','');
	if ~isempty(query), query = [' where ' query]; end
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help DBGetFigures">DBGetFigures</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'output',
			output = lower(varargin{i+1});
			if ~isastring(output,'figures','info','full','keys'),
				error('Incorrect value for property ''output'' (type ''help <a href="matlab:help DBGetFigures">DBGetFigures</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help DBGetFigures">DBGetFigures</a>'' for details).']);
	end
end

% Query database
switch output,
	case 'full',
		f = mym(['select fig,eid,name,comments,parameters,mfiles,code,date,user from figures' query]);
	case 'figures',
		f = mym(['select fig from figures' query]);
	case 'info',
		f = mym(['select eid,name,comments,parameters,mfiles,code,date,user from figures' query]);
	case 'keys',
		f = mym(['select eid,name from figures' query]);
end

% Make sure query results are not empty
if isempty(f),
	warning(['No figures match (' query ').']);
end

% Figures
if strcmp(output,'full') || strcmp(output,'figures'),

	for i = 1:length(f.fig),
		if isempty(f.fig{i}),
			warning(['Figure was stored as PNG only.']);
			f.fig{i} = [];
			continue;
		end
		blob = f.fig{i};
		if isastring(blob) && length(blob) > 4 && blob(1) == '/' && strcmp('.fig',blob(end-3:end)),
			% Blob is a path: load figure from external storage
			figName = blob;
			if ~exist(figName,'file'),
				error(['External storage file for (' x.eid{i} ',' x.name{i} ') is missing.']);
			else
				f.fig{i} = openfig(figName);
				set(f.fig{i},'visible','on');
			end
		else
			% Blob is the figure data: save to a temporary fig file, open it and delete it
			basename = tempname;
			figName = [basename '.fig'];
			file = fopen(figName,'wb');
			if file == -1,
				error(['Could not create temporary file for figure (' f.eid{i} ',' f.name{i} ').']);
			end
			fwrite(file,blob);
			fclose(file);
			f.fig{i} = openfig(figName);
			delete(figName);
			set(f.fig{i},'visible','on');
		end
	end

end

% Format code
if strcmp(output,'full') || strcmp(output,'info'),
	code = {};
	for i = 1:length(f.code),
		for j = 1:length(f.code{i}),
			code{i}{j} = char(f.code{i}{j})';
		end
	end
	f.code = code;
end
