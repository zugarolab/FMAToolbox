function DBAddVariable(v,eid,name,comments,parameters,mfiles)

%DBAddVariable - Save variable in current database.
%
% Save variable and code used to compute it.
%
%  USAGE
%
%    DBAddVariable(v,eid,name,comments,parameters,mfiles)
%
%    v              variable
%    eid            experiment ID (identifier string)
%    name           descriptive name
%    comments       comment string
%    parameters     parameters
%    mfiles         cell array of m-file names (used to compute the variable)
%
%  NOTE
%
%    The pair (eid,name) must uniquely define the variable (they are used as
%    primary keys in the database).
%
%  DEFAULTS
%
%    Use <a href="matlab:help DBSettings">DBSettings</a> to check the maximum lengths for fields eid, name, etc.
%    These defaults can be changed. See <a href="matlab:help Database">Database</a> for details.
%
%    By default, the data are stored in the file system, and only a reference
%    is held in the database, but for special use cases it is possible to
%    change this behavior. See <a href="matlab:help Database">Database</a> for details.
%
%  SEE
%
%    See also DBAddFigure, DBGetVariables.
%

% Copyright (C) 2007-2023 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global dbUser;
settings = DBSettings;
maxLength = settings.fields;
maxSize = settings.storage.maxInternalSize;

% Make sure MyM is installed and functional
CheckMyM;

% Check number of parameters
if nargin < 6,
	error('Incorrect number of parameters (type ''help <a href="matlab:help DBAddVariable">DBAddVariable</a>'' for details).');
end

% Check parameters
if ~ischar(eid),
	error('Incorrect EID string (type ''help <a href="matlab:help DBAddVariable">DBAddVariable</a>'' for details).');
end
if length(eid) > maxLength.eid,
	error('EID too long (type ''help <a href="matlab:help DBAddVariable">DBAddVariable</a>'' for details).');
end

if ~ischar(name),
	error('Incorrect figure name (type ''help <a href="matlab:help DBAddVariable">DBAddVariable</a>'' for details).');
end
if length(name) > maxLength.name,
	error('Name too long (type ''help <a href="matlab:help DBAddVariable">DBAddVariable</a>'' for details).');
end

if ~isempty(strfind(eid,'/')) || ~isempty(strfind(name,'/')),
        error('EIDs and names should not contain ''/'' (type ''help <a href="matlab:help DBAddVariable">DBAddVariable</a>'' for details).');
end

if ~ischar(parameters),
	error('Incorrect figure parameters (type ''help <a href="matlab:help DBAddVariable">DBAddVariable</a>'' for details).');
end
if length(parameters) > maxLength.parameters,
	error('Parameters too long (type ''help <a href="matlab:help DBAddVariable">DBAddVariable</a>'' for details).');
end
if ~ischar(comments),
	error('Incorrect comment string (type ''help <a href="matlab:help DBAddVariable">DBAddVariable</a>'' for details).');
end
if length(comments) > maxLength.comments,
	error('Comments too long (type ''help <a href="matlab:help DBAddVariable">DBAddVariable</a>'' for details).');
end
if isempty(comments),
	comments = '-';
end
if ~ischar(mfiles) && ~iscell(mfiles),
	error('Incorrect list of m-files (type ''help <a href="matlab:help DBAddVariable">DBAddVariable</a>'' for details).');
end
if ischar(mfiles),
	mfiles = {mfiles};
end

% M-Files
code = '';
for i = 1:length(mfiles),
	mfileName = which(mfiles{i});
	if isempty(mfileName),
		error(['M-File for function ''' mfiles{i} ''' not found.']);
	end
	fid = fopen(mfileName);
	code{i} = fread(fid);
	fclose(fid);
end

% Date and md5
d = datestr(now);
try
	md5 = CalcMD5(v);
catch
	warning('Could not compute MD5 (works only for numeric arrays)');
	md5 = 0;
end

% Check variable size: if too large (>0B by default), save to file and store path in database
s = whos('v');
if s.bytes > maxSize,
	% store outside DB
	storage = DBExternalStoragePath;
	targetDirectory = [storage '/' DBUse];
	if ~isdir(targetDirectory), mkdir(targetDirectory); end
	if ~isdir([targetDirectory '/variables']), mkdir([targetDirectory '/variables']); end
	matFile = [targetDirectory '/variables/' eid '-' name '.mat'];
	h = mym(['insert into variables (eid,name,comments,parameters,mfiles,code,date,user,md5,v) values ("{S}","{S}","{S}","{S}","{M}","{M}","{S}","{S}","{Si}","{M}")'],eid,name,comments,parameters,mfiles,code,d,dbUser,md5,matFile);
	save(matFile,'v');
else
	% store in DB
	h = mym(['insert into variables (eid,name,comments,parameters,mfiles,code,date,user,md5,v) values ("{S}","{S}","{S}","{S}","{M}","{M}","{S}","{S}","{Si}","{M}")'],eid,name,comments,parameters,mfiles,code,d,dbUser,md5,v);
end
