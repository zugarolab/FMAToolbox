function settings = DBSettings

%DBSettings - Get default settings for databases.
%
%  This is an 'internal' function used by FMAToolbox. You should not need
%  to use it, unless you are developping new functions for this toolbox.


% Copyright (C) 2023 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global SETTINGS;

% Defaults
settings.fields.eid = 50;
settings.fields.name = 100;
settings.fields.parameters = 50;
settings.fields.comments = 255;
settings.storage.maxInternalSize = 0;
settings.storage.externalPath = '';

% Check user settings

if isfield(SETTINGS,'dbExternalStoragePath'),
    % For backward compatibility
    settings.storage.externalPath = SETTINGS.dbExternalStoragePath;
end

% Database field lengths
if ~isfield(SETTINGS,'db'), return; end

if  isfield(SETTINGS.db,'fields'),
    if isfield(SETTINGS.db.fields,'eid'),
        settings.field.eid = SETTINGS.db.fields.eid;
    end
    if isfield(SETTINGS.db.fields,'name'),
        settings.fields.name = SETTINGS.db.fields.name;
    end
    if isfield(SETTINGS.db.fields,'parameters'),
        settings.fields.parameters = SETTINGS.db.fields.parameters;
    end
    if isfield(SETTINGS.db.fields,'comments'),
        settings.fields.comments = SETTINGS.db.fields.comments;
    end
end

% External storage
if isfield(SETTINGS.db,'storage'),
    if isfield(SETTINGS.db.storage,'maxInternalSize'),
        settings.storage.maxInternalSize = SETTINGS.db.storage.maxInternalSize;
    end
    if isfield(SETTINGS.db.storage,'externalPath'),
        settings.storage.externalPath = SETTINGS.db.storage.externalPath;
    end
end

