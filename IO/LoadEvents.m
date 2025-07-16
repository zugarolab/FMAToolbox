function events = LoadEvents(filename,opt)

%LoadEvents - Read events from file.
%
%  USAGE
%
%    events = LoadEvents(filename)
%
%    filename       event file name
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'compact'     if false (default), output is DESCRIBE; otherwise it is
%                   DESCRIBE
%    =========================================================================

% Copyright (C) 2004-2015 by Michaël Zugaro & (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    filename (1,:) char
    opt.compact (1,1) {mustBeLogical} = false
end

if ~isfile(filename)
	error(['File ''' filename ''' not found.']);
end

[~,~,extension] = fileparts(filename);
if extension ~= ".mat"
    % Read file into cell array
    file = fopen(filename,'r');
    if file == -1
	    error(['Cannot read ' filename ' (insufficient access rights?).']);
    end
    c = textscan(file,'%s','delimiter','\n');
    fclose(file);

    % parse cell array (extract time and messages using regular expressions)
    t = regexprep(c{1},'([^ \t]*).*','$1','once');
    events.time = cellfun(@(x) str2double(x)/1000,t); % convert to seconds
    events.description = regexprep(c{1},'[^ \t]*[ \t]*','','once');

    if opt.compact
        if any(cellfun(@(x) lower(x(1:9)),events.description(1:2:end),'UniformOutput',false) ~= "beginning") || ...
           any(cellfun(@(x) lower(x(1:3)),events.description(2:2:end),'UniformOutput',false) ~= "end")
            error('loadEvents:fileFormat','Wrong events file formnat.')
        end
        names = string(cellfun(@prepareName,events.description,'UniformOutput',false));
        if any(names(1:2:end) ~= names(2:2:end))
            error('loadEvents:fileFormat','Wrong events file formnat.')
        end
        names = names(1:2:end);

        % fuse intervals which are a repetition
        repetition = names(1:end-1) == names(2:end);
        times = {[events.time(1),events.time(2)]};
        for i = 2 : numel(names)
            if repetition(i-1)
                times{end} = [times{end};events.time(2*i-1),events.time(2*i)];
            else
                times{end+1,1} = [events.time(2*i-1),events.time(2*i)];
            end
        end
        events.description = names([true;~repetition]);
        events.time = times;
    end
else
    % Load CellExplorer file
    CE = load(filename);
    if isscalar(CE) && isscalar(fieldnames(CE)) % if only one variable was saved, use it
        CE = CE.(string(fieldnames(CE)));
    end
    if ~isfield(CE,'timestamps')
        error('Expecting field ''timestamps'' in .mat structure.')
    end
    events.time = CE.timestamps;
    % name of event
    name = extractBefore(filename,'.events.mat');
    
    if opt.compact
        % return compact output
        events.time = {events.time};
        events.description = name;
    else
        % standardize format
        % event name
        name = name(find(name == '.',1,'last')+1:end);
        name(1) = upper(name(1));
        if name(end) == 's'
            name = name(1:end-1);
        end

        % description of every time stamp
        events.description = repmat(string(name),size(events.time,2),1);
        if size(events.time,2) == 2
            events.description = events.description + [" start";" stop"];            
        end

        % add peak if present
        if size(events.time,2) == 2 && isfield(CE,'peaks')
            events.time = [events.time(:,1),CE.peaks,events.time(:,1)];
            events.description = [events.description(1);name+" peak";events.description(2)];
        end
        events.description = cellstr(repmat(events.description,size(events.time,1),1));
        events.time = events.time.';
        events.time = events.time(:);
    end
end
end

function name = prepareName(str)
    % words, first two are 'beginning of' or 'end of' and can be discarded
    str = strsplit(str);
    % parts of event name
    str = strsplit(str{end},'_');
    if all(isstrprop(str{end},'digit'))
        % remove last number
        str = str(1:end-1);
    end
    % remove 'bis' to correctly detect duplicates
    str{end} = erase(str{end},'bis');
    name = [str{:}];
end