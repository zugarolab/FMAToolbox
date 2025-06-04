function varargout = semplot(varargin,opt)

%semplot - Plot mean (line) +/- s.e.m. (shaded area) of a matrix "y"
%
%  USAGE
%
%    handles = semplot(x,y,color,<options>)
%
%    x               abscissae (optional, default is 1 : size(y,2))
%    y               ordinates, each column contains multiple values for an
%                    element of x
%    color           line color (optional, default is black)
%    <options>       optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties     Values
%    -------------------------------------------------------------------------
%     'mode'         use either 'sem' (default) or 'std' (standard deviation)
%                    to draw shaded area
%     'smooth'       Gaussian kernel size in number of samples (default is none)
%     'solid'        if false (default), shaded area is transparent
%     'faceColor'    shaded area color, default is same as color
%     'legend'       if 'off', plotted elements won't appear in legend
%                    (default = 'on'); value is also legend label if
%                    different from 'on'
%     'lineProp'     cell array of property-value pairs to set line properties
%                    (see MATLAB Line Properties)
%     'patchProp'    cell array of property-value pairs to set shaded area
%                    properties (see MATLAB Patch Properties)
%    =========================================================================
%
%  OUTPUT
%
%    handles        graphics handle for shaded area Patch object (optional)
%
%  SEE
%
%    See also shadedErrorBar

% Copyright (C) 2017 by Ralitsa Todorova, 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments (Repeating)
    varargin
end
arguments
    opt.mode (1,1) string {mustBeMember(opt.mode,["sem","std"])} = 'sem'
    opt.smooth (1,1) {mustBeNumeric,mustBeNonnegative} = 0
    opt.solid (1,1) {mustBeMember(opt.solid,[0,1])} = false
    opt.faceColor (:,:) = NaN
    opt.legend (1,1) string = 'on'
    opt.lineProp (:,1) cell = {}
    opt.patchProp (:,1) cell = {}
end

% retrocompatibility: syntax semplot(x,y,color,smooth,solid) used to be accepted
% varargin is thus necessary
[x,y,color,opt.smooth,opt.solid,opt.faceColor] = parseSemplot(varargin,opt);

% convert colors like 'k' to RGB triplet
try
    opt.faceColor = validatecolor(opt.faceColor);
catch ME
    throw(ME)
end
% set default values
if opt.solid
    opt.color = mean([opt.faceColor;1 1 1]);
    faceAlpha = 1;
else
    faceAlpha = 0.5;
end

% if only one value per ascissa is given
if isvector(y)
    handles = plot(x,Smooth(y,opt.smooth),'color',color);
    if opt.legend == "off"
        RemoveFromLegend(handles)
    elseif opt.legend ~= "on"
        handles.DisplayName = opt.legend;
    end
    if nargout > 0, varargout{1} = handles; end
    hold on
    return
end

% remove ascissae where all values are NaN
y_mean = mean(y,'omitmissing').';
bad = isnan(y_mean);
x = x(~bad);
y = y(:,~bad);
y_mean = y_mean(~bad);

% prepare patch coordinates and smooth
xx = [x;flipud(x)];
if opt.mode == "sem"
    y_area = nansem(y);
else
    y_area = std(y,'omitmissing');
end
yy = [Smooth(y_mean-y_area.',opt.smooth); Smooth(flipud(y_mean+y_area.'),opt.smooth)];
y_mean = Smooth(y_mean,opt.smooth);

% plot shaded area
handles = fill(xx,yy,opt.faceColor,'EdgeAlpha',0,'FaceAlpha',faceAlpha,opt.patchProp{:});
if opt.legend ~= "on"
    handles.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

% plot line
hold on
h = plot(x,y_mean,'color',color,'linewidth',2,opt.lineProp{:});
if opt.legend == "off"
    RemoveFromLegend(h)
elseif opt.legend ~= "on"
    h.DisplayName = opt.legend;
end

if nargout > 0, varargout{1} = handles; end

end

% --- helper functions ---

function RemoveFromLegend(h)
    % this f is a copy of the private f RemoveFromLegend of the Plot
    % folder, it won't be necessary once semplot is moved there

    arguments
        h (1,:)
    end

    % set property of graphic objects
    for hand = h
        hand.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end

    % get parent figure
    figHandle = ancestor(h(1),'figure');

    % if legends exist, refresh them to effectively hide entries
    legHandle = findobj(figHandle,'Type','Legend');
    for hand = legHandle.'
        hand.Visible = 'off';
        hand.Visible = 'on';
    end
end

function [x,y,color,smooth,solid,faceColor] = parseSemplot(args,opt)
    % automatic handling of Name-Value Arguments removes from varargin all valid property-value pairs* and
    % stores them in the struct opt
    % this function checks the validity of remaining arguments and assigns default valus
    %
    % * the only case when this is not true is if the user specified an invalid property name

    narg = numel(args);

    % minimun accepted number of arguments is 1
    if narg < 1
        error('semplot:arguments','Invalid argument list. Function requires 1 more input.')
    end

    % set default values
    color = [0,0,0];
    smooth = opt.smooth;
    solid = opt.solid;
    faceColor = opt.faceColor;

    % parse arguments
    x = args{1};
    ind = 3; % index of next element of varargin to parse
    if narg < 2
        % set default x
        y = x;
        x = (1:size(y,2)).';
    elseif ischar(args{2}) || isstring(args{2})
        % set default x, use given char color
        color = args{2};
        y = x;
        x = (1:size(y,2)).';
    else
        if ~isvector(x)
            error('semplot:xSize','X must be a vector');
        end
        if size(x,2) ~= 1
            x = x.';
        end
        y = args{2};
        if size(y,2) ~= length(x)
            y = y.';
            if size(y,2) ~= length(x)
                error('Y must have one column for each element in X');
            end
        end
        if narg > 2
            % use given color
            color = args{3};
            ind = 4;
        end
    end

    if any(isnan(faceColor),'all')
        % NaN signals to match value of color
        faceColor = color;
    end

    % if any string (even a valid property) is left in varargin, one of the  properties was invalid
    bad_property = 0; % index of invalid property
    if narg >= ind
        if isastring(args{ind})
            bad_property = ind;
        elseif ~isscalar(args{ind}) || ~isnumeric(args{ind}) || args{ind} < 0
            error('semplot:smoothValue','Incorrect value for property ''smooth'' (type ''help <a href="matlab:help semplot">semplot</a>'' for details).')
        end
        opt.smooth = args{ind};
        ind = ind + 1;
    end
    if ~bad_property && narg >= ind
        if isastring(args{ind})
            bad_property = ind;
        elseif ~isscalar(args{ind}) || ~islogical(args{ind})
            error('semplot:solidValue','Incorrect value for property ''solid'' (type ''help <a href="matlab:help semplot">semplot</a>'' for details).')
        end
        opt.solid = args{ind};
        ind = ind + 1;
    end

    % if an invalid property was detected, error accordingly
    % the following doesn't error if the bad property has a valid name (e.g., a valid property written twice), in which case the last control will error
    if bad_property
        for i = bad_property : 2 : numel(args)
            if ~ismember(args{i},fieldnames(opt))
                error('semplot:unknownArgument',"Unknown property '"+args{i}+"' (type 'help <a href=""matlab:help semplot"">semplot</a>' for details).");
            end
        end
    end

    % if any argument is left, it must be an invalid property
    if bad_property || narg >= ind
        error('semplot:unknownArgument',"Invalid value in position "+num2str(min(ind,bad_property))+" (type 'help <a href=""matlab:help semplot"">semplot</a>' for details).");
    end
end