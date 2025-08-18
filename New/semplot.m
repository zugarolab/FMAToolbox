function varargout = semplot(x,y,color,solid,smooth,opt)

%semplot - Plot mean (line) +/- s.e.m. (shaded area) of a matrix 'y'
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

arguments
    x
    y = []
    color = [0,0,0]
    solid = false
    smooth = 0
    opt.mode (1,1) string {mustBeMember(opt.mode,["sem","std"])} = "sem"
    opt.smooth {mustBeScalarOrEmpty} = []
    opt.solid {mustBeLogicalScalarOrEmpty} = []
    opt.faceColor (:,:) = NaN
    opt.legend (1,1) string = "on"
    opt.lineProp (:,1) cell = {}
    opt.patchProp (:,1) cell = {}
end

% retrocompatibility: syntax semplot(x,y,color,smooth,solid) used to be accepted
try
    [x,y,color,solid,smooth,opt] = parseSemPlot(x,y,color,solid,smooth,nargin,opt);
catch ME
    throw(ME)
end

% set default values
if solid
    opt.faceColor = mean([opt.faceColor;1 1 1]);
    faceAlpha = 1;
else
    faceAlpha = 0.5;
end

% if only one value per ascissa is given
if isvector(y)
    try
        if opt.charColor
            handles = plot(x,Smooth(y,smooth),color,'linewidth',2,opt.lineProp{:});
        else
            handles = plot(x,Smooth(y,smooth),'color',color,'linewidth',2,opt.lineProp{:});
        end
    catch ME
        % catch invalid color
        throw(ME)
    end
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
y_mean = mean(y,'omitnan').';
bad = isnan(y_mean);
x = x(~bad);
y = y(:,~bad);
y_mean = y_mean(~bad);

% prepare patch coordinates and smooth
xx = [x;flipud(x)];
if opt.mode == "sem"
    y_area = nansem(y);
else
    y_area = std(y,'omitnan');
end
yy = [Smooth(y_mean-y_area.',smooth); Smooth(flipud(y_mean+y_area.'),smooth)];
y_mean = Smooth(y_mean,smooth);

% plot shaded area
try
    handles = fill(xx,yy,opt.faceColor,'EdgeAlpha',0,'FaceAlpha',faceAlpha,opt.patchProp{:});
catch ME
    % catch invalid color
    throw(ME)
end
if opt.legend ~= "on"
    handles.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

% plot line
hold on
try
    if opt.charColor
        h = plot(x,y_mean,color,'linewidth',2,opt.lineProp{:});
    else
        h = plot(x,y_mean,'color',color,'linewidth',2,opt.lineProp{:});
    end
catch ME
    % catch invalid color
    throw(ME)
end
if opt.legend == "off"
    RemoveFromLegend(h)
elseif opt.legend ~= "on"
    h.DisplayName = opt.legend;
end

if nargout > 0, varargout{1} = handles; end

end

% --- helper functions ---

function [x,y,color,solid,smooth,opt] = parseSemPlot(x,y,color,solid,smooth,n,opt)

% this f identifies the used syntax and validates inputs
% allowed syntaxes:
%   semplot(y,<color>,<solid>,<smooth>,<opt>) % only if color is text (like 'r')
%   semplot(x,y,<color>,<solid>,<smooth>,<opt>)
% specifying  'solid',<value>  always overrides  <solid>, as for 'smooth'

if n < 2
    % syntax: semplot(y)
    % set default x
    y = x;
    x = (1 : size(y,2)).';
elseif ischar(y) || isstring(y)
    % syntax: semplot(y,color,...)
    % set default x, use given char color
    if n > 4
        error('semplot:nArgs','1 to 4 positional arguments allowed when omitting argument ''x''')
    end
    if n > 3, smooth = solid; end
    if n > 2, solid = color; end
    color = y;
    y = x;
    x = (1 : size(y,2)).';
else
    % syntax: semplot(x,y,...)
    if ~isvector(x)
        error('semplot:xSize','Argument ''x'' must be a vector');
    end
    if size(x,2) ~= 1
        x = x.';
    end
    if size(y,2) ~= length(x)
        y = y.';
        if size(y,2) ~= length(x)
            error('semplot:ySize','Argument ''y'' must have one column for each element in ''x''');
        end
    end
end

% syntax of plot is different if color is a LineSpec (as '-.r'), this flag allows to call plot correctly
opt.charColor = ischar(color) || isstring(color);

% property-value pairs for 'solid' and 'smooth' have precedence over positional arguments
if ~isempty(opt.solid)
    solid = opt.solid;
end
try
    solid = GeneralLogical(solid);
catch
    error('Invalid value for argument ''solid''. Value must be logical.')
end
if ~isempty(opt.smooth)
    smooth = opt.smooth;
end
if ~isscalar(smooth) || ~isnumeric(smooth) || smooth < 0
    error('Invalid value for argument ''smooth''. Value must be non-negative scalar.')
end

if any(isnan(opt.faceColor),'all')
    % NaN signals to match value of 'color'
    opt.faceColor = color;
end

end