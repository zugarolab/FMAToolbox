function handle = RasterPlot(spikes,varargin,opt)

%RasterPlot - Create raster plot from spike times and IDs
%
%  USAGE
%
%    handle = RasterPlot(spikes,height,varargin,<NV args>)
%
%    spikes         list of [timestamps id]
%    height         optional spike height, default is 1
%    varargin       all other arguments (except Name-Value ones) will be passed to plot
%    <NV args>      optional list of Name-Value arguments (see table below)
%
%    =========================================================================
%     Name          Values
%    -------------------------------------------------------------------------
%     'label'       legend label for spikes, default is default plot behavior
%                   (e.g., data1), specify '' not to add spikes to legend
%     'ax'          ax to plot on, default is default plot behavior
%    =========================================================================
%
%  OUTPUT
%
%    handle         Line handle
%
%  EXAMPLES
%
%    % plot spikes with Line properties and label
%    h = RasterPlot(spikes,'LineWidth',1.3,'Marker','*',label='spikes');
%
%    % plot spikes with given height, without adding them to legend
%    h = RasterPlot(spikes,1.5,label='');
%
%    % plot spikes with empty legend label, on existing axes
%    h = RasterPlot(spikes,label=' ',ax=fistAx);

% Copyright (C) 2018-2022 by Ralitsa Todorova & (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    spikes (:,2) {mustBeNumeric}
end
arguments (Repeating)
    varargin
end
arguments
    opt.label = missing
    opt.ax (1,1) matlab.graphics.axis.Axes = gca
end

% validate label
if ~all(ismissing(opt.label),'all') && ~isText(opt.label)
    error('Argument ''label'' must be string or missing.')
end

% default value
if nargin == 1 || ~isnumeric(varargin{1}) || ~isscalar(varargin{1})
    height = 1;
else
    height = varargin{1};
    varargin = varargin(2:end);
end

% validate Line properties
if mod(numel(varargin),2) ~= 0
    error('Invalid Line properties (type ''help <a href="matlab:help plot">plot</a>'' for details).')
end

times = spikes(:,1);
rows = spikes(:,2);

times = [times times nan(size(times))].';
rows =  [rows-0.45*height rows+0.45*height nan(size(rows))].';

if ismissing(opt.label)
    % plot with default label (e.g., data1)
    handle = plot(opt.ax,times(:),rows(:),'LineWidth',1,varargin{:});
else
    % plot with given label
    handle = plot(opt.ax,times(:),rows(:),'LineWidth',1,'DisplayName',opt.label,varargin{:});
end
if opt.label == ""
    % plot without adding line to legend
    RemoveFromLegend(handle)
end