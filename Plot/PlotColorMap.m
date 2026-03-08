function PlotColorMap(data,dimm,opt)

%PlotColorMap - Plot a color map.
%
%  Plot a color map (e.g. the firing field of a place cell).
%
%  USAGE
%
%    PlotColorMap(data,dimm,<options>)
%
%    data           data
%    dimm           optional luminance map
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'x'           abscissae
%     'y'           ordinates
%     'threshold'   dimm values below this limit are zeroed (default = 0.01)
%     'cutoffs'     lower and upper cutoff values ([] = autoscale, default; NaNs
%                   are also replaced by default values)
%     'hgamma'      gamma-like correction for hue (1 = no correction, default)
%     'bar'         draw a color bar (default = 'off'); if value isn't
%                   'on' it's used as label for the bar
%     'barProp'     cell array of ColorBar properties passed to the call of colorbar
%     'type'        either 'linear' or 'circular' (default 'linear')
%     'map'         colormap (default = Bright(100,'hgamma',hgamma,'type',type))
%     'ydir'        either 'normal' (default) or 'reverse' (useful when the
%                   x and y coordinates correspond to spatial positions,
%                   as video cameras measure y in reverse direction)
%     'piecewise'   if 'on' (default), set piecewise linear axis labels
%     'ax'          axes to plot on (default = gca)
%    =========================================================================
%
%  NOTE
%
%    The luminance map is used to dimm the color map. A single scalar value
%    is interpreted as a constant luminance map. If this parameter is not
%    provided, normal equiluminance is assumed (i.e. scalar value of 1).
%
%  EXAMPLE
%
%    fm = FiringMap(positions,spikes);      % firing map for a place cell
%    figure; PlotColorMap(fm.rate,fm.time);  % plot, dimming with occupancy map
%
%  SEE
%
%    See also FiringMap, PhaseMap, MTSpectrogram, PlotShortTimeCCG.

% Copyright (C) 2004-2012 by Michaël Zugaro
% (C) 2025 by Pietro Bozzo (graphics optimization)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    data (:,:) double
    dimm = ~isnan(data)
    opt.x (:,1) = 1 : size(data,2)
    opt.y (:,1) = 1 : size(data,1)
    opt.threshold (1,1) {mustBeNonnegative} = 0.01
    opt.cutoffs (1,:) = []
    opt.hgamma (1,1) = NaN
    opt.gamma (1,1) {mustBeNonnegative} = 1
    opt.bar string = "off"
    opt.barProp cell = {}
    opt.type string {mustBeMember(opt.type,["linear","circular"])} = "linear"
    opt.map = ''
    opt.ydir string {mustBeMember(opt.ydir,["normal","reverse"])} = "normal"
    opt.piecewise {mustBeGeneralLogical} = true
    opt.ax matlab.graphics.axis.Axes = gca
end

% Validate parameters
if ~isempty(opt.cutoffs) && (numel(opt.cutoffs) ~= 2 || opt.cutoffs(1) > opt.cutoffs(2))
    error('Incorrect value for property ''cutoffs'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
end
if ~isText(opt.map,'scalar',true) && ~(isnumeric(opt.map) && size(opt.map,2) == 3)
    error('Incorrect value for property ''map'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
end

% Special defaults
if isnan(opt.hgamma)
    opt.hgamma = 1 / opt.gamma;
elseif opt.hgamma < 0
    error('Incorrect value for property ''hgamma'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
end
default_cutoffs = [min(data,[],'all'),max(data,[],'all')];
if isempty(opt.cutoffs)
    opt.cutoffs = [NaN,NaN];
end
opt.cutoffs(isnan(opt.cutoffs)) = default_cutoffs(isnan(opt.cutoffs));
m = opt.cutoffs(1);
M = opt.cutoffs(2);
if m == M
    M = m + 1;
end
if isscalar(dimm)
    dimm = dimm * ones(size(data));
end
if strcmp(opt.map,'')
    opt.map = Bright(100,'hgamma',opt.hgamma,'type',opt.type);
end
fig = gcf;

% Plot data
data = squeeze(data);
dimm = squeeze(dimm);
p = imagesc(opt.ax,opt.x,opt.y,data,[m M]);
set(opt.ax,'color',[0 0 0]);
if any(dimm(:)~=1)
    alpha(p,1./(1+opt.threshold./(dimm+eps)));
end

% Set X and Y axes
set(opt.ax,'ydir',opt.ydir,'tickdir','out','box','off');
if opt.piecewise && length(opt.x) > 1
    PiecewiseLinearAxis(opt.x,'ax',opt.ax);
end
if opt.piecewise && length(opt.y) > 1
    PiecewiseLinearAxis(opt.y,'y','ax',opt.ax);
end

% Color map and bar
colormap(opt.ax,opt.map);
if ~strcmpi(opt.bar,'off')
    b = colorbar(opt.ax,'vert','TickDirection','out','FontSize',opt.ax.FontSize*.8,'Color',[0,0,0],'Box','off','LineWidth',opt.ax.LineWidth,opt.barProp{:});
    if ~strcmpi(opt.bar,'on')
        b.Label.String = opt.bar;
    end
    set(fig,'currentaxes',opt.ax);
end