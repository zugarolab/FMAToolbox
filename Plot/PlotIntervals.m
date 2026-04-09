function PlotIntervals(intervals,style,direction,opt)

%PlotIntervals - Plot lines, rectangles or bars next to an axis to show interval limits.
%
% Given a list of [start stop] intervals, draw each of them with either:
% 1. a green (red) line representing its beginning (end)
% 2. a rectangle
% 3. a bar next to the x or y axis
%
%  USAGE
%
%    PlotIntervals(intervals,<options>)
%
%    intervals      list of [start stop] intervals
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'style'       'bars' for lines, 'rectangles' (default), or 'axis'
%     'direction'   'h' for horizontal, or 'v' for vertical (default)
%     'color'       rectangle color ('rectangles' mode, default = grey)
%     'alpha'       rectangle transparency ('rectangles' mode, default = 0.5)
%     'ylim'        desired y-coordinates of the plotted areas ('v' mode)
%     'legend'      if 'off', plotted elements won't appear in legend
%                   (default = 'on'); value is also legend label if
%                   different from 'on'
%     'bottom'      if 'on' (default), lower plotted rectangles to bottom of
%                   of visual stack ('rectangles' mode)
%     'ax'          axis to plot on, default is default plot behavior
%    =========================================================================
%
%  SEE
%
%    See also ConsolidateIntervals, SubtractIntervals, ExcludeIntervals,
%    IntersectIntervals, InIntervals, Restrict, FindInInterval, CountInIntervals.

% Copyright (C) 2008-2013 by Gabrielle Girardeau & Michaël Zugaro, 
% (C) 2023-2025 by Ralitsa Todorova & Pietro Bozzo (graphics optimization)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    intervals (:,2) {mustBeNumeric}
    style (1,1) string = missing % kept for retrocompatibility
    direction (1,1) string = missing % kept for retrocompatibility
    opt.style (1,1) string {mustBeMember(opt.style,["rectangles","bars","axis"])} = "rectangles"
    opt.direction (1,1) string {mustBeMember(opt.direction,["h","v"])} = "v"
    opt.color = [0.9,0.9,0.9]
    opt.alpha (1,1) {mustBeInRange(opt.alpha,0,1)} = 0.5
    opt.ylim (1,2) = NaN
    opt.legend (1,1) string = missing
    opt.bottom {mustBeGeneralLogical} = true
    opt.ax (1,1) matlab.graphics.axis.Axes = gca
end

% validate input
if ismissing(style)
    style = opt.style;
else
    if ~ismember(style,["rectangles","bars","axis"])
        error('Incorrect value for property ''style'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).')
    end
end
if ismissing(direction)
    direction = opt.direction;
else
    if ~ismember(direction,["h","v"])
        error('Incorrect value for property ''direction'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).')
    end
end
color = validatecolor(opt.color);
if diff(opt.ylim) <= 0
    error('Incorrect value for property ''ylim'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).')
end
if style == "bars"
    if ismissing(opt.legend)
        opt.legend = 'on';
    end
    if ~ismember(opt.legend,["on","off"])
        error('Incorrect value for property ''legend'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
    end
end
bottom = GeneralLogical(opt.bottom);
if any(isnan(opt.ylim))
    opt.ylim = ylim(opt.ax);
end

if isempty(intervals)
    return
end

hold on
xLim = xlim(opt.ax);
if strcmp(style,'bars')
    for i = 1:size(intervals,1)
		if strcmp(direction,'v')
			plot(opt.ax,[intervals(i,1) intervals(i,1)],opt.ylim,'Color',[0 0.75 0],'HandleVisibility',opt.legend);
			plot(opt.ax,[intervals(i,2) intervals(i,2)],opt.ylim,'Color',[0.9 0 0],'HandleVisibility',opt.legend);
		else
			plot(opt.ax,xLim,[intervals(i,1) intervals(i,1)],'Color',[0 0.75 0],'HandleVisibility',opt.legend);
			plot(opt.ax,xLim,[intervals(i,2) intervals(i,2)],'Color',[0.9 0 0],'HandleVisibility',opt.legend);
		end
    end

elseif strcmp(style,'rectangles')
    if strcmp(direction,'v')
        x_coord = intervals(:,[1,1,2,2]).';
        y_coord = opt.ylim(ones(size(intervals,1),1),[1,2,2,1]).';
    else
        x_coord = xLim(ones(size(intervals,1),1),[1,1,2,2]).';
        y_coord = intervals(:,[1,2,2,1]).';
    end
    h = patch(opt.ax,x_coord,y_coord,color,'FaceAlpha',opt.alpha,'LineStyle','none');
    if strcmp(opt.legend,'off')
        RemoveFromLegend(h)
    elseif ~ismissing(opt.legend)
        h.DisplayName = opt.legend;
    end
    % if requested, lower rectangles to bottom of visual stack
    if bottom
        uistack(h,'bottom')
    end

else
    % get current axis location
    fig = gcf;
    pos = opt.ax.Position;
    intervals = [intervals,nan(size(intervals,1),1)].';
    all_axes = findall(fig,'type','axes');
    % plot bars in invisible axis next to either x or y axis
    bar_ax = [];
    if strcmp(direction,'v')
        % look for pre-existing bar axis
        for axi = all_axes.'
            if axi.Tag == "bar_ax_x" && all(axi.Position([1,3]) == pos([1,3]))
                bar_ax = axi;
            end
        end
        if isempty(bar_ax)
            % create new bar axis where ticks of primary axis are
            bar_ax = axes(fig,'Position',[pos(1),pos(2)-0.01,pos(3),0.01],'Color','none','XColor','none','YColor','none','XLim',opt.ax.XLim,'Tag','bar_ax_x'); hold on
            linkaxes([opt.ax,bar_ax],'x') % link horizontal zoom of invisible axis to primary axis
        end
        plot(bar_ax,intervals(:),ones(size(intervals(:)))*0.1,'Color',color,'LineWidth',4)
    else
        for axi = all_axes.'
            if axi.Tag == "bar_ax_y" && all(axi.Position([2,4]) == pos([2,4]))
                bar_ax = axi;
            end
        end
        if isempty(bar_ax)
            bar_ax = axes(fig,'Position',[pos(1)-0.01,pos(2),0.01,pos(4)],'Color','none','XColor','none','YColor','none','YLim',opt.ax.YLim,'Tag','bar_ax_y'); hold on
            linkaxes([opt.ax,bar_ax],'y')
        end
        plot(bar_ax,ones(size(intervals(:)))*0.1,intervals(:),'Color',color,'LineWidth',4)
    end
    % add bars to legend
    if ~ismissing(opt.legend) && opt.legend ~= "off"
        plot(opt.ax,nan,nan,'Color',color,'LineWidth',4,'DisplayName',opt.legend)
    end
    % reset current axes, also raising it to top of visual stack
    axes(opt.ax)
end