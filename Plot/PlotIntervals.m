function PlotIntervals(intervals,varargin)

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
%                   
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

% Default values
color = [0.9 0.9 0.9];
alpha = 0.5;
style = 'rectangles';
direction = 'v';
yLim = ylim;
legendValue = 'on';
bottom = 'on';

if nargin < 1
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
end
if size(intervals,2) ~= 2
  error('Incorrect list of intervals (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
end

% Backwards compatibility: previous versions used the syntax PlotIntervals(intervals,style,direction)
parsed = false;
if (nargin == 2 || nargin == 3) && isastring(lower(varargin{1}),'bars','rectangles','axis')
    style = lower(varargin{1});
    parsed = true;
end
if nargin == 3 && isastring(lower(varargin{2}),'h','v')
    direction = lower(varargin{2});
    parsed = true;
end

% Parse parameter list
if ~parsed
	for i = 1:2:length(varargin)
		if ~ischar(varargin{i})
			error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).']);
		end
		switch lower(varargin{i})
			case 'style'
				style = lower(varargin{i+1});
				if ~isastring(style,'bars','rectangles','axis')
					error('Incorrect value for property ''style'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
				end
			case 'direction'
				direction = lower(varargin{i+1});
				if ~isastring(direction,'h','v')
					error('Incorrect value for property ''direction'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
				end
			case 'color'
				color = lower(varargin{i+1});
				if ~isastring(color,'r','g','b','c','m','y','k','w') && ~isdvector(color,'#3','>=0','<=1')
					error('Incorrect value for property ''direction'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
				end
			case 'alpha'
                alpha = varargin{i+1};
                if ~isdscalar(alpha,'>=0','<=1')
                    error('Incorrect value for property ''alpha'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
                end
            case 'ylim'
                yLim = varargin{i+1};
                if ~isdvector(yLim,'<')
                    error('Incorrect value for property ''yLim'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
                end
            case 'legend'
                legendValue = lower(varargin{i+1});
            case 'bottom'
                bottom = lower(varargin{i+1});
                if islscalar(bottom), if bottom, bottom = 'on'; else, bottom = 'off'; end; end % accept boolean input for backwards compatibility
                if ~isastring(bottom,'on','off')
                    error('Incorrect value for property ''bottom'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
                end
            otherwise
				error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).']);
		end
	end
end

% Validate legend value
if isastring(style,'bars') && ~isastring(legendValue,'on','off')
    error('Incorrect value for property ''legend'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
end

hold on;
xLim = xlim;
if strcmp(style,'bars')
    for i = 1:size(intervals,1)
		if strcmp(direction,'v')
			plot([intervals(i,1) intervals(i,1)],yLim,'Color',[0 0.75 0],'HandleVisibility',legendValue);
			plot([intervals(i,2) intervals(i,2)],yLim,'Color',[0.9 0 0],'HandleVisibility',legendValue);
		else
			plot(xLim,[intervals(i,1) intervals(i,1)],'Color',[0 0.75 0],'HandleVisibility',legendValue);
			plot(xLim,[intervals(i,2) intervals(i,2)],'Color',[0.9 0 0],'HandleVisibility',legendValue);
		end
    end

elseif strcmp(style,'rectangles')
    if strcmp(direction,'v')
        x_coord = intervals(:,[1,1,2,2]).';
        y_coord = yLim(ones(size(intervals,1),1),[1,2,2,1]).';
    else
        x_coord = xLim(ones(size(intervals,1),1),[1,1,2,2]).';
        y_coord = intervals(:,[1,2,2,1]).';
    end
    h = patch(x_coord,y_coord,color,'FaceAlpha',alpha,'LineStyle','none');
    if strcmp(legendValue,'off')
        RemoveFromLegend(h)
    elseif ~strcmp(legendValue,'on')
        h.DisplayName = legendValue;
    end
    % if requested, lower rectangles to bottom of visual stack
    if strcmp(bottom,'on')
        uistack(h,'bottom')
    end

else
    % get current axis location
    ax = gca; fig = gcf;
    pos = ax.Position;
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
            bar_ax = axes(fig,Position=[pos(1),pos(2)-0.01,pos(3),0.01],Color='none',XColor='none',YColor='none',XLim=ax.XLim,Tag='bar_ax_x'); hold on
            linkaxes([ax,bar_ax],'x') % link horizontal zoom of invisible axis to primary axis
        end
        plot(bar_ax,intervals(:),ones(size(intervals(:)))*0.1,Color=color,LineWidth=4)
    else
        for axi = all_axes.'
            if axi.Tag == "bar_ax_y" && all(axi.Position([2,4]) == pos([2,4]))
                bar_ax = axi;
            end
        end
        if isempty(bar_ax)
            bar_ax = axes(fig,Position=[pos(1)-0.01,pos(2),0.01,pos(4)],Color='none',XColor='none',YColor='none',YLim=ax.YLim,Tag='bar_ax_y'); hold on
            linkaxes([ax,bar_ax],'y')
        end
        plot(bar_ax,ones(size(intervals(:)))*0.1,intervals(:),Color=color,LineWidth=4)
    end
    % add bars to legend
    if ~strcmp(legendValue,'off') && ~strcmp(legendValue,'on')
        plot(ax,nan,nan,Color=color,LineWidth=4,DisplayName=legendValue)
    end
    % reset current axes, also raising it to top of visual stack
    axes(ax)
end