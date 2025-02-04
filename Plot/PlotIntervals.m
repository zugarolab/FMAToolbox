function PlotIntervals(intervals,varargin)

%PlotIntervals - Plot vertical bars or rectangles to show interval limits.
%
% Given a list of intervals [start stop], draw a green vertical bar at
% the beginning of each interval, and a red vertical bar at the end of
% each interval or a grey rectangle representing the interval.
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
%     'style'       'bars' for colored vertical bars, or 'rectangles' for
%                   background shaded rectangles (default)
%     'direction'   'h' for horizontal, or 'v' for vertical (default)
%     'color'       rectangle color ('rectangles' mode, default = grey)
%     'alpha'       rectangle transparency ('rectangles' mode, default = 1)
%     'ylim'        desired y-coordinates of the plotted areas
%     'legend'      if 'off', plotted elements won't appear in legend
%                   (default = 'on')
%     'bottom'      if 'on' (default = 'on'), lower the newly plotted
%                   elements to the bottom of visual stack. Note that
%                   setting this to 'off' is faster, so consider this (and
%                   PlotIntervals before other plot funcitons) for large
%                   numbers of intervals.
%    =========================================================================
%

% Copyright (C) 2008-2013 by Gabrielle Girardeau & Michaël Zugaro, 
% (C) 2023-2025 by Ralitsa Todorova and Pietro Bozzo (graphics optimization) 
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
legend = 'on';
bottom = 'on';

if nargin < 1
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
end
if size(intervals,2) ~= 2
  error('Incorrect list of intervals (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
end

% Backward compatibility: previous versions used the syntax PlotIntervals(intervals,style,direction)
parsed = false;
if (nargin == 2 || nargin == 3) && isastring(lower(varargin{1}),'rectangles','bars')
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
				if ~isastring(style,'bars','rectangles')
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
                legend = lower(varargin{i+1});
                if ~isastring(legend,'on','off')
					error('Incorrect value for property ''legend'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
                end
            case 'bottom'
                bottom = lower(varargin{i+1});
                if islscalar(bottom), if bottom, bottom = 'on'; else, bottom = 'off'; end; end % accept boolean input for back compatibility
                if ~isastring(bottom,'on','off')
                    error('Incorrect value for property ''bottom'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
                end
            otherwise
				error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).']);
		end
	end
end

hold on;
xLim = xlim;
if strcmp(style,'bars')
	for i = 1:size(intervals,1)
		if strcmp(direction,'v')
			plot([intervals(i,1) intervals(i,1)],yLim,'Color',[0 0.75 0],'HandleVisibility',legend);
			plot([intervals(i,2) intervals(i,2)],yLim,'Color',[0.9 0 0],'HandleVisibility',legend);
		else
			plot(xLim,[intervals(i,1) intervals(i,1)],'Color',[0 0.75 0],'HandleVisibility',legend);
			plot(xLim,[intervals(i,2) intervals(i,2)],'Color',[0.9 0 0],'HandleVisibility',legend);
		end
	end
else
    nPreexistingObjects = numel(gca().Children); % to use later to restore lines order in plot
    for i = 1:size(intervals,1)
        if strcmp(direction,'v')
			dx = intervals(i,2)-intervals(i,1);
			dy = yLim(2)-yLim(1);
            patch(intervals(i,1)+[0,0,dx,dx],yLim(1)+[0,dy,dy,0],color,'FaceAlpha',alpha,'LineStyle','none','HandleVisibility',legend);
		else
			dx = xLim(2)-xLim(1);
			dy = intervals(i,2)-intervals(i,1);
            patch(xLim(1)+[0,0,dx,dx],intervals(i,1)+[0,dy,dy,0],color,'FaceAlpha',alpha,'LineStyle','none','HandleVisibility',legend);
        end
    end
    nObjects = numel(gca().Children);
    nNewObjects = nObjects-nPreexistingObjects;
    % if rectangles were plotted and if requested, lower them to bottom
    if nNewObjects>0 && strcmp(bottom,'on')
        ax = gca();
        order = [nNewObjects+1:nNewObjects+nPreexistingObjects,1:nNewObjects];
        ax.Children = ax.Children(order);
    end
end