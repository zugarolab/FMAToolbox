function PlotColorMap(data,dimm,varargin)

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

% Default values
cutoffs = [];
hgamma = 1;
gamma = 1;
hg = 0;
threshold = 0.01;
drawBar = false;
type = 'linear';
[y,x] = size(data);
x = 1:x; y = 1:y;
map = '';
ydir = 'normal';
piecewise = true;
ax = [];

% Check parameters
if nargin < 1
	error('Incorrect number of parameters (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
end
if nargin == 1
	dimm = 1;
end
if isa(dimm,'char')
	varargin = [dimm,varargin];
	dimm = 1;
end

% Parse parameter list
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).']);
    end
    switch lower(varargin{i})

        case 'threshold'
            threshold = varargin{i+1};
            if ~isdscalar(threshold,'>=0')
                error('Incorrect value for property ''threshold'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
            end
  		  case 'x'
    			  x = varargin{i+1};
      			if ~isdvector(x)
      				  error('Incorrect value for property ''x'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
      			end
  		  case 'y'
    	  		y = varargin{i+1};
    		  	if ~isdvector(y)
    			  	  error('Incorrect value for property ''y'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
    	  		end
  		  case 'cutoffs'
    			  cutoffs = varargin{i+1};
      			if ~isempty(cutoffs) && (~isvector(cutoffs) || numel(cutoffs) ~= 2 || cutoffs(1) > cutoffs(2))
      				  error('Incorrect value for property ''cutoffs'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
      			end
  		  case 'hgamma'
    			  hg = 1;
      			hgamma = varargin{i+1};
      			if ~isdscalar(hgamma,'>=0')
        				error('Incorrect value for property ''hgamma'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
    	  		end
  		  case 'gamma'
    			  gamma = varargin{i+1};
            if ~isdscalar(gamma,'>=0')
      			  	error('Incorrect value for deprecated property ''gamma'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
            end
  		  case 'bar'
      			barLabel = varargin{i+1};
            drawBar = ~strcmpi(varargin{i+1},'off');
      			if ~isastring(barLabel)
    	  		  	error('Incorrect value for property ''bar'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
    		  	end
  	  	case 'type'
  			    type = lower(varargin{i+1});
            if ~isastring(type,'linear','circular')
                error('Incorrect value for property ''type'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
            end
        case 'map'
            map = lower(varargin{i+1});
            if ~isText(map,'scalar',true) && ~(isnumeric(map) && size(map,2) == 3)
                error('Incorrect value for property ''map'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
            end
    		case 'ydir'
    	  		ydir = lower(varargin{i+1});
            if ~isastring(ydir,'normal','reverse')
    	  		  	error('Incorrect value for property ''ydir'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
            end
        case 'piecewise'
    			  piecewise = strcmpi(varargin{i+1},'on');
            if ~isastring(lower(varargin{i+1}),'on','off')
      		  		error('Incorrect value for property ''piecewise'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
            end
        case 'ax'
            ax = varargin{i+1};
            if ~isa(ax,'matlab.graphics.axis.Axes')
                error('Incorrect value for property ''ax'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
            end
        otherwise
    		  	error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).']);
  	end
end

data = double(data);

% Special defaults
x = x(:);
y = y(:);
if hg == 0
    hgamma = 1/gamma;
end
default_cutoffs = [min(data,[],'all'),max(data,[],'all')];
if isempty(cutoffs)
    cutoffs = [NaN,NaN];
end
cutoffs(isnan(cutoffs)) = default_cutoffs(isnan(cutoffs));
m = cutoffs(1);
M = cutoffs(2);
if m == M
    M = m + 1;
end
if isscalar(dimm)
    dimm = dimm*ones(size(data));
end
if strcmp(map,'')
    map = Bright(100,'hgamma',hgamma,'type',type);
end
if isempty(ax)
    ax = gca;
end
fig = gcf;

% Plot data
data = squeeze(data);
dimm = squeeze(dimm);
p = imagesc(ax,x,y,data,[m M]);
set(ax,'color',[0 0 0]);
if any(dimm(:)~=1)
    alpha(p,1./(1+threshold./(dimm+eps)));
end

% Set X and Y axes
set(ax,'ydir',ydir,'tickdir','out','box','off');
if piecewise && ~isempty(x) && length(x) ~= 1
    PiecewiseLinearAxis(x,'ax',ax);
end
if piecewise && ~isempty(y) && length(y) ~= 1
    PiecewiseLinearAxis(y,'y','ax',ax);
end

% Color map and bar
colormap(ax,map);
if drawBar
    b = colorbar(ax,'vert','TickDirection','out','FontSize',12,'Color',[0,0,0],'Box','off','LineWidth',1.7);
        if ~strcmpi(barLabel,'on'), b.Label.String = barLabel; end
    set(fig,'currentaxes',ax);
end