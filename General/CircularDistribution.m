function [dist,bin_centers,stats] = CircularDistribution(angles,varargin)

%CircularDistribution - Compute circular distribution and statistics.
%
%  USAGE
%
%    [dist,bin_centers,stats] = CircularDistribution(angles,<options>)
%
%    angles         angles in radians, which will be shifted in [0,2π) 
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'nBins'       number of bins (default = 100)
%     'smooth'      standard deviation of Gaussian kernel (default = 0)
%     'groups'      groups for multiple circular distributions (see below)
%     'normalize'   type of histogram normalization, can be 'probability' (default),
%                   'count', 'pdf' (see <a href="https://www.mathworks.com/help/matlab/ref/double.histcounts.html">histcounts</a> help for a complete list)
%    =========================================================================
%
%  OUTPUT
%
%    dist           circular distribution (one column per group)
%    bin_centers    centers of the angular bins
%    stats.m        mean angle (one per group)
%    stats.mode     distribution mode (one per group)
%    stats.r        mean resultant length (one per group)
%    stats.k        concentration (one per group)
%    stats.p        p-value for Rayleigh test (one per group)
%
%  NOTE
%
%    For multiple circular distributions, groups can be indicated in two different
%    manners:
%
%     - a vector of group IDs (one per angle)
%     - a logical matrix (one line per angle, one column per group), where
%       the element (i,j) is 1 iff angle i belongs to group j
%
%    The vector form is convenient when each angle can only belong to one group.
%    The matrix form is useful when a single angle can belong to multiple groups.
%
%  SEE
%
%    See also PlotCircularDistribution.

% Copyright (C) 2011-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
nBins = 100;
smooth = 0;
normalize = 'probability';

% Check number of parameters
if nargin < 1 || mod(length(varargin),2) ~= 0
  error('Incorrect number of parameters (type ''help <a href="matlab:help CircularDistribution">CircularDistribution</a>'' for details).');
end

% Check parameter size
if ~isdvector(angles)
	error('Incorrect angles (type ''help <a href="matlab:help CircularDistribution">CircularDistribution</a>'' for details).');
end
angles = angles(:);
groups = ones(size(angles));

% Parse parameter list
for i = 1:2:length(varargin)
  if ~ischar(varargin{i})
    error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CircularDistribution">CircularDistribution</a>'' for details).']);
  end
  switch lower(varargin{i})
    case 'nbins'
      nBins = varargin{i+1};
      if ~isiscalar(nBins,'>0')
        error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help CircularDistribution">CircularDistribution</a>'' for details).');
      end
    case 'smooth'
      smooth = varargin{i+1};
      if ~isdscalar(smooth,'>=0')
        error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help CircularDistribution">CircularDistribution</a>'' for details).');
      end
    case 'groups'
      groups = varargin{i+1};
      if ~isdvector(groups,'>0') && ~islmatrix(groups)
        error('Incorrect value for property ''groups'' (type ''help <a href="matlab:help CircularDistribution">CircularDistribution</a>'' for details).');
      end
      if isdvector(groups), groups = groups(:); end
      if length(angles) ~= size(groups,1)
        error('Phases and groups have different numbers of lines (type ''help <a href="matlab:help CircularDistribution">CircularDistribution</a>'' for details).');
      end
    case 'normalize'
      normalize = varargin{i+1};
      if ~isastring(normalize) || size(normalize,1) ~= 1
        error('Incorrect value for property ''normalize'' (type ''help <a href="matlab:help CircularDistribution">CircularDistribution</a>'' for details).');
      end
    otherwise
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CircularDistribution">CircularDistribution</a>'' for details).']);
  end
end

% Remap angles in [0,2π)
angles = mod(angles,2*pi);

% Angle bins
bin_edges = linspace(0,2*pi,nBins+1).';
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

% Groups: transform vector form into matrix form
if isdvector(groups)
	groupIDs = unique(groups);
	nGroups = max(groupIDs);
	g = groups;
	groups = false(length(g),nGroups);
	for i = 1:nGroups
		groups(g==i,i) = true;
	end
end

% Loop through groups
nGroups = size(groups,2);
dist = zeros(size(bin_centers,1),nGroups);
for i = 1 : nGroups
	% Distribution
	this_angles = angles(groups(:,i));
    counts = histcounts(this_angles,bin_edges,'Normalization',normalize);
    counts = Smooth(counts,smooth);
	dist(:,i) = counts;
	% Stats
	[stats.m(i),stats.r(i)] = CircularMean(this_angles);
	stats.k(i) = Concentration(this_angles);
	n = sum(groups(:,i));
	R = stats.r(i)*n;
	stats.p(i) = exp(sqrt(1+4*n+4*(n^2-R^2))-(1+2*n)); % Zar, Biostatistical Analysis, p. 617
	x = find(counts==max(counts));x = x(1);
	stats.mode(i) = bin_centers(x);
end