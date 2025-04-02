function V = LinearVelocity(X,smooth, varargin)

%LinearVelocity - Compute instantaneous linear velocity.
%
% Compute linear velocity for a time-varying vector X.
%
%  USAGE
%
%    V = LinearVelocity(X,smooth, <options>)
%
%    X              position <a href="matlab:help samples">samples</a>
%    smooth         optional standard deviation for Gaussian kernel used for
%                   differentiating, measured in number of samples
%                   (default = no smoothing)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'split'       if the data is split in time by an interval greater than maxGap
%                   (e.g. two subsessions), compute the linear velocity separately ('on')
%                   or handle them as contiguous points ('off', default).
%     'maxGap'      if a time gap between data points exceeds this threshold,
%                   the data is split in chunks and processed separately (in s, default 60)
%    =========================================================================
%  
%
%    See also AngularVelocity.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help LinearVelocity">LinearVelocity</a>'' for details).');
end
if nargin >= 2,
	if ~isdscalar(smooth,'>=0'),
		error('Incorrect smoothing stdev (type ''help <a href="matlab:help LinearVelocity">LinearVelocity</a>'' for details).');
	end
else
	smooth = 0;
end


% Default parameters
maxGap = 60;
split = 'off';
splitIndx = [1 size(X,1)];


% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+firstIndex) ' is not a property (type ''help LinearVelocity'' for details).']);
	end
	switch(lower(varargin{i})),
	    case 'maxgap',
		maxGap = varargin{i+1};
		if ~isdscalar(maxGap,'>0'),
			error('Incorrect value for property ''maxGap'' (type ''help LinearVelocity'' for details).');
        end
	    case 'split',
		split = varargin{i+1};
		if ~isastring(split,'on','off'),
			error('Incorrect value for property ''split'' (type ''help LinearVelocity'' for details).');
		end
	otherwise,
		error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help LinearVelocity'' for details).']);
	end
end


% Look for gaps in time and process blocks of data separately
% Namely, to split different subsessions and avoid the contamination of speed
% between subsessions at the borders (due to smoothing, the bordering data
% points influence eachother)
if strcmp(split,'on')
    gaps = diff(X(:,1));
    splitIndx = find(gaps>maxGap);
    
    if ~isempty(splitIndx),
        splitIndx = reshape(splitIndx,[],1);
        splitIndx = [[1; splitIndx+1] [splitIndx; size(X,1)]];
    else
        splitIndx = [1 size(X,1)];
    end
end


V = [];
for i = 1:size(splitIndx,1)
    DX = Diff(X(splitIndx(i,1):splitIndx(i,2),:),'smooth',smooth);
    Y = DX(:,2:3).*DX(:,2:3);
    N = sqrt(Y(:,1)+Y(:,2));
    V = [V; [X(splitIndx(i,1):splitIndx(i,2),1) N]];
end