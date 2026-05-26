function varargout = semplotGrouped(x,y,color,smooth,solid)

%semplot - plot mean (line) +/- s.e.m. (shaded area) of a matrix "y"
% semplot(x,y,color,smooth)
%
% Copyright (C) 2026 by Ralitsa Todorova
%
% SEE ALSO shadedErrorBar 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ~exist('solid','var')
    solid = false;
end

if ~exist('color','var')
    color = [0 0 0];
end

if ~exist('smooth','var')
    smooth = 0;
end

if numel(y)~=numel(x)
    error('Y should have one value for each element in X');
end
x = x(:); y = y(:);
bad = isnan(x) | isnan(y); x(bad) = []; y(bad) = [];
[u,~,index] = unique(x);

n = Accumulate(index,1);
m = Accumulate(index,y)./n;
% get sem
d = (y-m(index)).^2;
SEM = (sqrt(Accumulate(index,d)./n))./(sqrt(n-1));

xx = [u(:);flipud(u(:))];
yy = [Smooth(m-SEM,smooth,'type','c'); Smooth(flipud(m+SEM),smooth,'type','c')];
y = Smooth(m,smooth,'type','c');
handles = fill(xx,yy,color);

if solid
    set(handles,'FaceColor',mean([color;1 1 1]),'edgeAlpha',0);
else % transparent
    set(handles,'FaceAlpha',0.5,'edgeAlpha',0);
end

hold on;
plot(u,y,'color',color,'linewidth',2);

if nargout>0, varargout{1} = handles; end