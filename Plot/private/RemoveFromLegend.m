function RemoveFromLegend(h)

%RemoveFromLegend - Remove graphic objects from axis legend
%
%  USAGE
%
%    RemoveFromLegend(h,fig)
%
%    h              handle(s) to graphic object(s)

% Copyright (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    h (1,:)
end

% set property of graphic objects
for hand = h
    hand.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

% get parent figure
figHandle = ancestor(h(1),'figure');

% if legends exist, refresh them to effectively hide entries
legHandle = findobj(figHandle,'Type','Legend');
for hand = legHandle.'
    hand.Visible = 'off';
    hand.Visible = 'on';
end