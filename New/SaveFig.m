function SaveFig(filename)

%SaveFig - Save current figure in fig, svg, png, and pdf formats
%
% Copyright (C) 2018 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


set(gcf,'renderer','painters'); % otherwise svg can be saved as png if figure too big (uses OpenGL instead of painters)
[folder,name] = fileparts(filename);
filename = fullfile(folder,name); % gets rid of extension, if any
try saveas(gcf,[filename '.fig']); end
saveas(gcf,[filename '.png']);
saveas(gcf,[filename '.svg']);
% saveas(gcf,[filename '.pdf']); % These are almost always cut in a weird way; perform inkscape command to rpoduce pdf instead
try evalc(['!inkscape ' filename '.svg  --export-filename=' filename '.pdf']); end

