function indices = findmin(data,k)
%findmin - find local maximum with a  window length of k, if k<2, it
%calculates global maximum 
%
% (C) 2016-2025 by Ralitsa Todorova 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin<2
    [~,indices] = min(data);
else
    [~,order] = sort(data,'ascend');
    indices = order(1:k);
end