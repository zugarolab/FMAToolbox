function mustBeGeneralLogical(A)
%MUSTBEGENERALLOGICAL Validate that value is a generalized logical
%   MUSTBEGENERALLOGICAL(A) throws an error if A contains values different from: false, true, 0, 1, "off", "on".

% Copyright (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ischar(A)
    A = string(A);
end

if isstring(A)
    if ~all(A(A~="off") == "on")
        error('mustBeLogical:notLogical','Value must be logical.')
    end
elseif isnumeric(A)
    if ~all(A(A~=0) == 1)
        error('mustBeLogical:notLogical','Value must be logical.')
    end
elseif ~islogical(A)
    error('mustBeLogical:notLogical','Value must be logical.')
end