function value = GeneralLogical(a)

%GeneralLogical - Cast input to logical array
%
%  USAGE
%
%    value = GeneralLogical(a)
%
%    a               generalized logical: can take values in:
%                    [true, false, 1, 0, "on", "off"]
%
%  OUTPUT
%
%    value           converted input

% Copyright (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

value = a;

if ischar(a)
    a = string(a);
end

if isstring(a)
    if ~all(a(a~="off") == "on")
        error('mustBeLogical:notLogical','Value must be logical.')
    end
    value = a == "on";
elseif isnumeric(a)
    if ~all(a(a~=0) == 1)
        error('mustBeLogical:notLogical','Value must be logical.')
    end
    value = a ~= 0;
elseif ~islogical(a)
    error('mustBeLogical:notLogical','Value must be logical.')
end