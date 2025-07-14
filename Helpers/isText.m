function out = isText(value)

%isText - Test if value is character vector or string
%
%  USAGE
%
%    test = isText(value)
%
%    value          item to test
%
%  NOTE
%
%    This function doesn't reject empty character vectors or string, such
%    as '', "", strings().empty
%
%  SEE ALSO
%
%    See also isastring, isdmatrix, isdvector, isdscalar, isimatrix, isivector, isiscalar.
%

% Copyright (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

out = ischar(value) || isstring(value);