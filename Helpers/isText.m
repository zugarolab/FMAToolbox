function out = isText(value,opt)

%isText - Test if value is character vector or string
%
%  USAGE
%
%    test = isText(value)
%
%    value          item to test
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'scalar'      if 'on' (default = 'off'), reject non-scalar input,
%                   accepting values such as 'word', '', "word", ""
%
%  NOTE
%
%    This function doesn't reject empty character vectors or string, such
%    as '', "", string.empty
%
%  SEE
%
%    See also isastring, isdmatrix, isdvector, isdscalar, isimatrix, isivector, isiscalar

% Copyright (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    value
    opt.scalar {mustBeLogicalScalar} = false
end

out = ischar(value) || isstring(value);

opt.scalar = GeneralLogical(opt.scalar);
if out && opt.scalar
    if ischar(value)
        out = out && sum(size(value)>1)<2;
    else
        out = out && isscalar(value);
    end
end