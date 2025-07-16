function mustBeLogical(A)
%MUSTBELOGICAL Validate that value is logical
%   MUSTBELOGICAL(A) throws an error if A contains values different from 0 and 1.

if (~isnumeric(A)&&~islogical(A)) || ~all(A(A~=0)==1)
    error('mustBeLogical:notLogical','Value must be logical.')
end