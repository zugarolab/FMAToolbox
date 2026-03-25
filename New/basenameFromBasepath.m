function basename = basenameFromBasepath(basepath)
% basenameFromBasepath-gets basename from a basepath input to a function.  Uses name of the last
% directory on the path to give basename.  Takes care of the problem of the
% fact that if a fileseparator is on the end of the basepath, you have to
% take it off first.
%
% INPUTS
% basepath  is assumed to be a directory containing all data for a given 
%           recording session.  It is assumed organization is 
%           /whateverdirectory/basename/basename.xml and basepath is 
%           equivalent to whateverdirectory/basename/.  
%
% OUTPUTS
% basename  as specified above
%
% Brendon Watson, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if strcmp(basepath(end),filesep);
    basepath = basepath(1:end-1);
end
% if there are dots in the filename take that last folder
if contains(basepath,'.')
    parts = strsplit(basepath,filesep);
    basename = parts{end};
    return
end
[~,basename] = fileparts(basepath);
end
