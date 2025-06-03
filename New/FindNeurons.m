function [selectedNeurons,logicalLoc,missingNeuorns] = FindNeurons(toLookFor,whereToFind)

% FindNeurons - find [group unit] pair conserved between two matrices.
%
% Often helpful when needing to select a subgroup of neurons and we have a
% list of [group unit] IDs.
%
% toLookFor         2D matrix ([group unit]) with neurons to be found
% whereToFind       2D matrix ([group unit]) with all the units, this is
%                   where the neurons of toLookFor will be searched
%
% OUTPUT -----------------------------------------------------------------
%
% selectedNeurons   same format as the input but with only the neurons
%                   present in both matrices
% logicalLoc        logical vector of length(whereToFind) with True only
%                   where the target neurons are present
% missingNeuorns    logical vector of size(toLookFor) true where neurons
%                   where not found in both inputs
%
% (c) 2025 Federica LARENO FACCINI
%


[logicalLoc, b] = ismember(whereToFind,toLookFor, 'rows');
selectedNeurons = whereToFind(logicalLoc,:);

test = unique(b); test(test==0) = [];
missingNeuorns = (~ismember(1:length(toLookFor),test))';
end

