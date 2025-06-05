function [matched,logicalReference,missingElements] = MatchMatrices(reference,query)

% MatchMatrices - match matrices and find conserved rows between them.
%
% USAGE
%
%   [matched,logicalReference,missingElements] = MatchMatrices(reference,query)
%
%   reference           2D matrix where the elements of query will searched
%   query               2D matrix with elements to look for
%
%
% OUTPUT -----------------------------------------------------------------
%
%   matched             same as the query input but with only the rows that
%                       were matched with rows in the reference matrix
%   logicalReference	logical vector of length(reference) true only where
%                       the rows where matched to the query 
%   missingElements     logical vector of size(query) true where elements
%                       where NOT matched to any row in the reference
%
% Note
% 
% This function matches two matrices and finds rows (combination of all columns!)
% that appear identical between the two. For this reason the two inputs
% should have the same number of columns.
% It looks up a query matrix in a reference matrix and returns a new matrix (matched)
% with only the rows present in both elements.
% It also returns a vector (logicalReference) with size(reference), true
% where rows of the query where found in reference, else false.
% The third output is a logical vector (missingElements) of size(query) true where rows where
% NOT matched with any row in the reference.
%
% Often helpful when needing to select a subgroup of neurons and we have a
% list of [group unit] IDs.
%
% (c) 2025 Federica LARENO FACCINI
%


[logicalReference, b] = ismember(reference,query, 'rows');
matched = reference(logicalReference,:);

test = unique(b); test(test==0) = [];
missingElements = (ismember(1:length(query),test))';
end

