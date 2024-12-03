function [c,slope,rSquared,SSres,startstop] = SpearmanByColumn(A,B)

% transforms data into rank data, and then calls CorrByColumn
% each column in A is compared to each column in B
% e.g.
% c = SpearmanByColumn(sequenceList',repmat([1:size(sequenceList,2)],size(sequenceList,1),1)')

A = tiedrank(A);
B = tiedrank(B);

vargout = cell(1,nargout);

[vargout{:}] = CorrByColumn(A,B);

try
    c=vargout{1};
    slope=vargout{2};
    rSquared=vargout{3};
    SSres=vargout{4};
    startstop=vargout{5};
end
