function interpolated = fillnans(matrix,coefficients)

% coeffifcients are vertical and horizontal kernel size (like in Smooth):
% how much should each dimension weigh in the interpolation?

% This function is an usage of John D'Errico's method 2 in
% INPAINT_NANS (dowloaded from FileExchange on 01/02/2016)
% e-mail address: woodchips@rochester.rr.com
% which can do exactly the same interpolation process for
% coefficients  =  [1 1];
% INPAINT_NANS: in-paints over nans in an array
% usage: B = INPAINT_NANS(A)          % default method
% usage: B = INPAINT_NANS(A,method)   % specify method used
%
% Solves approximation to one of several pdes to
% interpolate and extrapolate holes in an array
%
% arguments (input):
%   A - nxm array with some NaNs to be filled in
%
%         This method uses a simple plate metaphor.
%         del^2 is used.
%         Extrapolation behavior is linear.
%
%         Note: This method has problems in 1-d, so this
%         method is disabled for vector inputs.


% I always need to know which elements are NaN,
% and what size the array is for any method
[nRows,nCols] = size(matrix);
matrix = matrix(:);
nm = nRows*nCols;
k = isnan(matrix(:));

if nargin<2,
    coefficients = [1 1];
end

% list the nodes which are known, and which will
% be interpolated
nan_list = find(k);
known_list = find(~k);

% convert NaN indices to (r,c) form
% nan_list =  = find(k) are the unrolled (linear) indices
% (row,column) form
[nr,nc] = ind2sub([nRows,nCols],nan_list);

% both forms of index in one array:
% column 1  =  =  unrolled index
% column 2  =  =  row index
% column 3  =  =  column index
nan_list = [nan_list,nr,nc];

% Direct solve for del^2 BVP across holes

% generate sparse array with second partials on row
% variable for each nan element, only for those nodes
% which have a row index > 1 or < n


% a 2-d case
L  =  find((nan_list(:,2) > 1) & (nan_list(:,2) < nRows));
nl = length(L);
if nl>0
    fda = sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
        repmat([1 -2 1]*coefficients(1),nl,1),nRows*nCols,nRows*nCols);
else
    fda = spalloc(nRows*nCols,n*nCols,size(nan_list,1)*5);
end

% 2nd partials on column index
L  =  find((nan_list(:,3) > 1) & (nan_list(:,3) < nCols));
nl = length(L);
if nl>0
    fda = fda+sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-nRows 0 nRows],nl,1), ...
        repmat([1 -2 1]*coefficients(2),nl,1),nRows*nCols,nRows*nCols);
end

% fix boundary conditions at extreme corners
% of the array in case there were nans there
if ismember(1,nan_list(:,1)),fda(1,[1 2 nRows+1]) = [-2 1 1];end
if ismember(nRows,nan_list(:,1)),fda(nRows,[nRows, nRows-1,nRows+nRows]) = [-2 1 1];end
if ismember(nm-nRows+1,nan_list(:,1)),fda(nm-nRows+1,[nm-nRows+1,nm-nRows+2,nm-nRows]) = [-2 1 1];end
if ismember(nm,nan_list(:,1)),fda(nm,[nm,nm-1,nm-nRows]) = [-2 1 1];end

% eliminate knowns
rhs = -fda(:,known_list)*matrix(known_list);

% and solve...
interpolated = matrix;
k = nan_list(:,1);
interpolated(k) = fda(k,k)\rhs(k);


% all done, make sure that B is the same shape as
% A was when we came in.
interpolated = reshape(interpolated,nRows,nCols);
end
