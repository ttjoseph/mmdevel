function [Y, I, J] = shake(X,dim)
% SHAKE - Randomize a matrix along a specific dimension
%   Y = SHAKE(X) randomizes the order of the elements in each column of the
%   2D matrix. For N-D matrices it randomizes along the first non-singleton
%   dimension.
%
%   SHAKE(X,DIM) randomizes along the dimension DIM.
%
%   [Y, I, J] = SHAKE(X) returns indices so that Y = X(I) and X = Y(J).
%
%   Example:
%     A = [1 2 3 ; 4 5 6 ; 7 8 9 ; 10 11 12] ; % see <SLM> on the FEX ...
%     Dim = 2 ;
%     B = shake(A,Dim)  % ->, e.g.
%      %  3     2     1
%      %  6     4     5
%      %  7     8     9
%      % 11    10    12%   
%     C = sort(B,Dim) % -> equals A!
%
%     The function of SHAKE can be thought of as holding a matrix and shake
%     in a particular direction (dimension), so that elements are getting
%     shuffled within that direction only.
%
%   See also RAND, SORT, RANDPERM
%   and RANDSWAP on the File Exchange

% for Matlab R13
% version 4.1 (may 2008)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% Created: dec 2005
% Revisions
% 1.1 : changed the meaning of the DIM. Now DIM==1 works along the rows, preserving
% columns, like in <sort>.
% 2.0 (aug 2006) : randomize along any dimension
% 2.1 (aug 2006) : output indices argument
% 3.0 (oct 2006) : new & easier algorithm
% 4.0 (dec 2006) : fixed major error in 3.0
% 4.1 (may 2008) : fixed error for scalar input

error(nargchk(1,2,nargin)) ;

if nargin==1, 
    dim = min(find(size(X)>1)) ;
elseif (numel(dim) ~= 1) || (fix(dim) ~= dim) || (dim < 1),
    error('Shake:DimensionError','Dimension argument must be a positive integer scalar.') ;
end

% we are shaking the indices
I = reshape(1:numel(X),size(X)) ;

if numel(X) < 2 || dim > ndims(X) || size(X,dim) < 2,    
    % in some cases, do nothing
else
    % put the dimension of interest first
    [I,ndim] = shiftdim(I,dim-1) ;
    sz = size(I) ;
    % reshape it into a 2D matrix
    % we'll randomize along rows
    I = reshape(I,sz(1),[]) ;
    [ri,ri] = sort(rand(size(I)),1) ;  % get new row indices
    ci = repmat([1:size(I,2)],size(I,1),1) ; % but keep old column indices
    I = I(sub2ind(size(I),ri,ci)) ; % retrieve values    
    % restore the size and dimensions
    I = shiftdim(reshape(I,sz),ndim) ;    
end

% re-index
Y = X(I) ; 

if nargout==3,
    J = zeros(size(X)) ;
    J(I) = 1:numel(J) ;
end


