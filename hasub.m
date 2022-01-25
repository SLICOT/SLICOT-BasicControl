function [US,UU] = hasub(S,T,U,V,select)
% HASUB  Selected stable/unstable invariant subspace of a
%        Hamiltonian matrix.
%    [US] = HASUB(S,T,U,V,SELECT) produces an orthonormal basis
%    US for the invariant subspace belonging to a selected set of
%    eigenvalues with negative real part. S, T, U and V must be
%    the matrices of a symplectic URV/periodic Schur decomposition
%    of a Hamiltonian matrix H as computed by HAURVPS. If -w(i)^2
%    is an eigenvalue corresponding to the i-th diagonal element
%    of S*T, then w(i) is an eigenvalue of H. This eigenvalue is
%    selected if the i-th element of the logical or real n-vector
%    SELECT is not zero.
%
%    [US,UU] = HASUB(S,T,U,V,SELECT) also produces an orthonormal
%    basis UU for the invariant subspace belonging to a selected
%    set of eigenvalues with positive real part.
%
%    See also HASTAB, HAURVPS.

%    This is a Matlab gateway routine to the HAPACK routine
%    MB03ZD.f.
%
%    RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%    Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%    Based on HAPACK package, http://www.tu-chemnitz.de/mathematik/hapack/
%
%    Revision
%    V. Sima, Oct. 2008, Mar. 2009, Nov. 2012.

if nargin < 5,
   error('Five input arguments required.');
end
%
% Check input arguments.
%
n = size(S,1);
if (size(S,2) ~= n),
   error('S is not a square matrix');
end
if (size(T,1) ~= n) && (size(T,2) ~= n),
   error('T is not an n-by-n matrix');
end
if (size(U,1) ~= 2*n) && (size(U,2) ~= 2*n),
   error('U is not a 2*n-by-2*n matrix');
end
if (size(V,1) ~= 2*n) && (size(V,2) ~= 2*n),
   error('V is not a 2*n-by-2*n matrix');
end
if min(size(select)) > 1 || length(select) ~= n,
   error('SELECT must be a vector of length n.');
end
selct = zeros(n,1);
if islogical(select),
   for ii = 1:n
      if select(ii),  selct(ii) = 1;  else  selct(ii) = 0;  end
   end
elseif isreal(select)
   selct = select;
else
   error('SELECT must be a logical or real vector.');
end

if nargout > 1,
   [US,UU] = hapack_haeig(7,S,T,...
                          U(1:n,1:n), U(1:n,n+1:2*n),...
                          V(1:n,1:n), V(1:n,n+1:2*n),...
                          selct);
else
        US = hapack_haeig(7,S,T,...
                          U(1:n,1:n), U(1:n,n+1:2*n),...
                          V(1:n,1:n), V(1:n,n+1:2*n),...
                          selct);
end
