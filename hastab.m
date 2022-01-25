function [US,UU,e] = hastab(A,QG,meth)
% HASTAB Complete stable/unstable invariant subspace of a
%        Hamiltonian matrix.
%    US = HASTAB(H) is an orthonormal basis of the stable
%    invariant subspace of a 2*n-by-2*n Hamiltonian matrix H,
%    i.e., the invariant subspace belonging to all eigenvalues
%    with negative real part. It is assumed that H has no
%    eigenvalues with zero real part. 
%
%    US = HASTAB(H,METH) performs the computation using the
%    method specified by METH as follows:
%         METH = 1: compute US from a set of N vectors and
%                   apply no scaling;
%         METH = 2: compute US from a set of N vectors and
%                   apply symplectic scaling in order to
%                   reduce the norm of H;
%         METH = 3: compute US from a set of 2*N vectors and
%                   apply no scaling;
%         METH = 4: compute US from a set of 2*N vectors and
%                   apply symplectic scaling.
%    By default, METH = 1. Note that METH = 2 and METH = 4
%    will generally yield non-orthonormal bases. Use
%    [US,r] = QR(US,0) to produce an orthonormal basis. In
%    some cases, METH = 3 and METH = 4 may result in more
%    accurately computed invariant subspaces, for more details
%    see [1].
%
%    [US,E] = HASTAB(H{,METH}) also returns a vector E
%    containing the eigenvalues of H.
%
%    [US,UU,E] = HASTAB(H{,METH}) also computes an orthonormal
%    basis UU for the unstable invariant subspace, i.e., the
%    invariant subspace belonging to the eigenvalues with positive
%    real part.
%
%    [US{,UU,E}] = HASTAB(A,QG{,METH}) requires A and QG contain
%    the matrix H in compressed format.
%
%    [1] P. BENNER, V. MEHRMANN and H. XU: A new method for
%    computing the stable invariant subspace of a real
%    Hamiltonian matrix,  J. Comput. Appl. Math., vol. 86, 1997,
%    pp 17-43.
%
%    See also HAEIG and HASUB.

%    This is a MATLAB gateway routine to the SLICOT routines
%    MB03XD.f and MB03ZD.f.
%
%    RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%    Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%    Based on HAPACK package, http://www.tu-chemnitz.de/mathematik/hapack/
%
%    Revision
%    V. Sima, Oct. 2008, Mar. 2009, Nov. 2012.

if nargin < 1,
   error('At least one input argument required.');
elseif (nargin == 1),
   meth = 1;
   compressed = 0;
elseif ( nargin == 2 && numel(QG) == 1 ),
   meth = QG;
   compressed = 0;
elseif ( nargin <= 3),
   compressed = 1;
   if nargin == 2,
      meth = 1;
   end
else
   error('Too many input arguments');
end

if ( nargout > 3 ),
   error('Too many output arguments');
end

if ~compressed,
   n = fix(size(A,1)/2);
   if (size(A,1) ~= 2*n) || (size(A,2) ~= 2*n),
      error('H is not a 2n-by-2n matrix');
   end
   [A,QG] = haconv(A);
else
   n = size(A,1);
   if (size(A,2) ~= n),
      error('A is not a square matrix');
   end
   if (size(QG,1) ~= n) || (size(QG,2) ~= n+1),
      error('QG is not an n-by-(n+1) matrix');
   end
end

if ( meth == 2 || meth == 4 ),
   balanc = 'both';
else
   balanc = 'permute';
end
[X,A,QG] = habalance(A,QG,balanc);

if ( meth <= 2 ),
   [U,V,S,T,e] = haurvps(A,QG,0);
   if nargout <= 2,
      US = hapack_haeig(9,S,T,U(1:n,1:n), U(1:n,n+1:2*n),V(1:n,1:n), V(1:n,n+1:2*n));
   else
      [US,UU] = hapack_haeig(9,S,T,U(1:n,1:n), U(1:n,n+1:2*n),V(1:n,1:n), V(1:n,n+1:2*n));
   end
else
   [U,V,S,T,G,e] = haurvps(A,QG);
   if nargout <= 2,
      US = hapack_haeig(9,S,T,G,U(1:n,1:n), U(1:n,n+1:2*n),V(1:n,1:n), V(1:n,n+1:2*n));
   else
      [US,UU] = hapack_haeig(9,S,T,G,U(1:n,1:n), U(1:n,n+1:2*n),V(1:n,1:n), V(1:n,n+1:2*n));
   end
end

US = X*US;
if ( nargout > 2 ),
   UU = X*UU;
elseif nargout == 2,
   UU = e;
end
