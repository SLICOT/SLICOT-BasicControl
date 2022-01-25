function [T,A,QG] = shbalance(A,QG,job)
% SHBALANCE  Symplectic scaling to improve eigenvalue accuracy.
%    [T,R] = SHBALANCE(W) computes a symplectic similarity
%    transformation T for a skew-Hamiltonian matrix W so that
%    R = T\W*T has, as nearly as possible, equal row and column
%    norms.  T is a permutation of a diagonal matrix whose
%    elements are signed integer powers of two so that the
%    balancing does not introduce any round-off error.
%
%    [T,Ar,QGr] = SHBALANCE(A,QG) requires A and QG contain the
%    matrix W in compressed format and returns Ar and QGr
%    containing the matrix R in compressed format.
%
%    R = SHBALANCE(W) ( [Ar,QGr] = SHBALANCE(A,QG) ) returns only
%    the balanced matrix R (in compressed format).
%
%    An additional parameter may be used to specify the operations
%    performed on W.
%    SHBALANCE(W,'permute') only permutes W to bring it closer to
%    triangular form.
%    SHBALANCE(W,'scale') only scales W to equalize its row and
%    column norms.
%    SHBALANCE(W,'both') does both, permuting and scaling. This is
%    the default.
%    Similarly for SHBALANCE(A,QG,'...').
%
%    See also BALANCE, HABALANCE.

%    This is a MATLAB gateway routine to the SLICOT routine
%    MB04DS.f.
%
%    RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%    Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%    Based on HAPACK package, http://www.tu-chemnitz.de/mathematik/hapack/
%
%    Revision
%    V. Sima, Oct. 2008, Mar. 2009.

if nargin < 1,
   error('At least one input argument required.');
elseif (nargin == 1),
   job = 'both';
   compressed = 0;
elseif (nargin == 2) && (~ischar(QG)),
   job = 'both';
   compressed = 1;
elseif nargin <= 2,
   job = QG;
   compressed = 0;
elseif nargin <= 3,
   compressed = 1;
else
   error('Too many input arguments');
end

if ( compressed && (nargout > 3) ) || ...
   (~compressed && (nargout > 2) ),
   error('Too many output arguments');
end

job = lower(job);

if ~(isequal(job,'both') || isequal(job,'permute') || isequal(job,'scale') ),
   error('String argument is an unknown option');
end
job = job(1);

if ~compressed,
   n = fix(size(A,1)/2);
   if (size(A,1) ~= 2*n) || (size(A,2) ~= 2*n),
      error('W is not a 2n-by-2n matrix');
   end
   [A,QG] = shconv(A);
end
n = size(A,1);
if (size(A,2) ~= n),
   error('A is not a square matrix');
end
if (size(QG,1) ~= n) || (size(QG,2) ~= n+1),
   error('QG is not an n-by-(n+1) matrix');
end

if nargout >= 2,
   [A,QG,scale,ilo] = slicot_sheig(3,A,QG,job);
   d = [ones(ilo-1,1);scale(ilo:n)];
   T = diag([d;1./d]);
   for i = ilo-1:-1:1,
      sci = scale(i);
      if ( sci > n ),
         T(:,[i, n+i]) = [T(:,n+i), -T(:,i)];
         sci = sci - n;
      end
      T(:,[i, sci, n+i, n+sci]) = T(:,[sci, i, n+sci, n+i]);
   end
   T = T';
else
   [T,A,scale,ilo] = slicot_sheig(3,A,QG,job);
end

if ~compressed,
   if nargout >= 2,
      A = shconv(A,QG);
   else
      T = shconv(T,A);
   end
elseif nargout == 2,
    T = A;  A = QG;
end
