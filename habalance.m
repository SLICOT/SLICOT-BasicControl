function [T,A,QG] = habalance(A,QG,job)
% HABALANCE  Symplectic scaling to improve eigenvalue accuracy.
%    [T,R] = HABALANCE(H) computes a symplectic similarity
%    transformation T for a Hamiltonian matrix H so that
%    R = T\H*T has, as nearly as possible, equal row and column
%    norms.  T is a permutation of a diagonal matrix whose
%    elements are signed integer powers of two so that the
%    balancing does not introduce any round-off error.
%    The matrix H may be complex.
%
%    [T,Ar,QGr] = HABALANCE(A,QG) requires A and QG contain the
%    matrix H in compressed format and returns Ar and QGr
%    containing the matrix R in compressed format.
%    The matrices A and QG may be complex.
%
%    R = HABALANCE(H) ( [Ar,QGr] = HABALANCE(A,QG) ) returns only
%    the balanced matrix R (in compressed format).
%
%    An additional parameter may be used to specify the operations
%    performed on H.
%    HABALANCE(H,'permute') only permutes H to bring it closer to
%    triangular form.
%    HABALANCE(H,'scale') only scales H to equalize its row and
%    column norms.
%    HABALANCE(H,'both') does both, permuting and scaling. This is
%    the default.
%    Similarly for HABALANCE(A,QG,'...').
%
%    See also BALANCE, SHBALANCE.

%    This is a MATLAB gateway routine to the SLICOT routines
%    MB04DD.f and MB04DZ.f.
%
%    RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%    Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%    Partly based on HAPACK package, http://www.tu-chemnitz.de/mathematik/hapack/
%
%    Revision
%    V. Sima, Oct. 2008, Mar. 2009, Oct. 2012.

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
      error('H is not a 2n-by-2n matrix');
   end
   [A,QG] = haconv(A);
else
   n = size(A,1); 
end
if (size(A,2) ~= n),
   error('A is not a square matrix');
end
if (size(QG,1) ~= n) || (size(QG,2) ~= n+1),
   error('QG is not an n-by-(n+1) matrix');
end

realH = isreal(A) && isreal(QG);
if ~realH,
   jbal = -1;
   if isequal(job,'b'),
      ibal = 3;
   elseif isequal(job,'s'),
      ibal = 2;
   elseif isequal(job,'p'),
      ibal = 1;
   else
      ibal = 0;
   end
end
if nargout >= 2,
   if realH,
      [A,QG,scale,ilo] = hapack_haeig(4,A,QG,job);
   else
      [ilo,scale,A,QG] = HaeigZ(A,QG,jbal,ibal);
   end
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
   if realH,
      [T,A,~,~] = hapack_haeig(4,A,QG,job);
   else
      [~,~,T,A] = HaeigZ(A,QG,jbal,ibal);
   end
end

if ~compressed,
   if nargout >= 2,
      A = haconv(A,QG);
   else
      T = haconv(T,A);   
   end
elseif nargout == 2,
    T = A;  A = QG;
end
