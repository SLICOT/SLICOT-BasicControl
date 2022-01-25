function e = haeig(A,QG,flag)
% HAEIG  Eigenvalues of a Hamiltonian matrix.
%    E = HAEIG(H) is a vector containing the eigenvalues of
%    a 2*n-by-2*n Hamiltonian matrix H = [A, G; Q, -A']. If H
%    is real the eigenvalues are ordered so that E(i) = -E(n+i)
%    and the leading n elements of E have nonpositive real part.
%
%    E = HAEIG(H,'nobalance') performs the computation without
%    symplectic balancing.
%
%    E = HAEIG(A,QG{,'nobalance'}) requires A and QG contain
%    the matrix H in compressed format.
%   
%    Note that HAEIG supports complex Hamiltonian matrices.
%
%    See also EIG, HASTAB, HAURVPS.

%    This is a MATLAB gateway routine to the SLICOT routines
%    MB03XD.f and MB03XZ.f.
%
%    RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%    Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%    Partly based on HAPACK package, http://www.tu-chemnitz.de/mathematik/hapack/
%
%    Revision
%    V. Sima, 2007, 2008, Mar. 2009, Oct. 2012.

if nargin < 1,
   error('At least one input argument required.');
elseif (nargin == 1),
   balanc = 'b';
   ibal   = 3;
   compressed = 0;
elseif (nargin == 2) && (ischar(QG)),
   if ~(isequal(lower(QG),'nobalance')),
      error('String argument is an unknown option');
   end
   balanc = 'n';
   ibal   = 0;
   compressed = 0;
elseif nargin <= 2,
   balanc = 'b';
   ibal   = 3;
   compressed = 1;
elseif nargin <= 3,
   if ~(isequal(lower(flag),'nobalance')),
      error('String argument is an unknown option');
   end
   balanc = 'n';
   ibal   = 0;
   compressed = 1;   
else
   error('Too many input arguments');
end

if ( nargout > 1 ),
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

if isreal(A) && isreal(QG),
   [wr,wi] = hapack_haeig(5,A,QG,balanc);
   e = wr+1i*wi;
   e = [e;-e];
else
   e = HaeigZ(A,QG,0,0,ibal);
end

