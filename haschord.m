function [U,A,QG] = haschord(U,A,QG,index)
% HASCHORD  Reorders Schur form of a Hamiltonian matrix.
%    [Uo,Ho] = haschord(Ui,Hi,Index) finds, for a given real
%    2*n-by-2*n Hamiltonian matrix Hi = [Ai, Gi; 0, -Ai']
%    with Ai in Real Schur form, an orthogonal symplectic matrix
%    U so that the first n eigenvalues appearing on the diagonal
%    of Hi are ordered according to the increasing values of the
%    array Index where the i-th element of Index corresponds to
%    the eigenvalue appearing as the element Hi(i,i). For a
%    2-by-2 diagonal block only the minimum of the corresponding
%    Index pair is considered. On output, Ho = U'*Hi*U is again
%    in Hamiltonian Schur form and Uo = Ui*U.
%
%    Example: With Index = [ 2 2 1 3 1 1 ] the fifth and third
%    eigenvalues of the 6-by-6 matrix Hi are reordered to the top
%    of the matrix. 
%
%    [Uo,Ho] = HASCHORD(Hi,Index) initializes Ui to the identity
%    matrix.
%
%    [Uo,Ao,QGo] = HASCHORD(Ui,Ai,QGi,Index) requires Ai and QGi
%    contain the matrix Hi in compressed format and returns Ao
%    and QGo containing the matrix Ho in compressed format.
%
%           [Ho] = HASCHORD(Hi,Index)
%    [Uo,Ao,QGo] = HASCHORD(Ai,QGi,Index)
%       [Ao,QGo] = HASCHORD(Ai,QGi,Index)
%
%    See also SCHORD.

%    This is a MATLAB gateway routine to the SLICOT routine
%    MB03TD.f.
%
%    RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%    Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%    Based on HAPACK package, http://www.tu-chemnitz.de/mathematik/hapack/
%
%    Revision
%    V. Sima, Oct. 2008, Mar. 2009, Nov. 2012.

if nargin < 2,
   error('At least two input arguments required.');
end

if  nargin<=2 || ( nargin==3 && nargout<=2 && size(A,1)==size(A,2) ),
   compressed = 0;
else
   compressed = 1;
end
%
% Check input arguments.
%
if ( nargin == 2 )
   n = floor(size(U,1) / 2);
   if (size(U,1) ~= 2*n) || (size(U,2) ~= 2*n),
      error('H is not a 2n-by-2n matrix');
   end
elseif (nargin == 3 && ~compressed),
   n = floor(size(A,1) / 2);
   if (size(A,1) ~= 2*n) || (size(A,2) ~= 2*n),
      error('H is not a 2n-by-2n matrix');
   end
   if (size(U,1) ~= 2*n) || (size(U,2) ~= 2*n),
      error('U is not a 2n-by-2n matrix');
   end
elseif (nargin == 3)
   n = size(U,1);
   if (size(U,2) ~= n),
      error('A is not a square matrix');
   end
   if (size(A,1) ~= n) && (size(A,2) ~= n+1),
      error('QG is not an n-by-(n+1) matrix');
   end
elseif (nargin == 4)
   n = size(A,1);
   if (size(A,2) ~= n),
      error('A is not a square matrix');
   end
   if (size(QG,1) ~= n) && (size(QG,2) ~= n+1),
      error('QG is not an n-by-(n+1) matrix');
   end
   if (size(U,1) ~= 2*n) || (size(U,2) ~= 2*n),
      error('U is not a 2n-by-2n matrix');
   end
else
   error('Too many input arguments.');
end
if ( nargout > 3 ) || ( ~compressed && nargout > 2 ),
   error('Too many output arguments.');
end
if nargin == 3,
   index = QG;
elseif nargin < 3,
   index = A;
end
if min(size(index)) > 1 || length(index) ~= 2*n,
   error('Index must be a vector of length 2*n.');
end
if ~isreal(index),
   error('Index must be a real vector.');
end
if size(index,1) > size(index,2),
   index = index';
end
if ~compressed,
   if nargin == 3,
      [A,QG] = haconv(A);
   else
      [A,QG] = haconv(U);
   end
end
initu = 1;
if nargout == 3 || ( ~compressed && nargout >= 2),
   compu = 1;
   if (nargin == 2) || (compressed && nargin == 3),
      U1 = eye(n);
      U2 = zeros(n);
   else
      U1 = U(1:n,1:n);
      U2 = U(1:n,n+1:2*n);
      initu = 0;
   end
else
   compu = 0;
end
if compressed && initu,
   QG = A;
   A = U;
end
%
% Make complex eigenvalues look equally good.
%
for j = 1:n-1,
   if A(j+1,j),
      index([j j+1]) = min(index([j j+1]));
      index([n+j n+j+1]) = min(index([n+j n+j+1]));
   end
end
%
% Do reordering. In each step the maximal number of eigenvalues that can
% be treated by one call of MB03TD is reordered.
%
while true,
   myindex = index;
   select  = zeros(1,n);
   lower   = zeros(1,n);
   for j = 1:n,
      [~,i] = min(myindex);
      if i > n,
         i = i-n;
         islower = 1;
      else
         islower = 0;
      end
      if i < n && any(select(i+1:n)),
         break
      else
         myindex([i, n+i]) = inf;
         select(i) = 1;
         lower(i) = islower;
      end
   end
   if all(select) && ~any(lower),
      break;
   end
   if compu,
      [U1,U2,A,QG] = hapack_haeig(6,U1,U2,A,QG,select,lower);
   else
      [A,QG] = hapack_haeig(6,A,QG,select,lower);
   end
   index = [index(n+1:2*n).*lower + index(1:n).*~lower, ...
            index(1:n).*lower + index(n+1:2*n).*~lower ];
   index = [index(select>0), index(~select), ...
            index(n+find(select)), index(n+find(~select))];
end
if compu,
   U = [ U1, U2; -U2, U1 ];
end
if ~compressed,
   if compu,
      A = haconv(A,QG);
   else
      U = haconv(A,QG);
   end
end
if compressed && ~compu,
   U = A;
   A = QG;
end
