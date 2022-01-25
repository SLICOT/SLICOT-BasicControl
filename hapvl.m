function [U,A,QG] = hapvl(A,QG)
% HAPVL  Paige-Van Loan (PVL) form of a Hamiltonian matrix H.
%    R = HAPVL(H) produces the PVL form R = [Ar, Gr; Qr, -Ar']
%    of the same dimension as H = [A, G; Q, -A'], where Ar is
%    in upper Hessenberg form and Qr is diagonal.
%
%    [U,R] = HAPVL(H) produces an orthogonal symplectic matrix U
%    so that H = U*R*U'.
%
%    [Ar,QGr] = HAPVL(A,QG) requires A and QG contain the matrix
%    H in compressed format and returns Ar and QGr containing the
%    matrix R in compressed format.
%
%    [U,Ar,QGr] = HAPVL(A,QG)
%
%    See also HESS.

%    This is a MATLAB gateway routine to the SLICOT routines
%    MB04PB.f and MB04WP.f.
%
%    RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%    Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%    Based on HAPACK package, http://www.tu-chemnitz.de/mathematik/hapack/
%
%    Revision
%    V. Sima, Nov. 2008, Mar. 2009.

if nargin < 1,
   error('At least one input argument required.');
elseif nargin == 1,
   n = floor(size(A,1) / 2);
   if (size(A,1) ~= 2*n) || (size(A,2) ~= 2*n),
      error('H is not a 2n-by-2n matrix');
   end
   [A,QG] = haconv(A);
   if nargout <= 1,
      [A,QG] = hapack_haeig(1,A,QG);
      U = haconv(A,QG);
   elseif nargout <= 2,
      [U1,U2,A,QG] = hapack_haeig(1,A,QG);
      U = [U1, U2; -U2, U1];
      A = haconv(A,QG);
   else
      error('Too many output arguments.');
   end
else
   n = size(A,1);
   if (size(A,2) ~= n),
      error('A is not a square matrix');
   end
   if (size(QG,1) ~= n) || (size(QG,2) ~= n+1),
      error('QG is not an n-by-(n+1) matrix');
   end
   if nargout <= 2,
      [U,A] = hapack_haeig(1,A,QG);
   elseif nargout <= 3,
      [U1,U2,A,QG] = hapack_haeig(1,A,QG);
      U = [U1, U2; -U2, U1];
   else
      error('Too many output arguments.');
   end
end
