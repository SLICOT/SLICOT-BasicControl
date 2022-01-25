function [U,V,S,T,G] = haurv(A,QG)
% HAURV  Symplectic URV form of a general 2n-by-2n matrix H.
%    R = HAURV(H) produces the symplectic URV form
%              R = [ T, G; ...
%                    0, S' ]
%    of the same dimension as H, where S is in upper Hessenberg
%    form and T is upper triangular.
%
%    [U,V,R] = HAURV(H) produces orthogonal symplectic matrices
%    U and V so that H = U*R*V'.
%
%    [S,T,G] = HAURV(A,QG) assumes that H is Hamiltonian,
%              H = [ A, G; ...
%                    Q, -A' ],
%    and requires A and QG contain the matrix H in compressed format.
%    It returns the reduced matrix R in terms of its submatrices
%    S, T and G.
%
%    [U,R] = HAURV(H)
%    [U,S,T,G] = HAURV(A,QG)
%    [U,V,S,T,G] = HAURV(A,QG)
%
%    See also HAPVL, HAURVPS, HESS.

%    This is a MATLAB gateway routine to the HAPACK routines
%    MB04TB.f and MB04WR.f.
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
   if nargout <= 1,
      [U] = hapack_haeig(2,A);
   elseif nargout <= 2,
      [U1,U2,V] = hapack_haeig(2,A);
      U = [U1, U2; -U2, U1];
   elseif nargout <= 3,
      [U1,U2,V1,V2,S] = hapack_haeig(2,A);
      U = [U1, U2; -U2, U1];
      V = [V1', V2'; -V2', V1'];
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
   if nargout <= 3,
      [U,V,S] = hapack_haeig(2,A,QG);
   elseif nargout <= 4,
      [U1,U2,V,S,T] = hapack_haeig(2,A,QG);
      U = [U1, U2; -U2, U1];
   elseif nargout <= 5,
      [U1,U2,V1,V2,S,T,G] = hapack_haeig(2,A,QG);
      U = [U1, U2; -U2, U1];
      V = [V1, V2'; -V2', V1];
   else
      error('Too many output arguments.');
   end
end
