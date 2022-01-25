function [U,V,S,T,G,w] = haurvps(A,QG,flag)
% HAURVPS  Symplectic URV/periodic Schur form of a Hamiltonian
%          matrix H.
%
%    [w] = HAURVPS(H) computes the (stable, if H is real) 
%    eigenvalues of H.
%
%    If H is real, then:
%    [R,w] = HAURVPS(H) produces the symplectic URV/periodic
%    Schur form
%              R = [ T, G; ...
%                    0, S' ]
%    of the same dimension as H, where S is in real Schur form and
%    T is upper triangular. The vector w contains the nonpositive
%    square roots of the eigenvalues of -S*T, which are the
%    stable eigenvalues of H.
%
%    [U,V,R,w] = HAURVPS(H) produces orthogonal symplectic
%    matrices U and V so that H = U*R*V'.
%
%    [S,T,G,w] = HAURVPS(A,QG) requires A and QG contain the
%    matrix H in compressed format and returns the reduced matrix
%    R in terms of its submatrices S, T and G.
%
%    [w] = HAURVPS(H)
%    [R,w] = HAURVPS(H)
%    [U,R,w] = HAURVPS(H)
%    [U,V,R,w] = HAURVPS(H)
%    [w] = HAURVPS(A,QG)
%    [S,T,G,w] = HAURVPS(A,QG)
%    [U,S,T,G,w] = HAURVPS(A,QG)
%    [U,V,S,T,G,w] = HAURVPS(A,QG)
%    [S,T,w] = HAURVPS(A,QG,0)
%    [U,S,T,w] = HAURVPS(A,QG,0)
%    [U,V,S,T,w] = HAURVPS(A,QG,0)
%
%    If H is complex, then:
%    [R,w] = HAURVPS(H) produces the structured complex Schur form of
%    the matrix i*He,
%           R = U'*(i*He)*U = [ S,  G; ...
%                               0, -S' ]
%    where S is upper triangular, G is Hermitian, and He is a real
%    skew-Hamiltonian matrix of order twice the order of H, built using
%    the real and imaginary parts of H. (See the SLICOT Library routine
%    MB03XZ, for details.) The vector w contains the eigenvalues of S,
%    which are also the eigenvalues of H.
%
%    [U,R,w] = HAURVPS(H) also produces the unitary symplectic matrix U.
%
%    [U,S,G,w] = HAURVPS(A,QG) requires A and QG contain the
%    matrix H in compressed format.
%
%    [Ur,E,w] = HAURVPS(H,'real') or  [Ur,Sr,Gr,w] = HAURVPS(A,QG,'real')
%    produces the structured real Schur form of the matrix He,
%           E = Ur'*He*Ur = [ Sr, Gr; ...
%                              0, Sr' ]
%    where Sr is in real Schur form, Gr is skew-symmetric, and Ur is real
%    orthogonal symplectic.
%
%    [w] = HAURVPS(H)
%    [R,w] = HAURVPS(H)
%    [U,R,w] = HAURVPS(H)
%    [w] = HAURVPS(A,QG)
%    [S,G,w] = HAURVPS(A,QG)
%    [U,S,G,w] = HAURVPS(A,QG)
%    [S,w] = HAURVPS(A,QG,0)
%    [U,S,w] = HAURVPS(A,QG,0)
%    [E,w] = HAURVPS(H,'real')
%    [Ur,E,w] = HAURVPS(H,'real')
%    [Sr,Gr,w] = HAURVPS(A,QG,'real')
%    [Ur,Sr,Gr,w] = HAURVPS(A,QG,'real')
%
%    See also HAEIG, HAURV, SCHUR.

%    This is a MATLAB gateway routine to the SLICOT routines MB03XD.f
%    and MB03XZ.f.
%
%    RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%    Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%    Partly based on HAPACK package, http://www.tu-chemnitz.de/mathematik/hapack/
%
%    Revisions:
%    V. Sima, November 2008, Mar. 2009, Oct. 2012.


if nargin < 1,
   error('At least one input argument required.');
end
if nargin > 3,
   error('Too many input arguments.');
end

if nargin > 1 && nargout > 1,
   compg = 1;
else
   compg = 0;
end

realFactHe = nargin == 2 && ischar(QG)   && ( ~isreal(A) || numel(A) == 0 ) || ...
             nargin == 3 && ischar(flag) && ( ~( isreal(A) && isreal(QG) ) || numel(A) == 0 );

if nargin == 3 && ~realFactHe,
   if flag == 0,
      compg = 0;
   else
      error('Use haurvps(A,QG,0) if G is not to be computed.');
   end
end
   
if nargin == 1 || ( nargin == 2 && realFactHe ),
   n = floor(size(A,1) / 2);
   if (size(A,1) ~= 2*n) || (size(A,2) ~= 2*n),
      error('H is not a 2n-by-2n matrix');
   end
else
   n = size(A,1);
   if (size(A,2) ~= n),
      error('A is not a square matrix');
   end
   if (size(QG,1) ~= n) || (size(QG,2) ~= n+1),
      error('QG is not an n-by-(n+1) matrix');
   end
end

if n == 0,
   U = [];
   if nargout > 1,
      V = [];
      if nargout > 2,
         S = [];
         if nargout > 3,
            T = [];
            if nargout > 4,
               G = [];
               if nargout > 5,
                  w = [];
               end
            end
         end
      end
   end
   return
end

if nargin == 1,
   realH  = isreal(A);
   [A,QG] = haconv(A);
   if realH,
      if nargout <= 1,
         [wr,wi] = hapack_haeig(5,A,QG,'b');
         U = wr + 1i*wi;
      elseif nargout <= 2,
         [Sx,Tx,Gx,wr,wi] = hapack_haeig(3,A,QG);
         U = [ Tx, Gx; zeros(n), Sx' ];
         V = wr + 1i*wi;
      elseif nargout <= 3,
         [U1,U2,Sx,Tx,Gx,wr,wi] = hapack_haeig(3,A,QG);
         U = [U1, U2; -U2, U1];
         V = [ Tx, Gx; zeros(n), Sx' ];
         S = wr + 1i*wi;
      elseif nargout <= 4,
         [U1,U2,V1,V2,Sx,Tx,Gx,wr,wi] = hapack_haeig(3,A,QG);      
         U = [U1, U2; -U2, U1];
         V = [V1, V2; -V2, V1];
         S = [ Tx, Gx; zeros(n), Sx' ];
         T = wr + 1i*wi;
      else
         error('Too many output arguments.');
      end
   else
      balanc = 3;  n2 = 2*n;
      if nargout <= 1,
         U = HaeigZ(A,QG,0,0,balanc);
      elseif nargout <= 2,
         [V,Tx,Gx] = HaeigZ(A,QG,2,0,balanc);
         U = [ Tx, Gx; zeros(n2), -Tx' ];
      elseif nargout <= 3,
         [S,Tx,Gx,U1,U2] = HaeigZ(A,QG,2,1,balanc);
         U = [U1, U2; -U2, U1];
         V = [ Tx, Gx; zeros(n2), -Tx' ];
      else
         error('Too many output arguments.');
      end
   end
elseif ~realFactHe
   realH = isreal(A) && isreal(QG);
   if realH,
      if compg,
         if nargout <= 4,
            [U,V,S,wr,wi] = hapack_haeig(3,A,QG);
            T = wr + 1i*wi;
         elseif nargout <= 5,
            [U1,U2,V,S,T,wr,wi] = hapack_haeig(3,A,QG);
            U = [U1, U2; -U2, U1];
            G = wr + 1i*wi;
         elseif nargout <= 6,
            [U1,U2,V1,V2,S,T,G,wr,wi] = hapack_haeig(3,A,QG);
            U = [U1, U2; -U2, U1];
            V = [V1, V2; -V2, V1];
            w = wr + 1i*wi;
         else
            error('Too many output arguments.');
         end
      else
         if nargout == 1,
            [wr,wi] = hapack_haeig(5,A,QG,'b');
            U = wr + 1i*wi;
         elseif nargout <= 3,
            [U,V,wr,wi] = hapack_haeig(3,A,QG);
            S = wr + 1i*wi;
         elseif nargout <= 4,
            [U1,U2,V,S,wr,wi] = hapack_haeig(3,A,QG);
            U = [U1, U2; -U2, U1];
            T = wr + 1i*wi;
         elseif nargout <= 5,
            [U1,U2,V1,V2,S,T,wr,wi] = hapack_haeig(3,A,QG);
            U = [U1, U2; -U2, U1];
            V = [V1, V2; -V2, V1];
            G = wr + 1i*wi;
         else
            error('Too many output arguments.');
         end
      end
   else
      balanc = 3; 
      if compg,
         if nargout <= 3,
            [S,U,V] = HaeigZ(A,QG,2,0,balanc);
         elseif nargout <= 4,
            [T,V,S,U1,U2] = HaeigZ(A,QG,2,1,balanc);
            U = [U1, U2; -U2, U1];
         else
            error('Too many output arguments.');
         end
      else
         if nargout <= 1,
            U = HaeigZ(A,QG,0,0,balanc);
         elseif nargout <= 2,
            [V,U] = HaeigZ(A,QG,1,0,balanc);
         elseif nargout <= 3,
            [S,V,U1,U2] = HaeigZ(A,QG,1,1,balanc);
            U = [U1, U2; -U2, U1];
         else
            error('Too many output arguments.');
         end
      end
   end
else
   % Complex case, real factorization only.
   balanc = 3; 
   if nargin == 2,
      [A,QG] = haconv(A);  n2 = 2*n;
      if nargout <= 2,
         [V,Tx,Gx] = HaeigZ(A,QG,0,0,balanc);
         U = [ Tx, Gx; zeros(n2), Tx' ];
      elseif nargout <= 3,
         [S,Tx,Gx,U1,U2] = HaeigZ(A,QG,0,1,balanc);
         U = [U1, U2; -U2, U1];
         V = [ Tx, Gx; zeros(n2), Tx' ];
      else
         error('Too many output arguments.');
      end
   else
      if nargout <= 3,
         [S,U,V] = HaeigZ(A,QG,0,0,balanc);
      elseif nargout <= 4,
         [T,V,S,U1,U2] = HaeigZ(A,QG,0,1,balanc);
         U = [U1, U2; -U2, U1];
      else
         error('Too many output arguments.');
      end
   end
end
