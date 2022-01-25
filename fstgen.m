function [G,L,R,X] = fstgen(T,B)
%FSTGEN  Factor a symmetric positive definite (block) Toeplitz matrix BT,
%        computes the generator of inv(BT), and/or solve associated linear  
%        systems using the Cholesky factor of inv(BT), given the first 
%        (block) row / column T of BT.
% 
%        G = FSTGEN(T)  computes the generator G of inv(BT).
%
%        [G,L] = FSTGEN(T)  also computes the Cholesky factor L
%        of inv(BT), i.e., L is lower / upper triangular, so that
%        L'*L = inv(BT) / L*L' = inv(BT).
%
%        [G,L,R] = FSTGEN(T)  also computes the Cholesky factor R, i.e., 
%        R is upper / lower triangular, so that R'*R = BT / R*R' = BT.
%
%        [G,L,X]   = FSTGEN(T,B)  also solve the linear systems 
%        [G,L,R,X] = FSTGEN(T,B)       X*BT = B / BT*X = B.
%
%        See also FSTOEP, FSTCHOL, FSTSOL, FSTUPD
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 20-12-2000. 
%        Revised 03-03-2009.
%

ni = nargin;  nout = nargout;
%
if ni == 1,
   %
   % Factorization only.
   %
   if nout == 1,
      G = fstoep(2,T);
   elseif nout == 2,
      [G,L] = fstoep(4,T);
   elseif nout == 3,
      [G,L,R] = fstoep(5,T);
   else
      error('Improper number of output arguments.')
   end
elseif ni == 2,
   %
   % Factorization and solution.
   %
   if nout == 3,
      %
      % The third output argument below is actually the solution X.
      %
      [G,L,R] = fstoep(4,T,B);
   elseif nout == 4,
      [G,L,R,X] = fstoep(5,T,B);
   else
      error('Improper number of output arguments.')
   end
else
   error('Improper number of input arguments.')
end
%
% end fstgen
