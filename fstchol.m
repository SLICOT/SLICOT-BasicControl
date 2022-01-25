function [R,G,L,X] = fstchol(T,B)
%FSTCHOL Factor a symmetric positive definite (block) Toeplitz matrix BT
%        and solve associated linear systems, given the first (block) 
%        row / column T of BT.
% 
%        R = FSTCHOL(T)  computes the Cholesky factor R, i.e., R is
%        upper / lower triangular, so that R'*R = BT / R*R' = BT.
%
%        [R,G] = FSTCHOL(T)  also computes the generator G of inv(BT).
%
%        [R,G,L] = FSTCHOL(T)  also computes the Cholesky factor L
%        of inv(BT), i.e., L is lower / upper triangular, so that
%        L'*L = inv(BT) / L*L' = inv(BT).
%
%        [R,X]     = FSTCHOL(T,B)  also solve the linear systems 
%        [R,G,X]   = FSTCHOL(T,B)       X*BT = B / BT*X = B.
%        [R,G,L,X] = FSTCHOL(T,B)  
%
%        See also FSTOEP, FSTGEN, FSTSOL, FSTUPD
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
      R = fstoep(1,T);
   elseif nout == 2,
      [G,R] = fstoep(3,T);
   elseif nout == 3,
      [G,L,R] = fstoep(5,T);
   else
      error('Improper number of output arguments.')
   end
elseif ni == 2,
   %
   % Factorization and solution.
   %
   if nout == 2,
      %
      % The second output argument below is actually the solution X.
      %
      [R,G] = fstoep(1,T,B);
   elseif nout == 3,
      %
      % The third output argument below is actually the solution X.
      %
      [G,R,L] = fstoep(3,T,B);
   elseif nout == 4,
      [G,L,R,X] = fstoep(5,T,B);
   else
      error('Improper number of output arguments.')
   end
else
   error('Improper number of input arguments.')
end
%
% end fstchol
