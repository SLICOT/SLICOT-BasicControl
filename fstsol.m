function [X] = fstsol(T,B)
%FSTSOL  Solve linear systems X*BT = B / BT*X = B, where BT is a 
%        symmetric positive definite (block) Toeplitz matrix,
%        given the first (block) row / column T of BT.
% 
%        X = FSTSOL(T,B)  solves the linear system 
%        X*BT = B / BT*X = B.
%
%        See also FSTOEP, FSTCHOL, FSTGEN, FSTUPD
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 20-12-2000. 
%        Revised: Mar. 2009. 
%

ni = nargin;  nout = nargout;
%
if ni == 2,
   %
   % Factorization and solution.
   %
   if nout == 1,
      X = fstoep(11,T,B);
   else
      error('Improper number of output arguments.')
   end
else
   error('Improper number of input arguments.')
end
%
% end fstsol
