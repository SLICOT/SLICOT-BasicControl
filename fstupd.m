function [R,G,H,CS,L] = fstupd(T,R,G,H,CS,L)
%FSTUPD  Factor and/or update a factorization of a symmetric positive 
%        definite (block) Toeplitz matrix BT or BTA and solve associated
%        linear systems, given the first (block) row / column T of BT or
%        [ T Ta ] / [ T  ] of BTA.
%                   [ Ta ] 
%
%        [R,G,H,CS] = FSTUPD(T)  computes the Cholesky factor R, i.e., 
%        R is upper / lower triangular, so that R'*R = BT / R*R' = BT,
%        the generator G of inv(BT), and deliver the information needed
%        to update the Cholesky factorization.
%
%        [R,G,H,CS,L] = FSTUPD(T)  also computes the Cholesky factor L 
%        of inv(BT), i.e., L is lower / upper triangular, so that
%        L'*L = inv(BT) / L*L' = inv(BT).
%
%        [Ru,Gu,Hu,CSu]    = FSTUPD(Ta,R,G,H,CS)    or  
%        [Ru,Gu,Hu,CSu,Lu] = FSTUPD(Ta,R,G,H,CS,L)  updates the Cholesky  
%        factor R, the generator G, the information needed to further
%        update the Cholesky factorization, and possibly the Cholesky
%        factor L, given additional blocks of data, Ta, and previous 
%        factorization results.
%        Note: the output parameters G or Gu above are appropriate for
%              further updating. To get the true generator, the submatrix 
%              (2,1) / (1,2) of G or Gu must be set to zero.
%
%        [Ru,X] = FSTUPD(Ta,R,G,H,CS,B)  also computes the solution of
%        X*BTA = B / BTA*X = B.
%
%        [Ru,Gu,X] = FSTUPD(Ta,R,G,H,CS,B)  also computes the updated 
%        generator.
%
%        See also FSTOEP, FSTCHOL, FSTGEN, FSTSOL
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 20-12-2000. 
%        Revised: Mar. 2009. 
%

ni = nargin;  nout = nargout;
%
if ni == 1,
   %
   % Initial factorization.
   %
   if nout == 2,
      [G,R] = fstoep(6,T);
   elseif nout == 3,
      [G,R,H] = fstoep(6,T);
   elseif nout == 4,
      [G,R,H,CS] = fstoep(6,T);
   elseif nout == 5,
      [G,L,R,H,CS] = fstoep(7,T);
   else
      error('Improper number of output arguments.')
   end
elseif ni >= 5 && nout >= 4,
   %
   % Factorization updating.
   %
   if nout == 4 || ni == 5,
      [G,R,H,CS] = fstoep(10,T,H,CS,G,R);
   elseif nout == 5,
      [G,L,R,H,CS] = fstoep(10,T,H,CS,G,R,L);
   else
      error('Improper number of output arguments.')
   end
elseif ni == 6 && nout <= 3,
   %
   % Factorization updating and solution.
   %
   if nout == 2,
      %
      % The second output argument below is actually the solution X.
      %
      [R,G] = fstoep(8,T,H,CS,G,R,L);
   elseif nout == 3,
      %
      % The third output argument below is actually the solution X.
      %
      [G,R,H] = fstoep(9,T,H,CS,G,R,L);
   else
      error('Improper number of output arguments.')
   end
else
   error('Improper number of input arguments.')
end
%
% end fstupd
