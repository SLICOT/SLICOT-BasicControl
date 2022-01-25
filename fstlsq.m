function [X,Y] = fstlsq(TC,TR,B,C)
%FSTLSQ  Solve linear least-squares problems min(B-T*X) or find the
%        minimum norm solution of T'*Y = C where T is a (block)
%        Toeplitz matrix with full column rank, given the first block
%        column TC and the first block row TR of T.
%
%        X = FSTLSQ(TC,TR,B)  computes a matrix X such that the
%        Euclidean norms of the columns of B-T*X are minimized.
%
%        [X,Y] = FSTLSQ(TC,TR,B,C)  also computes a matrix Y such that
%        T'*Y = C and the Euclidean norms of the columns of Y are
%        minimized.
%
%        Y = FSTLSQ(TC,TR,[],C) only computes the minimum norm solutions
%        of T'*Y = C.
%
%        See also FSTMUL, FSTQR.
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        D. Kressner 31-07-2002.
%        Revised V. Sima 17-08-2002, 03-03-2009.
%

ni = nargin;  nout = nargout;
%
if ni == 3,
   if nout <= 1,
      X = fstoeq(10,TC,TR,B);
   else
      error('Improper number of output arguments.')
   end
elseif ni == 4,
   if (isempty(B)) && nout <= 1,
      X = fstoeq(11,TC,TR,C);
   elseif nout == 2,
      [X,Y] = fstoeq(12,TC,TR,B,C);
   else
      error('Improper number of output arguments.')
   end
else
   error('Improper number of input arguments.')
end
%
% end fstlsq
