function [X] = fstmul(TC,TR,B)
%FSTMUL  Compute the matrix-vector products X = T*B for a block Toeplitz
%        matrix T, given the first block column TC and the first block
%        row TR of T.
%
%        X = FSTMUL(TC,TR,B)
%
%        See also FSTLSQ, FSTQR.
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        D. Kressner 01-08-2002.
%        Revised -
%

ni = nargin;  nout = nargout;
%
if ni ~= 3,
   error('Improper number of input arguments')
end
if nout > 1,
   error('Improper number of output arguments')
end

X = fstoeq(13,TC,TR,B);

%
% end fstmul
