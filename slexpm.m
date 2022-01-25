function [F,mdig,idig] = slexpm(A,delta,scale,ordPad)
% SLEXPM Computes the exponential of a matrix using a diagonal Pade
%        approximant with scaling and squaring.
%
%        F = SLEXPM(A)  computes the matrix exponential, F = exp(A).
%
%        F = SLEXPM(A,DELTA)  computes F = exp(A*DELTA), where DELTA
%        is a real scalar.
%
%        [F,MDIG,IDIG] = SLEXPM(A,DELTA,SCALE,ORDPAD)  has additional
%        input and output arguments.
%
%        SCALE is an integer scalar indicating whether or not the matrix
%        should be diagonally scaled:
%        SCALE = 0 :  do not scale the matrix;
%        SCALE = 1 :  diagonally scale the matrix, i.e., replace A by
%                     D*A*D**(-1), where D is a diagonal matrix chosen to
%                     make the rows and columns of A more equal in norm.
%        Default:  SCALE = 1.
%
%        ORDPAD is an integer scalar specifying the order of the
%        diagonal Pade approximant. 1 <= ORDPAD <= 15. In the absence
%        of further information, ORDPAD should be set to 9.
%        Default:  ORDPAD = 9.
%
%        MDIG is the minimal number of accurate digits in the 1-norm
%        of exp(A*DELTA).
%
%        IDIG is the number of accurate digits in the 1-norm of
%        exp(A*DELTA) at 95% confidence level.
%
%        See also SLMEXP, SLEXPE, SLEXPI
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, April 2003.
%
%        Revisions:
%        March 2009.
%

ni   = nargin;
nout = nargout;
%
if min( ni, nout ) < 1,
    error(['Usage: F             = SLEXPM(A)',  sprintf('\n'),...
           '       [F,MDIG,IDIG] = SLEXPM(A,DELTA,SCALE,ORDPAD)'])
end
% 
if ni == 1,
    delta  = 1;
    scale  = 1;
    ordPad = 9;
elseif ni == 2,
    scale  = 1;
    ordPad = 9;
elseif ni == 3,
    ordPad = 9;
end
%
if nout == 1,
    F = slmexp( 1, A, delta, scale, ordPad );
elseif nout == 2,
    [ F, mdig ] = slmexp( 1, A, delta, scale, ordPad );
else
    [ F, mdig, idig ] = slmexp( 1, A, delta, scale, ordPad );
end
%
% end slexpm
