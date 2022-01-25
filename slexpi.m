function [F,H] = slexpi(A,delta,tol)
% SLEXPI Computes the exponential of a matrix and its integral using
%        a Pade approximation of the integral.
%
%        F = SLEXPI(A)  computes the matrix exponential, F = exp(A).
%
%        F = SLEXPI(A,DELTA)  computes F = exp(A*DELTA), where DELTA
%        is a real scalar.
%
%        [F,H] = SLEXPI(A,DELTA,TOL)  also computes the integral
%        of the matrix exponential and has additional input arguments.
%
%        TOL is a tolerance to be used in determining the order of
%        the Pade approximation to H(t), where t is a scale factor
%        determined internally.
%        Default:  TOL = sqrt(eps).
%
%        H contains an approximation to H(DELTA), where
%           H(t) =  Integral[F(s) ds] from s = 0 to s = t
%        and F(s) = exp(A*s).
%
%        See also SLMEXP, SLEXPE, SLEXPM
%
%        Comments: a Pade approximation to H(t) for some small value of t
%        (where 0 < t <= delta) is used, and then F(t) is calculated from
%        H(t). Finally, the results are re-scaled to give F(delta) and
%        H(delta).
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
    error(['Usage: F     = SLEXPI(A)',  sprintf('\n'),...
           '       [F,H] = SLEXPI(A,DELTA,TOL)'])
end
% 
if ni == 1,
    delta  = 1;
    tol = sqrt( eps );
elseif ni == 2,
    tol = sqrt( eps );
end
%
if nout == 1,
    F = slmexp( 2, A, delta, tol );
else
    [ F, H ] = slmexp( 2, A, delta, tol );
end
%
% end slexpi
