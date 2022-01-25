function [beta,omega] = slstabr(A,tol,interval)
%SLSTABR  Complex stability radius of a matrix.
%
%        BETA = SLSTABR(A)  returns beta(A), an estimate of the
%        2-norm distance from a real matrix A to the nearest 
%        complex matrix with an eigenvalue on the imaginary axis.
%        The number beta(A) is the minimum of the smallest
%        singular value of the matrix (A - jwI), where I is the
%        identity matrix and j**2 = -1, and the minimum is taken
%        over all real w. The computed BETA is an upper bound
%        for beta(A). If all eigenvalues of A lie in the open 
%        left half complex plane, then beta(A) is the distance 
%        to the nearest unstable complex matrix, i.e., the
%        complex stability radius.
%
%        [BETA,OMEGA] = SLSTABR(SYS,TOL)  has additional input
%        and output arguments:
%
%        TOL specifies the accuracy with which beta(A) is to be
%        calculated. If TOL is less than eps, than eps is used
%        instead.
%        Default: TOL = eps.
%
%        OMEGA is the value of w such that the smallest singular
%        value of (A - jwI) equals beta(A).
%
%        [LOW,HIGH] = SLSTABR(A,TOL,INTERVAL)  for INTERVAL = 1, 
%        returns a lower and an upper bound for beta(A).
%        OMEGA is not returned.
%
%        INTERVAL specifies the desired results:
%        INTERVAL = 0 :  compute an upper bound BETA on beta(A);
%        INTERVAL = 1 :  compute upper and lower bounds on beta(A).
%        Setting INTERVAL = 1 usually implies a faster run.
%        Default: INTERVAL = 0.
%
%        If INTERVAL = 1, TOL specifies the accuracy with which
%        LOW and HIGH approximate beta(A). If TOL is less than
%        sqrt(eps), then that value is used instead.
%        The recommended value is TOL = 9, which gives an estimate
%        of beta(A) correct to within an order of magnitude.
%        Default: TOL = sqrt(eps).

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, Dec. 2002.
%
%        Revisions:
%        -
%

if nargin < 1,
    error('Usage: BETA = SLSTABR(A,TOL,INTERVAL)')
elseif nargin < 2,
   tol = 0;
   interval = 0;
elseif nargin < 3,
   interval = 0;
end

if interval > 0
   job = 3;
else
   job = 4;
end

[beta,omega] = Hnorm(job,A,tol);
%
% end slstabr
