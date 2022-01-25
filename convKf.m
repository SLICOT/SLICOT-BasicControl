function [P1,Rin,K,rcnd] = convKf( sys,P,Q,R,tol )
%CONVKF  Computes one recursion of the conventional Kalman filter
%        equations.
%
%        P1 = CONVKF(SYS,P,Q,R)  computes the updated value (at
%        instant i), P1, of the state covariance matrix, P, at
%        instant i-1. 
%
%        SYS is an ss object containing the system matrices at
%        instant i. Matrix D is not used.
%
%        Q is the input (process) noise covariance matrix at 
%        instant i.
%
%        R is the the output (measurement) noise covariance at
%        instant i.
%
%        [P1,RIN,K,RCND] = CONVKF(SYS,P,Q,R,TOL)  has additional
%        input and output parameters.
%
%        TOL is a tolerance value used to test for near singularity of
%        the matrix Rinov. If TOL > 0, then the given value of TOL is
%        used as a lower bound for the reciprocal condition number of
%        that matrix; a matrix whose estimated condition number is less
%        than 1/TOL is considered to be nonsingular. If TOL <= 0,
%        then ps*eps, is used instead, where ps is the product of
%        the matrix dimensions, and eps is the machine precision.
%        Default: TOL = 0.
%
%        RIN is the square root (left Cholesky factor) of the
%        covariance matrix of the innovations, Rinov, at instant i.
%
%        K is the Kalman filter gain matrix K at instant i.
%
%        RCND is an estimate of the reciprocal of the condition number
%        (in the 1-norm) of the matrix Rinov.
%
%        See also KFILTUPD, SRCF, SRIF
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Apr-1-2004.
%        Revised Jul-14-2004, Mar-2-2009.
%

%
ni = nargin;  nout = nargout;
if ni < 4 || nout == 0,  
    error( ['Usage: [P1]            = CONVKF(SYS,P,Q,R)', sprintf('\n'), ...
            '       [P1,RIN,K,RCND] = CONVKF(SYS,P,Q,R,TOL)' ] );
end
%
if ~isa( sys, 'ss' ),
    error( 'SYS must be an ss system object' );
end
A = sys.a;  B = sys.b;  C = sys.c;
if ni == 4,  tol = 0;  end
[ P1, K, Rin, rcnd ] = Kfiltupd( 0, P, A, B, C, Q, R, tol );
%
% end convKf
