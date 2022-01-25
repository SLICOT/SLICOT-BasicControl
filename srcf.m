function [S1,Rin,K,rcnd] = srcf( sys,S,Q,R,Hess,tol )
%SRCF    Computes a combined measurement and time update of one
%        iteration of the time-varying or time-invariant Kalman filter
%        in the Square Root Covariance Form.
%
%        S1 = SRCF(SYS,S,Q,R)  computes the updated value (at
%        instant i), S1, of the square root (left Cholesky factor),
%        S, of the state covariance matrix at instant i-1.  
%        The strict upper triangular part of S is not used.
%
%        SYS is an ss object containing the following matrices (at
%        instant i): the state transition matrix A of the system; the
%        input matrix B (or the product B*sqrt(Q), if the parameter
%        Q is empty); the output matrix C. Matrix D is not used.
%
%        Q is the covariance square root (left Cholesky factor) of
%        the input (process) noise covariance matrix at instant i.
%        The strict upper triangular part is not used.
%        If Q is specified as [ ], it is assumed that the left 
%        Cholesky factor postmultiplied SYS.b.
%
%        R is the square root (left Cholesky factor) of the output
%        (measurement) noise covariance matrix at instant i.
%        The strict upper triangular part is not used.
%
%        [S1,RIN,K,RCND] = SRCF(SYS,S,Q,R,HESS,TOL)  has additional
%        input and output parameters.
%
%        HESS specifies whether the system is in observer Hessenberg 
%        form, option useful in the time-invariant case.
%        If HESS = 1, it is assumed that SYS.b contains the matrix
%        product B*sqrt(Q) and (SYS.a,SYS.c) is in lower observer
%        Hessenberg form. Otherwise, it is assumed that the system
%        is general (possibly time-varying), and the product B*sqrt(Q)
%        is not input. Use KFILTUPD for more generality.
%        Default: HESS = 0.
%        
%        TOL is a tolerance value used to test for near singularity of
%        the matrix RIN. If TOL > 0, then the given value of TOL is
%        used as a lower bound for the reciprocal condition number of
%        that matrix; a matrix whose estimated condition number is less
%        than 1/TOL is considered to be nonsingular. If TOL <= 0,
%        then ps*eps, is used instead, where ps is the product of
%        the matrix dimensions, and eps is the machine precision.
%        Default: TOL = 0.
%
%        RIN is the square root (left Cholesky factor) of the
%        covariance matrix of the innovations at instant i.
%
%        K is the Kalman filter gain matrix K at instant i.
%
%        RCND is an estimate of the reciprocal of the condition number
%        (in the 1-norm) of the matrix RIN.
%
%        See also KFILTUPD, CONVKF, SRIF
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Apr-1-2004.
%        Revised Jul-14-2004, Mar-03-2009.
%

%
ni = nargin;  nout = nargout;
if ni < 4 || nout == 0,  
    error( ['Usage: [S1]            = SRCF(SYS,S,Q,R)', sprintf('\n'), ...
            '       [S1,RIN,K,RCND] = SRCF(SYS,S,Q,R,HESS,TOL)' ] );
end
%
if ~isa( sys, 'ss' ),
    error( 'SYS must be an ss system object' );
end
A = sys.a;  B = sys.b;  C = sys.c;
if ni == 4,  
    Hess = 0;  
    tol  = 0;  
elseif ni == 5,  
    tol  = 0;  
end
if nout < 3,  opt(1) = 0;  else  opt(1) = 1;  end
if Hess == 0,  
    task = -1;  opt(2) = 0;
else  
    task = 1;   opt(2) = 1;
end
if isempty( Q ),
    opt(2) = 1;
    [ S1, K, Rin, rcnd ] = Kfiltupd( task, opt, S, A, B, C, R, tol );
else  
    [ S1, K, Rin, rcnd ] = Kfiltupd( task, opt, S, A, B, C, Q, R, tol );
end
%
% end srcf
