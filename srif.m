function [Sinv1,Qin,X1,E,rcnd] = srif( sys,Sinv,Qinv,Rinv,X,Rinvy,Z,Hess,tol )
%SRIF    Computes a combined measurement and time update of one
%        iteration of the time-varying or time-invariant Kalman filter
%        in the Square Root Information Form.
%
%        SINV1 = SRIF(SYS,SINV,QINV,RINV,X,RINVY)  computes the
%        updated value (at instant i+1), SINV1, of the inverse of the
%        square root (right Cholesky factor), SINV, of the state
%        covariance matrix (hence the information square root)
%        at instant i. The strict lower triangular part of SINV
%        is not used.
%
%        SYS is an ss object containing the following matrices: 
%        the inverse AINV of the state transition matrix of the system
%        at instant i; the input matrix B (or the product AINV*B) at 
%        instant i; the output matrix C (or the product RINV*C, if the
%        parameter RINV is empty) at instant i+1. Matrix D is not used.
%
%        QINV is the inverse of the covariance square root (right Cholesky
%        factor) of the process noise (hence the information square root)
%        at instant i. The strict lower triangular part is not used.
%
%        RINV is the inverse of the covariance square root (right Cholesky
%        factor) of the output noise (hence the information square root)
%        at instant i+1. The strict lower triangular part is not used.
%        If RINV is specified as [ ], it is assumed that the right 
%        Cholesky factor premultiplied SYS.c.
%
%        X is the estimated filtered state at instant i.
%
%        RINVY is the product of RINV and the measured output vector at
%        instant i+1.
%
%        [SINV1,QIN,X1,E,RCND] = SRIF(SYS,SINV,QINV,RINV,X,RINVY,Z,...
%                                     HESS,TOL) 
%        has additional input and output parameters.
%
%        Z is the mean value of the state process noise at instant i.
%        The state process noise postmultiplies the matrix B.
%        Default: Z = zeros(m,min(m,1)), where m is the number of
%                 columns of B.
%        
%        HESS specifies whether the system is in controller Hessenberg 
%        form, option useful in the time-invariant case.
%        If HESS = 1, it is assumed that SYS.b contains the matrix
%        product AINV*B and (SYS.a,SYS.b) is in controller Hessenberg form.
%        Otherwise, it is assumed that the system is general (possibly
%        time-varying), and the product AINV*B is not input. Use KFILTUPD
%        for more generality.
%        Default: HESS = 0.
%        
%        TOL is a tolerance value used to test for near singularity of
%        the matrix SINV1. If TOL > 0, then the given value of TOL is
%        used as a lower bound for the reciprocal condition number of
%        that matrix; a matrix whose estimated condition number is less
%        than 1/TOL is considered to be nonsingular. If TOL <= 0,
%        then ps*eps, is used instead, where ps is the product of
%        the matrix dimensions, and eps is the machine precision.
%        Default: TOL = 0.
%
%        QIN is the inverse of the covariance square root (right
%        Cholesky factor) of the process noise (hence the information
%        square root) at instant i.
%
%        X1 is the estimated filtered state at instant i+1.
%
%        E is the estimated error at instant i+1.
%
%        RCND is an estimate of the reciprocal of the condition number
%        (in the 1-norm) of SINV1.
%
%        See also KFILTUPD, CONVKF, SRCF
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Apr-1-2004.
%        Revised Jul-14-2004, Mar-03-2009.
%

%
ni = nargin;  nout = nargout;
if ni < 6 || nout == 0,  
    error( ['Usage: [SINV1]               = SRIF(SYS,SINV,QINV,RINV,X,RINVY)', sprintf('\n'), ...
            '       [SINV1,QIN,X1,E,RCND] = SRIF(SYS,SINV,QINV,RINV,X,RINVY,Z,HESS,TOL)' ] );
end
%
if ~isa( sys, 'ss' ),
    error( 'SYS must be an ss system object' );
end
Ainv = sys.a;  B = sys.b;  C = sys.c;
m = size( B,2 );
if ni == 6,  
    Z    = zeros( m,min(m,1) );  
    Hess = 0;  
    tol  = 0;  
elseif ni == 7,  
    Hess = 0;  
    tol  = 0;  
elseif ni == 8,  
    tol  = 0;  
end
if nout < 3,  opt(1) = 0;  else  opt(1) = 1;  end
if Hess == 0,  
    task = -2;  opt(2) = 0;  lopt = 3;
else  
    task = 2;  lopt = 2;
end
if isempty( Rinv ),
    opt(lopt) = 1;  
    [ Sinv1, X1, Qin, E, rcnd ] = Kfiltupd( task, opt, Sinv, Ainv, ...
                                            B, C, Qinv, X, Rinvy, Z, tol );
else  
    opt(lopt) = 0;  
    [ Sinv1, X1, Qin, E, rcnd ] = Kfiltupd( task, opt, Sinv, Ainv, ...
                                            B, C, Qinv, Rinv, X, Rinvy, Z, tol );
end
%
% end srif
