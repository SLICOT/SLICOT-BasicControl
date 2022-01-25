function [rcnd,sep,ferr,T,U] = carecond( A,Q,G,X,job,flag )
%CARECOND Computes an estimate of the reciprocal of the condition
%        number of a continuous-time algebraic Riccati equation.
% 
%        RCND = CARECOND(A,Q,G,X)  computes an estimate of the reciprocal
%        of the condition number of the Riccati equation
%
%           A'*X + X*A + Q - X*G*X = 0,                              (1)
%
%        where Q, G, and X are symmetric matrices.
%
%        [RCND,SEP,FERR] = CARECOND(A,Q,G,X,JOB,FLAG)  has additional
%        input and output arguments:
%
%        JOB specifies the calculation to be performed:
%        JOB = 1 : compute reciprocal condition number RCND and
%                  separation SEP;
%        JOB = 2 : compute the error bound FERR only; then, RCND
%                  and SEP must not be specified as output parameters;
%        JOB = 3 : compute reciprocal condition number, separation, and
%                  error bound.
%        Default:  JOB = 1.
%
%        FLAG is a vector containing options:
%        FLAG(1) specifies whether or not a real Schur form T of the
%        matrix Ac, Ac = A - G*X (or Ac = A - X*G, if FLAG(2) = 1),
%        is supplied on entry:
%        FLAG(1) = 0 :  a real Schur form of Ac should be computed.
%        FLAG(1) = 1 :  a real Schur form T of Ac is already available
%                       in A, and Q, G, and X contain the corresponding
%                       transformed matrices, i.e., U'*Q*U, U'*G*U, and 
%                       U'*X*U, respectively, where U is the orthogonal 
%                       matrix reducing Ac to the real Schur form T, 
%                       T = U'*Ac*U.
%        FLAG(2) specifies whether estimates for an equation (1) with
%        A replaced by A' should be computed:
%        FLAG(2) = 0 :  compute estimates for the equation (1).
%        FLAG(2) = 1 :  compute estimates for an equation (1) with
%                       A replaced by A'.
%        FLAG(3) specifies which part of the symmetric matrices Q and G
%        is to be used:
%        FLAG(3) = 0 :  upper triangular part;
%        FLAG(3) = 1 :  lower triangular part.
%        Default:       FLAG(1:3) = [0,0,0].
%
%        [RCND,SEP,FERR,T,U] = CARECOND(A,Q,G,X,3)  also returns the
%        computed real Schur form T of Ac, and the orthogonal
%        transformation matrix U. Similar commands may be used for 
%        JOB = 1 or JOB = 2. If FLAG(1) = 1, then T and U are not returned. 
%
%        If FLAG(1) = 1 or JOB = 1, CARECOND solves reduced Lyapunov
%        equations in the iterative estimation process, without updating
%        the right-hand sides and solutions. This scheme is very fast.
%        The corresponding Riccati equation is 
%              _   _     _   _ _ _
%           T'*X + X*T + Q + X*G*X = 0,
%
%        (or that with T replaced by T', if FLAG(2) = 1), where
%           _           _           _
%           X = U'*X*U, Q = U'*Q*U, G = U'*G*U. 
%
%        See also ARECOND
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, March 2003.
%        Revised: March 2009.
%

ni = nargin;
%
if ni < 4,
    error('Usage: RCND = CARECOND(A,Q,G,X,JOB,FLAG)')
elseif ni == 4,
    job  = 1;
    flag = [0;0;0];
elseif ni == 5,
    flag = [0;0;0];
elseif ni == 6,
    sz = length( flag );
    if sz < 3,  flag = flag(:);  flag(sz+1:3) = zeros( 3-sz,1 );  end 
end
%
eq = -2;  fact = flag(1);  reduced = fact;
if job == 1,
    reduced = 1;
    % Below, ferr returns T.
    if fact == 1,
        [rcnd,sep]      = arecond(eq,job,reduced,fact,A,Q,G,X,flag(2:3));
    else
        [rcnd,sep,ferr] = arecond(eq,job,reduced,fact,A,Q,G,X,flag(2:3));
    end
elseif job == 2,
    % Below, rcnd, sep, and ferr return ferr, T, and U, respectively.
    if fact == 1,
        rcnd            = arecond(eq,job,reduced,fact,A,Q,G,X,flag(2:3));
    else
        [rcnd,sep,ferr] = arecond(eq,job,reduced,fact,A,Q,G,X,flag(2:3));
    end
else
    if fact == 1,
        [rcnd,sep,ferr]     = arecond(eq,job,reduced,fact,A,Q,G,X,flag(2:3));
    else
        [rcnd,sep,ferr,T,U] = arecond(eq,job,reduced,fact,A,Q,G,X,flag(2:3));
    end
end
%
% end carecond
