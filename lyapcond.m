function [rcnd,sep,ferr,T,U] = lyapcond( A,C,X,job,flag )
%LYAPCOND Computes an estimate of the reciprocal of the condition
%        number of a continuous-time Lyapunov equation.
% 
%        RCND = LYAPCOND(A,C,X)  computes an estimate of the reciprocal
%        of the condition number of the Lyapunov equation
%
%           A'*X + X*A = C,                                        (1)
%
%        where C and X are symmetric matrices.
%
%        [RCND,SEP,FERR] = LYAPCOND(A,C,X,JOB,FLAG)  has additional
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
%        matrix A is supplied on entry:
%        FLAG(1) = 0 :  a real Schur form of A should be computed.
%        FLAG(1) = 1 :  a real Schur form T of A is already available
%                       in A, and C and X contain the corresponding
%                       transformed matrices, i.e., U'*C*U and U'*X*U,
%                       respectively, where U is the orthogonal matrix
%                       reducing A to the real Schur form T, T = U'*A*U.
%        FLAG(2) specifies whether estimates for an equation (1) with
%        A replaced by A' should be computed:
%        FLAG(2) = 0 :  compute estimates for the equation (1).
%        FLAG(2) = 1 :  compute estimates for an equation (1) with
%                       A replaced by A'.
%        FLAG(3) specifies which part of the symmetric matrix C is to
%        be used:
%        FLAG(3) = 0 :  upper triangular part;
%        FLAG(3) = 1 :  lower triangular part.
%        Default:       FLAG(1:3) = [0,0,0].
%
%        [RCND,SEP,FERR,T,U] = LYAPCOND(A,C,X,3)  also returns the
%        computed real Schur form T of A, and the orthogonal
%        transformation matrix U. Similar commands may be used for 
%        JOB = 1 or JOB = 2. If FLAG(1) = 1, then T and U are not returned. 
%
%        If FLAG(1) = 1 or JOB = 1, LYAPCOND solves reduced Lyapunov
%        equations in the iterative estimation process, without updating
%        the right-hand sides and solutions. This scheme is very fast.
%        The corresponding Lyapunov equation is 
%              _   _     _
%           T'*X + X*T = C,
%
%        (or that with T replaced by T', if FLAG(2) = 1), where
%           _           _
%           X = U'*X*U, C = U'*C*U. 
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
if ni < 3,
    error('Usage: RCND = LYAPCOND(A,C,X,JOB,FLAG)')
elseif ni == 3,
    job  = 1;
    flag = [0;0;0];
elseif ni == 4,
    flag = [0;0;0];
elseif ni == 5,
    sz = length( flag );
    if sz < 3,  flag = flag(:);  flag(sz+1:3) = zeros( 3-sz,1 );  end 
end
%
eq = -1;  fact = flag(1);  reduced = fact;
if job == 1,
    reduced = 1;
    % Below, ferr returns T.
    if fact == 1,
        [rcnd,sep]      = arecond(eq,job,reduced,fact,A,C,X,flag(2:3));
    else
        [rcnd,sep,ferr] = arecond(eq,job,reduced,fact,A,C,X,flag(2:3));
    end
elseif job == 2,
    % Below, rcnd, sep, and ferr return ferr, T, and U, respectively.
    if fact == 1,
        rcnd            = arecond(eq,job,reduced,fact,A,C,X,flag(2:3));
    else
        [rcnd,sep,ferr] = arecond(eq,job,reduced,fact,A,C,X,flag(2:3));
    end
else
    if fact == 1,
        [rcnd,sep,ferr]     = arecond(eq,job,reduced,fact,A,C,X,flag(2:3));
    else
        [rcnd,sep,ferr,T,U] = arecond(eq,job,reduced,fact,A,C,X,flag(2:3));
    end
end
%
% end lyapcond
