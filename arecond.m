% ARECOND.F - MEX-function for estimating the conditioning of
%             Lyapunov and algebraic Riccati equations using SLICOT
%             routines SB03QD, SB03SD, SB02QD, and SB02SD.
%
%   If fact = 0,
%   [(rcond(,sep))(,ferr)(,T(,U))] =
%            arecond(eq,job,reduced,fact,A,C,X(,flag))       |eq| = 1;
%   [(rcond(,sep))(,ferr)(,T(,U))] =
%            arecond(eq,job,reduced,fact,A,Q,G,X(,flag))     |eq| = 2;
%   If fact = 1 and reduced = 1,
%   [(rcond(,sep))(,ferr)] =
%            arecond(eq,job,reduced,fact,T,C,X(,flag))       |eq| = 1;
%   [(rcond(,sep))(,ferr)] =
%            arecond(eq,job,reduced,fact,T,Q,G,X(,flag))     |eq| = 2;
%   If fact = 1 and reduced = 0,
%   [(rcond(,sep))(,ferr)] =
%            arecond(eq,job,reduced,fact,A,T,U,C,X(,flag))   |eq| = 1;
%   [(rcond(,sep))(,ferr)] =
%            arecond(eq,job,reduced,fact,A,T,U,Q,G,X(,flag)) |eq| = 2.
%
%   1. ARECOND estimates the conditioning and compute an error bound on
%   the solution of the real continuous-time Lyapunov matrix equation
%
%         op(A)'*X + X*op(A) = C,                                    (1)
%
%   or of the real discrete-time Lyapunov matrix equation
%
%         op(A)'*X*op(A) - X = C,                                    (2)
%
%   where op(A) = A or A' (A**T) and C is symmetric (C = C').
%
%   2. ARECOND estimates the conditioning and compute an error bound on
%   the solution of the real continuous-time matrix algebraic Riccati
%   equation
%
%         op(A)'*X + X*op(A) + Q - X*G*X = 0,                        (3)
%
%   or of the real discrete-time matrix algebraic Riccati equation
%                                 -1
%         X = op(A)'*X*(I_n + G*X)  *op(A) + Q,                      (4)
%
%   where Q and G are symmetric. Let Ac denote the closed-loop matrix,
%   hence, in the continuous-time case,
%
%         Ac = A - G*X,          if op(A) = A,  or
%         Ac = A - X*G,          if op(A) = A',
%
%   and, in the discrete-time case, 
%
%        Ac = inv(I_n + G*X)*A,  if op(A) = A,  or
%        Ac = A*inv(I_n + X*G),  if op(A) = A'.
%
%   Riccati equations (3) and (4) are equivalent to the following
%   reduced equations, respectively,
%               _   _         _   _ _ _
%        op(T)'*X + X*op(T) + Q + X*G*X = 0,
%        _          _                _ _ _         _
%        X = op(T)'*X*op(T) + op(T)'*X*G*X*op(T) + Q,
%         _           _               _
%   where X = U'*X*U, Q = U'*Q*U, and G = U'*G*U, with U the
%   orthogonal matrix reducing Ac to a real Schur form, T = U'*Ac*U.
%   Similar, simpler formulas stand for Lyapunov equations (1) and (2).
%
%   Description of input parameters:
%   eq     - integer option to indicate the equation type:
%            = -1 : continuous-time Lyapunov equation (1);
%            =  1 : discrete-time Lyapunov equation (2);
%            = -2 : continuous-time Riccati equation (3);
%            =  2 : discrete-time Riccati equation (4).
%   job    - integer option to indicate the calculation to be performed:
%            =  1 : compute the reciprocal condition number and
%                   separation;
%            =  2 : compute the error bound only;
%            =  3 : compute the reciprocal condition number, the
%                   separation, and the error bound.
%   reduced- integer option specifying whether or not the original
%            Lyapunov equations should be solved in the iterative
%            estimation process (also for (3) and (4)), as follows:
%            =  0 :  solve the original Lyapunov equations, updating
%                    the right-hand sides and solutions with the
%                    matrix U, e.g., RHS <-- U'*RHS*U;
%            =  1 :  solve reduced Lyapunov equations only, without
%                    updating the right-hand sides and solutions.
%                    This scheme is faster, but sometimes slightly
%                    less accurate.
%   fact   - integer option specifying whether or not the real Schur
%            factorization of the matrix A or Ac is supplied on entry,
%            as follows:
%            =  0 :  the Schur factorization of A (if |eq| = 1) or Ac
%                    (if |eq| = 2) will be computed and the factors can
%                    be stored in the output parameters T and U (if
%                    reduced = 0);
%            =  1 :  the input parameters T and U (if reduced = 0)
%                    contain the factors from the real Schur
%                    factorization of the matrix A (if |eq| = 1) or Ac
%                    (if |eq| = 2).
%   A      - if fact = 0 or reduced = 0, the real n-by-n system state
%            matrix A.
%   T      - if fact = 1, a real n-by-n Schur form of A, if |eq| = 1,
%            or Ac, if |eq| = 2.
%   U        If fact = 1 and reduced = 0, the n-by-n orthogonal matrix
%            U which reduced the matrix A or Ac to a real Schur form T.
%   C      - real symmetric n-by-n right-hand side matrix of the
%            original equation, if reduced = 0, or of the reduced
%            equation (with matrix T), if reduced = 1.
%   X      - real symmetric n-by-n solution matrix of original equation,
%            if reduced = 0, or of the reduced equation (with matrix T),
%            if reduced = 1.
%   Q      - real symmetric n-by-n state weighting matrix of the Riccati
%            equation. Matrix Q must correspond to that of the reduced
%            equation, if reduced = 1.
%   G      - real symmetric n-by-n matrix of the Riccati equation.
%            Matrix G must correspond to that of the reduced equation,
%            if reduced = 1.
%   flag   - (optional) vector containing options.
%            flag(1) specifies the form of op(A) to be used, as follows:
%                    = 0 : op(A) = A;
%                    = 1 : op(A) = A'.
%            flag(2) specifies which part of the symmetric matrix C
%                    (if |eq| = 1), or the symmetric matrices Q and G
%                    (if |eq| = 2), is to be used, as follows:
%                    = 0 : upper triangular part;
%                    = 1 : lower triangular part.
%            Default:      flag(1:2) = [0,0].
%
%   Description of output parameters:
%   rcond  - (optional) if job = 1 or job = 3, an estimate of the
%            reciprocal condition number of the Lyapunov or Riccati
%            equation.
%   sep    - (optional) if job = 1 or job = 3, an estimate of the
%            separation, i.e., of the quantity
%               sep(op(A),-op(A)'),   if eq = -1;
%               sepd(op(A),op(A)'),   if eq =  1;
%               sep(op(Ac),-op(Ac)'), if eq = -2;
%               sepd(op(Ac),op(Ac)'), if eq =  2.
%   ferr   - (optional) if job = 2 or job = 3, an estimated forward
%            error bound for the solution X. If Xtrue is the true
%            solution, ferr bounds the magnitude of the largest entry in
%            (X - Xtrue) divided by the magnitude of the largest entry
%            in X.
%   T      - (optional) if fact = 0, the n-by-n real Schur form matrix
%            corresponding to A or Ac.
%   U      - (optional) if fact = 0 and reduced = 0, the n-by-n real
%            orthogonal matrix which reduced the matrix A or Ac to the
%            real Schur form T.

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Katholieke Univ. Leuven, Belgium, March 2003.
%
% Revisions:
%   -
%
