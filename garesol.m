% GARESOL.F - MEX-function for solving descriptor algebraic Riccati
%             equations using SLICOT routine SG02AD.
%
%       [X,F,ev,db,rcond1] = garesol(A,E,Q,R,B,L,flag)
%         [X,ev,db,rcond1] = garesol(A,E,Q,R,B,L,flag)
%         [X,ev,db,rcond1] = garesol(A,E,Q,G,flag)
%   [Z,scale,ev,db,rcond1] = garesol(A,E,Q,R,B,L,flag)
%   [Z,scale,ev,db,rcond1] = garesol(A,E,Q,G,flag)
%
%   GARESOL solves the continuous-time algebraic Riccati equations
%                                              -1
%      0 = Q + A'*X*E + E'*X*A - (L + E'*X*B)*R  *(L + E'*X*B)' ,  (1a)
%
%   or the discrete-time algebraic Riccati equations
%                                                 -1
%      E'*X*E = A'*X*A - (L + A'*X*B)*(R + B'*X*B)  (L + A'*X*B)' + Q;
%                                                                  (1b)
%   or, in alternate forms,
%
%                 -1
%      0 = Q - L*R  *L' + A1'*X*E + E'*X*A1 - E'*X*G*X*E,          (2a)
%
%   or
%                      -1                     -1
%      E'*X*E = Q - L*R  *L' + A1'*X*(I + G*X)  *A1,               (2b)
%
%                -1                 -1
%   where G = B*R  *B', A1 = A - B*R  *L'. The optimal feedback gain
%   matrix is
%
%           -1
%      F = R  (L+E'XB)' ,        for (1a),
%
%     and
%                  -1
%      F = (R+B'XB)  (L+A'XB)' , for (2a).
%
%   The function can also be used to compute an orthogonal basis of the
%   subspace which generates the solution.
%
%   The generalized Schur method is used on the extended pencils
%
%             [ A   0   B ]      [ E   0   0 ]
%           a [ Q   A'  L ] -  b [ 0  -E'  0 ] ,                   (3a)
%             [ L'  B'  R ]      [ 0   0   0 ]
%   or
%
%             [ A   0   B ]      [ E   0   0 ]
%           a [ Q  -E'  L ] -  b [ 0  -A'  0 ] ,                   (3b)
%             [ L'  0   R ]      [ 0  -B'  0 ]
%
%   to solve the equations (1a) or (1b).
%
%   Description of input parameters:
%   A      - real N-by-N system state matrix.
%   E      - real N-by-N descriptor system matrix.
%   Q      - normally, real symmetric N-by-N state weighting matrix.
%            If flag(5) <> 0, array Q stores the P-by-N factor C of Q,
%            Q = C'*C.
%   R      - normally, real symmetric M-by-M input weighting matrix.
%            If R is an input parameter and flag(6) <> 0, array R stores
%            the P-by-M factor D of R, R = D'*D.
%   B      - real N-by-M input matrix.
%   L      - real N-by-M coupling matrix.
%   G      - real symmetric N-by-N matrix.
%   flag   - (optional) vector of length 10 containing options:
%               flag(1)  = 0 : solve the continuous-time equation (1a);
%                              otherwise, solve the discrete-time
%                              equation (1b).
%               flag(2)  = 0 : compute the stabilizing solution;
%                              otherwise, compute the anti-stabilizing
%                              solution.
%               flag(3)  = 0 : compute both the solution X and the
%                              feedback gain matrix F;
%                        > 0 : compute the solution only;
%                        < 0 : compute the related subspace.
%               flag(4)      : tolerance to check the singularity of
%                              the matrix pencil. If tol <= 0, tol will
%                              be replaced by EPS, the machine epsilon.
%               flag(5)  = 0 : matrix Q is not factored; otherwise Q is
%                              assumed factored as Q = C'*C, where C is
%                              P-by-N.
%               flag(6)  = 0 : matrix R is not factored; otherwise R is
%                              assumed factored as R = D'*D, where D is
%                              P-by-M.
%               flag(7)  = 0 : the upper triangles of matrices Q and G,
%                              or Q and R, are stored; otherwise, lower
%                              triangles are stored.
%               flag(8)  = 0 : matrix L is zero; otherwise, L is given.
%               flag(9)  = 0 : use a scaling strategy (for given R, B);
%                              otherwise, do not use a scaling strategy.
%               flag(10) = 0 : iterative refinement should not be used
%                              for solving the system of algebraic
%                              equations giving the solution matrix X;
%                              otherwise, use iterative refinement.
%            Default:    flag(1:4) = [0,0,1,0,0,0,0,1,0,0].
%
%   Description of output parameters:
%   X      - N-by-N real symmetric solution of Riccati equations (1)
%            or (2).
%   F      - M-by-N real state feedback matrix, returned only when R
%            and B are input arguments, the number of output parameters 
%            is at least 2, and X is computed.
%   Z      - 2N-by-N real matrix with orthonormal columns, the basis of
%            the subspace related to the solution X.
%   scale  - the scaling factor used internally, which should multiply 
%            the submatrix Z(N+1:2N,:) to recover X from Z.
%   ev     - complex vector of length N.
%            the "a" part of the associated eigenvalues in (3).
%   db     - real vector of length N containing the "b" part of
%            the associated eigenvalues in (3), i.e., the associated
%            eigenvalues with respect to X are ev(k)/db(k), k=1,...,N.
%   rcond1 - estimate of the reciprocal of the 1-norm condition number
%            of the N-th order system of algebraic equations from which
%            the solution matrix X is obtained.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Katholieke Univ. Leuven, Belgium, March 2003.
%
% Revisions:
%   -
%
