% ARESOLC.F - MEX-function for solving algebraic Riccati equations
%             using SLICOT routines SB02RD, SB02ND, SB02MT and SB02OD.
%             SB02RD is a more accurate version of SB02MD, also 
%             computing condition numbers for Riccati equations.
%
%   [X(,F),ev(,db),rcond1(,acc)] = aresolc(method,A,Q,R,B,L,flag)
%       [X,ev(,db),rcond1(,acc)] = aresolc(method,A,Q,G,flag)
%       [Z,scale,ev(,db),rcond1] = aresolc(method,A,Q,R,B,L,flag)
%       [Z,scale,ev(,db),rcond1] = aresolc(method,A,Q,G,flag)
%
%   method = 1:        [X,F,ev,rcond1] = aresolc(1,A,Q,R,B,L,flag)
%                    [X,ev,rcond1,acc] = aresolc(1,A,Q,R,B,L,flag)
%                    [X,ev,rcond1,acc] = aresolc(1,A,Q,G,flag)
%                  [Z,scale,ev,rcond1] = aresolc(1,A,Q,R,B,L,flag)
%                  [Z,scale,ev,rcond1] = aresolc(1,A,Q,G,flag)
%   method = 2:     [X,F,ev,db,rcond1] = aresolc(2,A,Q,R,B,L,flag)
%                     [X,ev,db,rcond1] = aresolc(2,A,Q,R,B,L,flag)
%                     [X,ev,db,rcond1] = aresolc(2,A,Q,G,flag)
%               [Z,scale,ev,db,rcond1] = aresolc(2,A,Q,R,B,L,flag)
%               [Z,scale,ev,db,rcond1] = aresolc(2,A,Q,G,flag)
%   method = 3:     [X,F,ev,db,rcond1] = aresolc(3,A,Q,R,B,L,flag)
%                     [X,ev,db,rcond1] = aresolc(3,A,Q,R,B,L,flag)
%               [Z,scale,ev,db,rcond1] = aresolc(3,A,Q,R,B,L,flag)
%         
%   ARESOLC solves the continuous-time algebraic Riccati equations
%                                      -1
%      0 = Q + A'*X + X*A - (L + X*B)*R  *(L + X*B)' ,             (1a)
%
%   or the discrete-time algebraic Riccati equations
%                                            -1
%      X = A'*X*A - (L + A'*X*B)*(R + B'*X*B)  (A'*X*B + L)' + Q;  (1b)
%
%   or, in alternate forms,
%
%      0 = Q + A'*X + X*A - X*G*X,                                 (2a)
%
%   or
%                            -1
%      X = Q + A'*X*(I + G*X)  *A.                                 (2b)
%
%   The function can also be used to compute an orthogonal basis of the
%   subspace which generates the solution.
%
%   When method = 1, the Schur method is used on Hamiltonian matrices
%                _      _
%             [  A     -G  ]         [ A  -G ]
%             [  _      _  ]   or    [       ]                     (3a)
%             [ -Q     -A' ]         [-Q  -A']
%
%   or symplectic matrices
%
%        [  _ -1       _ -1  _    ]       [   -1         -1     ]
%        [  A          A   * G    ]       [  A          A  *G   ]
%        [ _  _-1   _    _ _ -1 _ ]   or  [    -1          -1   ]  (3b)
%        [ Q* A     A' + Q*A   *G ]       [ Q*A    A' + Q*A  *G ]
%
%         _          -1     _          -1     _      -1
%   where A = A - B*R  *L', Q = Q - L*R  *L', G = B*R  *B'.
%   
%   When method = 2, the generalized Schur method is used for discrete-
%   time problems, on the pencil
%               _                  _
%             [ A   0 ]      [ I   G  ]
%           a [ _     ]  - b [     _  ],                           (4a)
%             [-Q   I ]      [ 0   A' ]
%         _  _     _
%   where A, Q and G are defined above, or on the pencil
%
%             [ A   0 ]      [ I   G  ]
%           a [       ]  - b [        ],                           (4b)
%             [-Q   I ]      [ 0   A' ]
%
%   for (1b) or (2b).
% 
%   When method = 3, the generalized Schur method is used on the
%   extended pencils
%
%             [ A  0  B ]      [ I   0  0 ]
%           a [ Q  A' L ] -  b [ 0  -I  0 ],                       (5a)
%             [ L' B' R ]      [ 0   0  0 ]     
%   or
%
%             [ A  0  B ]      [ I  0  0 ]
%           a [ Q -I  L ] -  b [ 0 -A' 0 ],                        (5b)
%             [ L' 0  R ]      [ 0 -B' 0 ]
%
%   to solve the equations (1a) or (1b).
%
%   Description of input parameters:
%   method - integer option to indicate the method to be used:
%            = 1 : use ordinary Schur method on matrix (3a) or (3b).
%            = 2 : use generalized Schur method on pencil (4a) or (4b).
%            = 3 : use generalized Schur method on pencil (5a) or (5b).
%   A      - real N-by-N system state matrix.
%   Q      - real symmetric N-by-N state weighting matrix.
%   R      - real symmetric M-by-M input weighting matrix.
%            When method = 1, or 2, R must be nonsingular.
%   B      - real N-by-M input matrix.
%   L      - real N-by-M coupling matrix.
%   G      - real symmetric N-by-N matrix.
%   flag   - (optional) vector containing options and parameters.
%            method = 1 : flag has length 4
%                flag(1) = 0 : solve the continuous-time equation
%                              (1a) or (2a); otherwise, solve the
%                              discrete-time equation (1b) or (2b).
%                flag(2) = 0 : compute the stabilizing solution;
%                              otherwise, compute the anti-stabilizing
%                              solution.
%                flag(3) = 0 : use a scaling strategy;
%                              otherwise, do not use a scaling strategy.
%                flag(4) = 0 : compute both the solution X and the
%                              feedback gain matrix F;
%                        = 1 : compute the solution X only;
%                        = 2 : compute the solution X and the accuracy
%                              estimates (separation, reciprocal 
%                              condition number and forward error bound);
%                        < 0 : compute the invariant subspace.
%            Default:    flag(1:4) = [0,0,0,1].
%
%            method = 2 : flag has length 3
%                flag(1) = 0 : compute the stabilizing solution;
%                              otherwise, compute the anti-stabilizing
%                              solution.
%                flag(2) = 0 : compute both the solution X and the
%                              feedback gain matrix F;
%                        > 0 : compute the solution X only;
%                        < 0 : compute the deflating subspace.
%                flag(3)     : tolerance to check the singularity of
%                              the matrix pencil. If tol <= 0, tol will
%                              be replaced by EPS, the machine epsilon.
%            Default:    flag(1:3) = [0,1,0].
%
%            method = 3 : flag has length 4
%                flag(1) = 0 : solve the continuous-time equation (1a);
%                              otherwise, solve the discrete-time
%                              equation (1b).
%                flag(2) = 0 : compute the stabilizing solution;
%                              otherwise, compute the anti-stabilizing
%                              solution.
%                flag(3) = 0 : compute both the solution X and the
%                              feedback gain matrix F;
%                        > 0 : compute the solution only;
%                        < 0 : compute the related subspace.
%                flag(4)     : tolerance to check the singularity of
%                              the matrix pencil. If tol <= 0, tol will
%                              be replaced by EPS, the machine epsilon.
%            Default:    flag(1:4) = [0,0,1,0].
%
%   Description of output parameters:
%   X      - N-by-N real symmetric solution of Riccati equations (1)
%            or (2).
%   F      - M-by-N real state feedback matrix, returned only when R
%            and B are input arguments, NLHS >= 2 and X is computed.
%   Z      - 2N-by-N real matrix with orthonormal columns, the basis of
%            the subspace related to the solution X.
%   scale  - the scaling factor used (internally), which should multiply
%            the submatrix Z(N+1:2N,:) to recover X from Z.
%   ev     - complex vector of length N.
%            method = 1   : the associated eigenvalues.
%            method = 2,3 : the "a" part of the associated eigenvalues
%                           (see (4) or (5)).
%   db     - real vector of length N containing the "b" part of
%            the associated eigenvalues (see (4) or (5)), i.e., the
%            associated eigenvalues with respect to X are ev(k)/db(k), 
%            k=1,...,N.
%   rcond1 - estimate of the reciprocal of the 1-norm condition number
%            of the N-th order system of algebraic equations from which
%            the solution matrix X is obtained.
%   acc    - real vector of length 3 containing estimates of the 
%            separation, reciprocal condition number and forward error 
%            bound, repectively. These results can be obtained only if
%            meth = 1 and flag(4) = 2.
%
%   Comments:
%   For METH = 2 or 3, this function behaves as ARESOL. 
% 
%   See also SLCAREGS, SLCARESC, SLDAREGS, SLDAREGSV, SLDARESC

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1999.
%
% Revisions:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2002.
%
