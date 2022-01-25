%LINMEQ  MEX-function for solving Sylvester and Lyapunov equations
%        using SLICOT routines SB04MD, SB04ND, SB04PD, SB04QD, 
%        SB04RD, SB03MD, and SB03OD.
%
%   [X(,sep)] = linmeq(task,A,(B,)C,flag,trans(,schur))
%
%       [X] = linmeq(1,A,B,C,flag,trans,schur)
%   [X,sep] = linmeq(2,A,C,flag,trans)
%       [X] = linmeq(2,A,C,flag,trans)
%       [X] = linmeq(3,A,C,flag,trans)
%         
%   LINMEQ solves various Sylvester and Lyapunov matrix equations:
%
%        op(A)*X + X*op(B) = C,                           (1a)
%
%        op(A)*X*op(B) + X = C,                           (1b)
%
%        op(A)'*X + X*op(A) = C,                          (2a)
%
%        op(A)'*X*op(A) - X = C,                          (2b)
%
%        op(A)'*(op(X)'*op(X)) + (op(X)'*op(X))*op(A) =
%                              -  op(C)'*op(C),           (3a)
%
%        op(A)'*(op(X)'*op(X))*op(A) - op(X)'*op(X) =
%                                    - op(C)'*op(C),      (3b)
%
%   where op(M) = M, or M'.
%
%   Description of other input parameters:
%   task  - integer option to determine the equation type:
%           = 1 : solve the Sylvester equation (1a) or (1b);
%           = 2 : solve the Lyapunov equation (2a) or (2b);
%           = 3 : solve for the Cholesky factor op(X) the Lyapunov
%                 equation (3a) or (3b).
%   flag  - (optional) integer vector of length 3 or 2 containing
%           options.
%           task = 1 : flag has length 3
%                flag(1) = 0 : solve the continuous-time equation (1a);
%                              otherwise, solve the discrete-time
%                              equation (1b).
%                flag(2) = 1 : A is (quasi) upper triangular;
%                          2 : A is upper Hessenberg;
%                              otherwise, A is in general form.
%                flag(3) = 1 : B is (quasi) upper triangular;
%                          2 : B is upper Hessenberg;
%                              otherwise, B is in general form.
%           task = 2 : flag has length 2
%                flag(1) = 0 : solve continuous-time equation (2a);
%                              otherwise, solve discrete-time
%                              equation (2b).
%                flag(2) = 1 : A is (quasi) upper triangular;
%                              otherwise, A is in general form.
%           task = 3 : flag has length 2
%                flag(1) = 0 : solve continuous-time equation (3a);
%                              otherwise, solve discrete-time
%                              equation (3b).
%                flag(2) = 1 : A is (quasi) upper triangular;
%                              otherwise, A is in general form.
%           Default:    flag(1) = 0, flag(2) = 0 (, flag(3) = 0).
%   trans - (optional) integer specifying a transposition option.
%           trans = 0 : solve the equations (1) - (3) with op(M) = M.
%                   1 : solve the equations (1) - (3) with op(M) = M'.
%                   2 : solve the equations (1) with op(A) = A';  op(B) = B;
%                   3 : solve the equations (1) with op(A) = A;   op(B) = B'.
%           Default:    trans = 0.
%   schur - (optional) integer specifying whether the Hessenberg-Schur
%           or Schur method should be used.
%           Available for task = 1. 
%           schur = 1 : Hessenberg-Schur method (one matrix is reduced
%                       to Schur form).
%           schur = 2 : Schur method (two matrices are reduced to Schur
%                       form).
%           Default:    schur = 1.
%
%   Description of other output parameters:
%   sep   - (optional) estimator of Sep(op(A),-op(A)') for (2.a) or
%           Sepd(A,A') for (2.b).
%
%   Comments
%   1. For equation (1a) or (1b), when schur = 1, the Hessenberg-Schur 
%      method is used, reducing one matrix to Hessenberg form and the
%      other one to a real Schur form.
%      Otherwise, both matrices are reduced to real Schur forms.
%      If one or both matrices are already reduced to Schur/Hessenberg
%      forms, this could be specified by flag(2) and flag(3).
%      For general matrices, the Hessenberg-Schur method could be 
%      significantly more efficient than the Schur method.
%   2. For equation (2a) or (2b), matrix C is assumed symmetric.
%   3. For equation (3a) or (3b), matrix A must be stable or convergent,
%      respectively.
%   4. For equation (3a) or (3b), the computed matrix X is the Cholesky
%      factor of the solution, i.e., the real solution is op(X)'*op(X),
%      where X is an upper triangular matrix.
%
% See also SLDSYL, SLLYAP, SLSTEI, SLSTLY, SLSTST, SLSYLV
%

%   RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   H. Xu, TU Chemnitz, FR Germany, Dec. 1998.
%
%   Revisions:
%   V. Sima, Katholieke Univ. Leuven, Belgium, May 1999, May, Sep. 2000.
%   D. Sima, University of Bucharest, Romania, May 2000.
