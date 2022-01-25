%GENLEQ  MEX-function for solving generalized linear matrix equations
%        using SLICOT routines SB04OD, SG03AD, and SG03BD.
%
%   [X,(Y,)dif] = genleq(task,A(,D,B),E,C,(F,)flag,trans)
%
%   [X,Y,dif] = genleq(1,A,D,B,E,C,F,flag,trans)
%       [X,Y] = genleq(1,A,D,B,E,C,F,flag,trans)
%     [X,dif] = genleq(2,A,D,B,E,C,flag,trans)
%         [X] = genleq(2,A,D,B,E,C,flag,trans)
%     [X,sep] = genleq(3,A,E,C,flag,trans)
%         [X] = genleq(3,A,E,C,flag,trans)
%         [X] = genleq(4,A,E,C,flag,trans)
%         
%   GENLEQ solves generalized Sylvester and Lyapunov linear matrix
%   equations:
%
%       (  A*X - Y*B = C,
%       <                                                          (1.a)
%       (  D*X - Y*E = F,
%
%       (  A'*X + D'*Y = C,
%       <                                                          (1.b)
%       (  X*B' + Y*E' = -F,
%
%       A*X*B - D*X*E = C,                                           (2)
%
%       op(A)'*X*op(E) + op(E)'*X*op(A) = C,                       (3.a)
%
%       op(A)'*X*op(A) - op(E)'*X*op(E) = C,                       (3.b)
%
%       op(A)'*op(X)'*op(X)*op(E) + op(E)'*op(X)'*op(X)*op(A) =
%                                               -op(C)'*op(C),     (4.a)
%
%       op(A)'*op(X)'*op(X)*op(A) - op(E)'*op(X)'*op(X)*op(E) =
%                                               -op(C)'*op(C),     (4.b)
%
%   where op(M) = M, if trans = 0, and op(M) = M', if trans <> 0, and
%   X in the equations (4) is upper triangular.
%
%   Description of other input parameters:
%   task  - integer option to determine the equation type
%           = 1 : solve the linear matrix equation pairs (1.a) or (1.b).
%           = 2 : solve the linear matrix equation (2).
%           = 3 : solve the linear matrix equation (3.a) or (3.b).
%           = 4 : solve the linear matrix equation (4.a) or (4.b).
%   flag  - (optional) integer vector of length 2 containing options
%           task = 1, 2 :
%                flag(1) = 1 : (A,D) is in generalized Schur form;
%                              otherwise, (A,D) is in general form.
%                flag(2) = 1 : (B,E) is in generalized Schur form;
%                              otherwise, (B,E) is in general form.
%           task = 3, 4 :
%                flag(1) = 0 : solve the continuous-time equation (x.a);
%                              otherwise, solve the discrete-time
%                              equation (x.b), where x is 3 or 4.
%                flag(2) = 1 : (A,E) is in generalized Schur form;
%                              otherwise, (A,E) is in general form.
%           Default:      flag(1) = 0, flag(2) = 0.
%   trans - (optional) integer specifying a transposition option
%           trans = 0 : solve the equations (1.a), (2); or solve (3)
%                       or (4) with op(M) = M.
%                       otherwise, solve the "transposed" equations
%                       (1.b); or solve (3) or (4) with op(M) = M'.
%           Default:      trans = 0.
%
%   Description of other output parameters:
%   dif   - (optional) estimator of Dif[(A,D),(B,E)].
%           dif is not computed for (1.b), i.e., for task = 1 with
%           trans <> 0.
%   sep   - (optional) estimator of sep(A,E).
%
%   Comments
%   1. Currently there is no available SLICOT routine for equations (2).
%   2. For task = 4, the pencil (A,E) must be stable, i.e., all
%      eigenvalues must have negative real parts, for (4.a), or moduli
%      strictly less than 1, for (4.b).
%
%   See also SLGELY, SLGESG, SLGEST, SLGSLY, SLGSST
%

%   RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   H. Xu, TU Chemnitz, FR Germany, Dec. 1998.
%
%   Revisions:
%   V. Sima, Katholieke Univ. Leuven, Belgium, May 1999.
