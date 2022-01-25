% AREBENCH.F - Gateway function to generate the benchmark examples
%             for algebraic Riccati equations using SLICOT routines
%             BB01AD, for continuous-time case, and BB02AD, for
%             discrete-time case.
%
% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
%
% Matlab call:
%   [X,A,G,Q(,B,C)(,S)(,parval)] = arebench(dico,nr1,nr2 ...
%                                           (,flag,param,opt,filnam))
%
% Purpose:
%   To generate benchmark examples for continuous-time algebraic
%   Riccati equations (CARE)
%
%     0  =  Q + A'X + XA - XGX,
%
%   where A,G,Q,X are real N-by-N matrices, Q and G are symmetric and
%   may be given in factored form
%                   -1 T                         T
%      (I)   G = B R  B  ,           (II)   Q = C W C ,
%   where C is P-by-N, W P-by-P, B N-by-M, and R M-by-M, where W
%   and R are symmetric,
%
%   or
%
%   to generate benchmark examples for the discrete-time algebraic
%   Riccati equations (DARE)
%            T            T               T    -1  T       T
%     0  =  A X A - X - (A X B + S) (R + B X B)  (B X A + S ) + Q,
%
%   where the matrices Q and R are symmetric and Q may be given in
%   factored form (II).
%   If R is nonsingular and S = 0, the DARE can be rewritten as
%                  T             -1
%     0  =  X  -  A X (I_n + G X)  A  -  Q,
%
%   where I_n is the N-by-N identity matrix and G is factored as in (I).
%
% Input parameters:
%   dico   - integer option to indicate continuous- or discrete-time:
%            = 1 : continuous-time example.
%            = 2 : discrete-time example.
%   nr1    - the group of examples:
%            = 1 : parameter-free problems of fixed size.
%            = 2 : parameter-dependent problems of fixed size.
%            = 3 : parameter-free problems of scalable size.
%            = 4 : parameter-dependent problems of scalable size.
%   nr2    - the number of the example in group nr1.
%            Let NEXi be the number of examples in group i. Currently,
%            NEX1 =  6, NEX2 = 9, NEX3 = 2, NEX4 = 4, for
%                                                     continuous-time,
%            NEX1 = 13, NEX2 = 5, NEX3 = 0, NEX4 = 1, for discrete-time.
%            1 <= nr1 <= 4;
%            1 <= nr2 <= NEXi , where i = nr1.
%   flag   - (optional) vector containing options:
%               flag(1) = 1  : G is returned.
%               flag(1) = 0  : G is returned in factored form, i.e.,
%                              B and R from (I) are returned.
%               flag(2) = 1  : Q is returned.
%               flag(2) = 0  : Q is returned in factored form, i.e.,
%                              C and W from (II) are returned.
%            If dico = 1, it has size 2; if dico = 2, it has size 3:
%               flag(3) = 1  : The coefficient matrix S of the DARE
%                              is returned in array S.
%               flag(3) = 0  : The coefficient matrix S of the DARE
%                              is not returned.
%   param  - (optional) real vector of parameter values:
%            = epsilon, if dico = 1, nr1 = 2;
%            = [q,r], if dico = 1, for 4.1 (i.e., nr1 = 4, nr2 = 1);
%            = [a,b,c,beta1,beta2,gamma1,gamma2], if dico = 1, for 4.2;
%            = [mu,delta,kappa], if dico = 1, nr1 = 4, nr2 = 3;
%            = epsilon, if dico = 2, nr1 = 2, nr2 = 2:4;
%            = [tau,D,K,r], if dico = 2, nr1 = 2, nr2 = 5;
%            = R, if dico = 2, nr1 = 2, nr2 = 1 or nr1 = 4, nr2 = 1;
%   opt    - (optional) integer vector for additional
%            options/parameters:
%            if dico = 1,
%            = 1, (when nr1 = 2, nr2 = 9): generates the CARE for
%                  optimal state feedback (default);
%            = 2, (when nr1 = 2, nr2 = 9): generates the
%                  Kalman filter CARE.
%            = the number of vehicles, if nr1 = 3, nr2 = 1.
%            = the order of the matrix A, for 3.2, 4.1 or 4.2.
%            = the dimension of the second-order system, i.e., the
%              order of the stiffness matrix for 4.3 or 4.4.
%            if dico = 2,
%            = the order of the output matrix A, for 4.1.
%   filnam - (optional) if nr1 = nr2 = 4 and opt differs from the
%            default dimension of the second-order system (211), a
%            character string containing the name of the data file
%            to be used.
% Output parameters:
%   X      - real symmetric solution of Riccati equation (an exact
%            solution is available for the continuous-time examples
%            1.1, 1.2, 2.1, 2.3-2.6, 3.2, and for the discrete-time
%            examples 1.1, 1.3, 1.4, 2.1, 2.3-2.5, 4.1;
%            otherwise, X is an empty array).
%   A      - real n-by-n coefficient matrix.
%   G      - real symmetric n-by-n matrix.
%            If flag(1) = 1, then a non-factored coefficient matrix G
%            of the ARE is returned.
%            If flag(1) = 0, then array G contains the
%            'control weighting matrix' R from (I).
%   Q      - real symmetric n-by-n matrix.
%            If flag(2) = 1, then a non-factored coefficient matrix Q
%            of the ARE is returned.
%            If flag(2) = 0, then array Q contains the
%            'output weighting matrix' W from (II).
%   B      - real n-by-m input matrix.
%            If flag(1) = 0, this array contains the coefficient
%            matrix B of the ARE.
%   C      - real p-by-n output matrix.
%            If flag(2) = 0, this array contains the coefficient
%            matrix C of the ARE.
%   S      - if (dico = 2, flag(3) = 1), then this array contains the
%            coefficient matrix S of the DARE.
%   parval - the values used for the problem parameters:
%            either the values given in param, or default values,
%            if param is not used.

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   D. Sima, University of Bucharest, August 2001.
%
% Revisions:
%   V. Sima, Research Institute for Informatics, Nov. 2001, Oct. 2004,
%   Apr. 2009, Dec. 2012, May 2016, Jan. 2017.
%
