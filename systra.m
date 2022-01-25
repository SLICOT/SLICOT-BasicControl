%SYSTRA  MEX-function for computing system transformations with scaling,
%        block diagonal decomposition, or Schur form reduction of the
%        state matrix, using SLICOT routines TB01ID, TB01KD, TB01LD,
%        and TB01WD.
%
%   [A,B,C(,num,ev,U)] = systra(task,A,B,C(,job,par))
%   [A,B,C(,num,ev,U)] = systra(task,A,B,C(,flag,par))
%
%          [A,B,C,U] = systra(1,A,B,C,job,par)
%       [A,B,C,ev,U] = systra(2,A,B,C)
%   [A,B,C,num,ev,U] = systra(3,A,B,C,flag,par)
%   [A,B,C,num,ev,U] = systra(4,A,B,C,flag,par)
%
%   SYSTRA applies a specified similarity transformation to a given
%                     [ A  B ]
%   system matrix S = [      ]. The possible transformations are:
%                     [ C  0 ]
%
%                                                        [ D  0 ]
%   1. Do balance on the matrix S with a diagonal matrix [      ];
%                                                        [ 0  I ]
%
%   2. Reduce A to Schur form using an orthogonal matrix U, and
%      transform B and C accordingly,
%
%                                    [* . . . * |* . *]--
%                                    [  .     . |.   .] |
%      [ U'| 0 ][ A | B ][ U | 0 ]   [    .   . |.   .] N
%      [---|---][---|---][---|---] = [      . . |.   .] |  ;       (1)
%      [ 0 | I ][ C | 0 ][ 0 | I ]   [        * |* . *]--
%                                    [----------|-----]
%                                    [* . . . * |     ]-|-
%                                    [.       . |  0  ] P
%                                    [* . . . * |     ]-|-
%                                     |-- N --|  |-M-|
%
%   3. Reduce A to real Schur form with eigenvalues in required order,
%      using an orthogonal matrix U, and transform B and C accordingly,
%
%                                   |--[* . * * . * |* . *]--
%                                  num [  . . .   . |* . *] |
%                                   |--[    * * . * |.   .] |
%      [ U'| 0 ][ A | B ][ U | 0 ]     [      * . * |.   .] N
%      [---|---][---|---][---|---] =   [        . . |.   .] |  ;   (2)
%      [ 0 | I ][ C | 0 ][ 0 | I ]     [          * |* . *]--
%                                      [------------|-----]
%                                      [* . . . . * |     ]-|-
%                                      [.         . |  0  ] P
%                                      [* . . . . * |     ]-|-
%                                       |--- N ---|  |-M-|
%   4. Reduce A to a block diagonal form with eigenvalues in required
%      order, using a nonsingular matrix U, and transform B and C
%      accordingly,
%
%                                     |--[* . *       |* . *]--
%                                    num [  . .       |* . *] |
%         -1                          |--[    *       |.   .] |
%      [ U  | 0 ][ A | B ][ U | 0 ]      [      * . * |.   .] N
%      [----|---][---|---][---|---] =    [        . . |.   .] | .  (3)
%      [ 0  | I ][ C | 0 ][ 0 | I ]      [          * |* . *]--
%                                        [------------|-----]
%                                        [* . . . . * |     ]-|-
%                                        [.         . |  0  ] P
%                                        [* . . . . * |     ]-|-
%                                         |--- N ---|  |-M-|
%
%   Description of other input parameters:
%   task  - integer option to indicate which transformation is needed:
%           = 1 : do balancing;
%           = 2 : compute the transformation (1);
%           = 3 : compute the transformation (2);
%           = 4 : compute the transformation (3).
%   job   - (for task = 1 only: optional) integer option parameter:
%           job = 1  : balancing involves the matrix A only.
%           job = 2  : balancing involves the matrices A and B.
%           job = 3  : balancing involves the matrices A and C.
%           otherwise, balancing involves all matrices A, B, and C.
%           default:  job = 0;
%   flag  - (optional) integer vector containing options.
%           task = 1, 2 : flag is not used.
%           task = 3, 4 : flag is a vector of length 2
%              flag(1) = 0 : the system is continuous-time (default),
%                            otherwise, the system is discrete-time.
%              flag(2) = 0 : reorder the stable eigenvalues of A on the
%                            top left diagonal block (default),
%                            otherwise, reorder the unstable eigenvalues
%                            on the top left diagonal block.
%   par   - (optional) real parameter specifying the following values:
%           task = 1    : the maximum allowed reduction in the 1-norm of
%                         S if zero rows or columns are encountered;
%                         default:  par = 10;
%           task = 2    : not used.
%           task = 3, 4 : the stability or instability boundary for the
%                         eigenvalues of A of interest. For the
%                         discrete-time case, par >= 0 represents the
%                         boundary value for the moduli of eigenvalues.
%             default:    -sqrt(epsilon_machine) for continuous-time;
%                      1.0-sqrt(epsilon_machine) for discrete-time.
%
%   Description of other output parameters:
%   num   - (optional) integer counting for the number of the required
%           eigenvalues, or, equivalently, the size of the diagonal
%           block containing the required eigenvalues. If task = 1, 2,
%           num is not available.
%   ev    - (optional) complex vector of length N, containing the
%           eigenvalues of A, possibly in the required order.
%           When task = 1, ev is not available.
%
%   Comments
%   1. For task = 1, when JOB = 1, 2, or 3, the scaling factors are
%      computed using only A, A and B, or A and C, respectively, but
%      the transformation is applied to the entire system matrix S.
%   2. For task = 1, U is N-by-1 storing the diagonal elements of D.
%      For task = 2, 3, U is orthogonal.
%      For task = 4, U is nonsingular.
%      Matrix U is returned optionally.
%
% See also SLSBAL, SLSDEC, SLSORSF, SLSRSF
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   H. Xu, TU Chemnitz, FR Germany, Dec. 1998
%
% Revisions:
%   V. Sima, Katholieke Univ. Leuven, Belgium, May 1999.
