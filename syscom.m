%SYSCOM  MEX-function for computing controllability, observability,
%        and minimality forms of a given system using SLICOT routines
%        AB07MD, TB01PD, TB01UD, and TB01ZD.
%
%   [A,B,C(,N,U,sizes)] = syscom(task,A,B,C,tol(,bal))
%
%   [A,B,C,Nc,U,sizes] = syscom(1,A,B,C,tol)
%   [A,B,C,No,U,sizes] = syscom(2,A,B,C,tol)
%           [Ar,Br,Cr] = syscom(3,A,B,C,tol,bal)
%
%   SYSCOM transforms the matrix triple (A, B, C) to various staircase
%   forms showing controllability, observability or minimality:
%
%   a) Controllability (block) form
%
%                                |-  Nc -|
%                           ---[ * . . . * * . . * | * ]
%                           |  [ *       . .     . | 0 ]
%                           Nc [   .     . .     . | . ]
%                           |  [     .   . .     . | . ]
%       [ U'AU | U'B ]      ---[       * * * . . * | . ]
%       [------|-----] =       [           * . . * | . ];          (1)
%       [  CU  |  0  ]         [           .     . | . ]
%                              [           * . . * | 0 ]
%                              [-------------------|---]
%                              [ * . . . * * . . * | 0 ]
%
%   b) Observability (block) form
%
%                                |-  No -|
%                           ---[ * *               | * ]
%                           |  [ .   .             | * ]
%                           No [ .     .           | . ]
%                           |  [ .       *         | . ]
%       [ U'AU | U'B ]      ---[ * . . . *         | . ]
%       [------|-----] =       [ * . . . * * . . * | . ];          (2)
%       [  CU  |  0  ]         [ .       . .     . | . ]
%                              [ * . . . * * . . * | * ]
%                              [-------------------|---]
%                              [ * 0 . . . . . . 0 | 0 ]
%
%   c) Minimality (block) form
%
%                            |--Nm --|
%                        --[ * . . . *         * . . . * | * ]
%                        | [ *       .         .       . | . ]
%                        Nm[   .  Ar .         .       . | Br]
%                        | [     .   .         .       . | . ]
%                        --[       * *         * . . . * | * ]
%                          [ * . . . * * . . * * . . . * | * ]
%          [ Ar | Br ]     [ .       . .     . .       . | . ]
%          [----|----] =   [ .       . .     . .       . | . ].    (3)
%          [ Cr | 0  ]     [ * . . . * * . . * * . . . * | * ]
%                          [                   * . . . * | 0 ]
%                          [                   .       . | . ]
%                          [                   * . . . * | 0 ]
%                          [-----------------------------|---]
%                          [* . Cr . * 0 . . 0 * . . . * | 0 ]
%
%   Description of other input parameters:
%   task  - integer option to determine which form is computed:
%           = 1 : compute the controllability form (1);
%           = 2 : compute the observabilty form (2);
%           = 3 : compute the minimal realization subsystem (Ar,Br,Cr)
%                 from (3).
%   tol   - (optional) real tolerance value used for rank decisions,
%           as a lower bound on the reciprocal condition numbers, or
%           for checking controllability or observability for
%           single-input or single-output systems, respectively.
%           Default:  tol = N*N*machine_epsilon, for rank decisions;
%                     tol = N*machine_epsilon*max(norm(A),norm(BC)), for
%                           checking controllability or observability,
%                           with BC = B or BC = C, respectively.
%   bal   - (for task = 3 only: optional) integer indicating whether
%           or not to balance the system triple before computing the
%           minimal realization subsystem.
%           bal = 0 : use balancing;
%                     otherwise, do not use balancing.
%           Default:  bal = 0.
%
%   Description of other output parameters:
%   Nc    - (optional) order of the controllable subsystem when
%           task = 1.
%   No    - (optional) order of the observable subsystem when
%           task = 2.
%   sizes - (optional) integer vector indicating the sizes of the
%           diagonal blocks of A corresponding to controllable or
%           observable parts, respectively. The length of sizes
%           indicates the associated controllable or observable index.
%
% See also SLCONF, SLMINR, SLOBSF, SLSBAL, SYSTRA
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   H. Xu, TU Chemnitz, FR Germany, Dec. 1998.
%
% Revisions:
%   V. Sima, Katholieke Univ. Leuven, Belgium, May 1999.
