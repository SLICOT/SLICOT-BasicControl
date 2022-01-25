% DEADBEAT.F - MEX-function for constructing the minimum norm
%              feedback matrix F to perform "deadbeat control" on
%              a (A,B)-pair, using SLICOT routines AB01OD and SB06ND.
%
%   [(F)(,Ao,Bo(,kstair),U(,V),scale)] = 
%                        deadbeat(A,B(,IStair(,kstair)(,tol),
%                                     WithU(,U1),WithV,bal,scale))
%
%   DEADBEAT constructs the minimum norm feedback matrix F performing
%   "deadbeat control" on a (A,B)-pair of a state-space model.
%   The (A,B)-pair must be controllable, unless the uncontrollable
%   part has all poles in the origin. Optionally, the orthogonal
%   canonical form (also called "staircase" form) of the (A,B)-pair,
%   or the staircase form with upper triangular stairs (also called
%   upper staircase form), can be returned.
%
%   Description of input parameters:
%   A      - the n-by-n state matrix A.
%            If IStair = -1 or IStair = 0, A contains the initial
%            matrix Ai of the original system.
%            If IStair = -2 or IStair = 1, A must be in the staircase
%            form, U1'*Ai*U1, as produced by SLICOT Library routine
%            AB01ND, or by this MEX-file with option IStair = -1.
%            If IStair = 2, A must be the transformed state-space
%            matrix U1'*Ai*U1 of the (A,B)-pair with triangular stairs,
%            as produced by SLICOT Library routine AB01OD (with option
%            STAGES = 'A'), or by this MEX-file with option IStair = -2.
%   B      - the n-by-m input matrix B.
%            If IStair = -1 or IStair = 0, B contains the initial
%            matrix Bi of the original system.
%            If IStair = -2 or IStair = 1, B must be in the staircase
%            form, U1'*Bi, as produced by SLICOT Library routine AB01ND,
%            or by this MEX-file with option IStair = -1.
%            If IStair = 2, B must be the transformed triangular input
%            matrix U1'*Bi*V of the (A,B)-pair as produced by SLICOT
%            Library routine AB01OD (with option STAGES = 'A'), or by
%            this MEX-file with option IStair = -2.
%   IStair - (optional) scalar indicating whether the (A,B)-pair is
%            already in the staircase form with or without triangular
%            stairs, or if such a form should be returned:
%            =-2 :  (A,B)-pair is in the staircase form and it should
%                   be returned in the upper staircase form;
%            =-1 :  (A,B)-pair is general and it should be returned in
%                   the staircase form;
%            = 0 :  (A,B)-pair is general;
%            = 1 :  (A,B)-pair is in the staircase form;
%            = 2 :  (A,B)-pair is in the upper staircase form.
%            Default: IStair = 0.
%            If IStair < 0, the deadbeat feedback matrix is not found.
%   kstair - (optional) if IStair = -2 or IStair >= 1, integer vector
%            containing the dimensions of each "stair", or, also, the
%            orders of the diagonal blocks of the controllable part
%            of A, Acont (see Method). If IStair = -1 or IStair = 0,
%            kstair must not be specified as an input parameter.
%   tol    - (optional) if IStair = -1 or IStair = 0, real scalar
%            containing the tolerance to be used in rank determination
%            when transforming (A, B). If tol > 0, then this value is
%            used as a lower bound for the reciprocal condition number;
%            a (sub)matrix whose estimated condition number is less than
%            1/tol is considered to be of full rank.  If tol <= 0, then
%            n*n*eps is used instead, where eps is the machine
%            precision. If IStair = -2 or IStair >= 1, tol must not be
%            specified as an input parameter.
%            Default: tol = 0.
%   WithU  - (optional) scalar indicating whether the matrix U should be
%            computed and returned:
%            = 0 :  do not form and return U;
%            = 1 :  form and return U.
%            Default: WithU = 1.
%            Note: The matrix U must be computed in order to find the
%            deadbeat feedback matrix in terms of the original state-
%            space coordinates.
%   U1     - (optional) if IStair = -2 or IStair >= 1 and WithU = 1, the
%            given n-by-n matrix U1.
%   WithV  - (optional) if IStair = -2 or IStair = 0 or IStair = 1,
%            scalar indicating whether the matrix V should be computed
%            and returned:
%            = 0 :  do not form and return V;
%            = 1 :  form and return V.
%            Default: WithV = 1.
%            Note: If IStair = -2 or IStair = 0 or IStair = 1, the
%            matrix V must be computed in order to find the deadbeat
%            feedback matrix in terms of the original state-space
%            coordinates.
%   bal    - (optional) if IStair >= -1, integer indicating whether
%            the (A,B)-pair should be balanced (if IStair = -1 or
%            IStair = 0), or the previously computed scaling factors
%            should be used (if IStair = 1 or IStair = 2).
%            = 0 :  use balancing;
%            = 1 :  do not use balancing.
%            Default: bal = 0.
%            If IStair = -2, bal must not be specified as an input
%            parameter.
%   scale  - (optional) if IStair > 0 and bal = 0, real n-vector
%            containing the previously computed scaling factors.
%            If IStair <= 0, or bal = 1, scale must not be specified as
%            an input parameter.
%            Note: If IStair > 0 and bal = 0, the vector scale must be
%            given in order to find the deadbeat feedback matrix in
%            terms of the original state-space coordinates, if the
%            original system was scaled.
%
%   Description of output parameters:
%   F      - (optional) if IStair >= 0, the m-by-n deadbeat feedback
%            matrix F.
%            If WithU = 0 (and WithV = 0 for IStair = 0 or IStair = 1),
%            F contains the deadbeat feedback matrix Fo, in terms of the
%            reduced system coordinates (see Method).
%            If WithU = 1 (and WithV = 1 for IStair = 0 or IStair = 1),
%            F contains the deadbeat feedback matrix in terms of the
%            original state-space coordinates.
%            The deadbeat feedback matrix is scaled (F <- F*inv(D),
%            where D = diag(scale)), if balancing was performed and
%            the vector scale is available.
%            If IStair < 0, F must not be specified as an output
%            parameter.
%   Ao     - (optional) if IStair < 0, the n-by-n matrix U'*A*U, in
%            staircase form, if IStair = -1, or in upper staircase form,
%            if IStair = -2.
%            If IStair >= 0, the n-by-n matrix U'*A*U + U'*B*V*Fo.
%            A and B are the scaled matrices (A <- inv(D)*A*D, 
%            B <- inv(D)*B), if balancing was performed.
%   Bo     - (optional) if IStair = -1, the n-by-m matrix U'*B in the
%            staircase form.
%            If IStair = -2 or IStair >= 0, the n-by-m matrix U'*B*V, in
%            upper staircase form.
%            B is the scaled matrix if balancing was performed.
%   kstair - (optional) if IStair = -1 or IStair = 0, integer vector
%            containing the dimensions of each "stair". If IStair = -2
%            or IStair >= 1, kstair must not be specified as an output
%            parameter.
%   U      - (optional) if WithU = 1, the n-by-n matrix U.
%            If IStair = -2 or IStair >= 1, U is the product of the
%            input matrix U1 and the computed state-space transformation
%            matrix.
%   V      - (optional) if IStair = -2 or IStair = 0 or IStair = 1 and
%            WithV = 1, the m-by-m matrix V. If IStair = -1 or
%            IStair = 2, V must not be specified as an output parameter.
%   scale  - (optional) if IStair = -1 or IStair = 0 and bal = 0, the
%            n-vector of scaling factors.
%            If IStair = -2 or IStair >= 1, or bal = 1, scale must not
%            be specified as an output parameter.
%
%   Method:
%   If IStair = 0, the matrices A and B are reduced using (and
%   optionally accumulating) state-space and input-space transformations
%   U1 and V respectively, such that the pair of matrices
% 
%      Ac = U1'*A*U1,    Bc = U1'*B*V
%
%   are in upper "staircase" form. Specifically,
%
%           [ Acont     *    ]         [ Bcont ]
%      Ac = [                ],   Bc = [       ],                    (1)
%           [   0    Auncont ]         [   0   ]
%
%      and
%
%              [ A11 A12  . . .  A1,p-1 A1p ]         [ B1 ]
%              [ A21 A22  . . .  A2,p-1 A2p ]         [ 0  ]
%              [  0  A32  . . .  A3,p-1 A3p ]         [ 0  ]
%      Acont = [  .   .   . . .    .     .  ],   Bc = [ .  ],        (2)
%              [  .   .     . .    .     .  ]         [ .  ]
%              [  .   .       .    .     .  ]         [ .  ]
%              [  0   0   . . .  Ap,p-1 App ]         [ 0  ]
%
%   where the blocks  B1, A21, ..., Ap,p-1  have full row ranks and
%   p is the controllability index of the pair.  The size of the
%   block Auncont is equal to the dimension of the uncontrollable
%   subspace of the pair (A, B).  The first stage of the reduction,
%   the "forward" stage, accomplishes the reduction to the orthogonal
%   canonical form. The matrix V is an identity matrix for this stage.
%   The blocks B1, A21, ..., Ap,p-1 are further reduced in a second,
%   "backward" stage to upper triangular form using RQ factorization.
%
%   If IStair = -1, the forward step only is performed, and the matrices
%   Ac and Bc in the staircase form (1) and (2), with identity V, are
%   returned in Ao and Bo, respectively.
%
%   If IStair = -2 or IStair = 1, the matrices A and B are assumed to be
%   given in the staircase form (1) and (2), with identity V, and the
%   backward step only is performed.
%
%   If IStair < 0, the deadbeat feedback matrix is not computed.
%
%   If IStair = 2, the matrices A and B are assumed to be given in the
%   staircase form (1) and (2) with upper triangular stairs, and the 
%   deadbeat feedback matrix is directly computed.
%
%   The deadbeat feedback matrix F can only be computed if the submatrix
%   Auncont in (1) is either empty or zero. The controllable subsystem
%   (Acont,Bcont) is procesed as described below, where (Acont,Bcont) is
%   redefined as (A,B).
% 
%   Starting from the (A,B)-pair in "staircase form" with "triangular"
%   stairs, dimensions kstair(i+1) x kstair(i):
%
%                  [ B1 | A11   A12 * . .  *  ]
%                  [    | A21   A22   .    .  ]
%                  [    |    .      .   .  .  ]
%    [ B | A ]  =  [    |      .      .       ]
%                  [    |        .      .  *  ]
%                  [ 0  |   0                 ]
%                  [    |          Ap,p-1 App ]
%
%   where the i-th diagonal block of A has dimension kstair(i), for
%   i = 1,2,...,p, the feedback matrix F is constructed recursively in
%   p steps. In each step a unitary state-space transformation U2 and a
%   part of Fo (in the reduced system coordinates) are updated in order
%   to achieve the final form:
%
%                             [ 0   A12    *   . . .  *    ]
%                             [                .      .    ]
%                             [     0    A23     .    .    ]
%                             [         .     .            ]
%   [ U2'*A*U2 + U2'*B*Fo ] = [           .     .     *    ] .
%                             [             .              ]
%                             [               .     Ap-1,p ]
%                             [                            ]
%                             [                       0    ]
%
%   If IStair = -1 or IStair = 0, and bal = 0, the matrices A and B are
%   preliminarily scaled, A <- inv(D)*A*D, B <- inv(D)*B. For the whole
%   procedure (IStair = 0), U = U1*U2, and F = V*Fo*U'*inv(D).
%
%   Comments
%   The (A,B)-pair must be controllable, unless the uncontrollable
%   part submatrix Auncont is zero.
%   The eigenvalues of the closed-loop matrix, A+B*F, could be far from
%   zero (especially for large values of n - m), since the deadbeat
%   problem could be ill-conditioned.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Sep. 2003.
%
% Revisions:
%   -
%
