#include "fintrf.h"
C DEADBEAT.F - Gateway function for constructing the minimum norm
C              feedback matrix F to perform "deadbeat control" on
C              a (A,B)-pair, using SLICOT routines AB01OD and SB06ND.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [(F)(,Ao,Bo(,kstair),U(,V),scale)] =
C                        deadbeat(A,B(,IStair(,kstair)(,tol),
C                                     WithU(,U1),WithV,bal,scale))
C
C Purpose:
C   To construct the minimum norm feedback matrix F performing
C   "deadbeat control" on a (A,B)-pair of a state-space model.
C   The (A,B)-pair must be controllable, unless the uncontrollable
C   part has all poles in the origin. Optionally, the orthogonal
C   canonical form (also called "staircase" form) of the (A,B)-pair,
C   or the staircase form with upper triangular stairs (also called
C   upper staircase form), can be returned.
C
C Input parameters:
C   A      - the n-by-n state matrix A.
C            If IStair = -1 or IStair = 0, A contains the initial
C            matrix Ai of the original system.
C            If IStair = -2 or IStair = 1, A must be in the staircase
C            form, U1'*Ai*U1, as produced by SLICOT Library routine
C            AB01ND, or by this MEX-file with option IStair = -1.
C            If IStair = 2, A must be the transformed state-space
C            matrix U1'*Ai*U1 of the (A,B)-pair with triangular stairs,
C            as produced by SLICOT Library routine AB01OD (with option
C            STAGES = 'A'), or by this MEX-file with option IStair = -2.
C   B      - the n-by-m input matrix B.
C            If IStair = -1 or IStair = 0, B contains the initial
C            matrix Bi of the original system.
C            If IStair = -2 or IStair = 1, B must be in the staircase
C            form, U1'*Bi, as produced by SLICOT Library routine AB01ND,
C            or by this MEX-file with option IStair = -1.
C            If IStair = 2, B must be the transformed triangular input
C            matrix U1'*Bi*V of the (A,B)-pair as produced by SLICOT
C            Library routine AB01OD (with option STAGES = 'A'), or by
C            this MEX-file with option IStair = -2.
C   IStair - (optional) scalar indicating whether the (A,B)-pair is
C            already in the staircase form with or without triangular
C            stairs, or if such a form should be returned:
C            =-2 :  (A,B)-pair is in the staircase form and it should
C                   be returned in the upper staircase form;
C            =-1 :  (A,B)-pair is general and it should be returned in
C                   the staircase form;
C            = 0 :  (A,B)-pair is general;
C            = 1 :  (A,B)-pair is in the staircase form;
C            = 2 :  (A,B)-pair is in the upper staircase form.
C            Default: IStair = 0.
C            If IStair < 0, the deadbeat feedback matrix is not found.
C   kstair - (optional) if IStair = -2 or IStair >= 1, integer vector
C            containing the dimensions of each "stair", or, also, the
C            orders of the diagonal blocks of the controllable part
C            of A, Acont (see Method). If IStair = -1 or IStair = 0,
C            kstair must not be specified as an input parameter.
C   tol    - (optional) if IStair = -1 or IStair = 0, real scalar
C            containing the tolerance to be used in rank determination
C            when transforming (A, B). If tol > 0, then this value is
C            used as a lower bound for the reciprocal condition number;
C            a (sub)matrix whose estimated condition number is less than
C            1/tol is considered to be of full rank.  If tol <= 0, then
C            n*n*eps is used instead, where eps is the machine
C            precision. If IStair = -2 or IStair >= 1, tol must not be
C            specified as an input parameter.
C            Default: tol = 0.
C   WithU  - (optional) scalar indicating whether the matrix U should be
C            computed and returned:
C            = 0 :  do not form and return U;
C            = 1 :  form and return U.
C            Default: WithU = 1.
C            Note: The matrix U must be computed in order to find the
C            deadbeat feedback matrix in terms of the original state-
C            space coordinates.
C   U1     - (optional) if IStair = -2 or IStair >= 1 and WithU = 1, the
C            given n-by-n matrix U1.
C   WithV  - (optional) if IStair = -2 or IStair = 0 or IStair = 1,
C            scalar indicating whether the matrix V should be computed
C            and returned:
C            = 0 :  do not form and return V;
C            = 1 :  form and return V.
C            Default: WithV = 1.
C            Note: If IStair = -2 or IStair = 0 or IStair = 1, the
C            matrix V must be computed in order to find the deadbeat
C            feedback matrix in terms of the original state-space
C            coordinates.
C   bal    - (optional) if IStair >= -1, integer indicating whether
C            the (A,B)-pair should be balanced (if IStair = -1 or
C            IStair = 0), or the previously computed scaling factors
C            should be used (if IStair = 1 or IStair = 2).
C            = 0 :  use balancing;
C            = 1 :  do not use balancing.
C            Default: bal = 0.
C            If IStair = -2, bal must not be specified as an input
C            parameter.
C   scale  - (optional) if IStair > 0 and bal = 0, real n-vector
C            containing the previously computed scaling factors.
C            If IStair <= 0, or bal = 1, scale must not be specified as
C            an input parameter.
C            Note: If IStair > 0 and bal = 0, the vector scale must be
C            given in order to find the deadbeat feedback matrix in
C            terms of the original state-space coordinates, if the
C            original system was scaled.
C
C Output parameters:
C   F      - (optional) if IStair >= 0, the m-by-n deadbeat feedback
C            matrix F.
C            If WithU = 0 (and WithV = 0 for IStair = 0 or IStair = 1),
C            F contains the deadbeat feedback matrix Fo, in terms of the
C            reduced system coordinates (see Method).
C            If WithU = 1 (and WithV = 1 for IStair = 0 or IStair = 1),
C            F contains the deadbeat feedback matrix in terms of the
C            original state-space coordinates.
C            The deadbeat feedback matrix is scaled (F <- F*inv(D),
C            where D = diag(scale)), if balancing was performed and
C            the vector scale is available.
C            If IStair < 0, F must not be specified as an output
C            parameter.
C   Ao     - (optional) if IStair < 0, the n-by-n matrix U'*A*U, in
C            staircase form, if IStair = -1, or in upper staircase form,
C            if IStair = -2.
C            If IStair >= 0, the n-by-n matrix U'*A*U + U'*B*V*Fo.
C            A and B are the scaled matrices (A <- inv(D)*A*D,
C            B <- inv(D)*B), if balancing was performed.
C   Bo     - (optional) if IStair = -1, the n-by-m matrix U'*B in the
C            staircase form.
C            If IStair = -2 or IStair >= 0, the n-by-m matrix U'*B*V, in
C            upper staircase form.
C            B is the scaled matrix if balancing was performed.
C   kstair - (optional) if IStair = -1 or IStair = 0, integer vector
C            containing the dimensions of each "stair". If IStair = -2
C            or IStair >= 1, kstair must not be specified as an output
C            parameter.
C   U      - (optional) if WithU = 1, the n-by-n matrix U.
C            If IStair = -2 or IStair >= 1, U is the product of the
C            input matrix U1 and the computed state-space transformation
C            matrix.
C   V      - (optional) if IStair = -2 or IStair = 0 or IStair = 1 and
C            WithV = 1, the m-by-m matrix V. If IStair = -1 or
C            IStair = 2, V must not be specified as an output parameter.
C   scale  - (optional) if IStair = -1 or IStair = 0 and bal = 0, the
C            n-vector of scaling factors.
C            If IStair = -2 or IStair >= 1, or bal = 1, scale must not
C            be specified as an output parameter.
C
C Method:
C   If IStair = 0, the matrices A and B are reduced using (and
C   optionally accumulating) state-space and input-space transformations
C   U1 and V respectively, such that the pair of matrices
C
C      Ac = U1'*A*U1,    Bc = U1'*B*V
C
C   are in upper "staircase" form. Specifically,
C
C           [ Acont     *    ]         [ Bcont ]
C      Ac = [                ],   Bc = [       ],                    (1)
C           [   0    Auncont ]         [   0   ]
C
C      and
C
C              [ A11 A12  . . .  A1,p-1 A1p ]         [ B1 ]
C              [ A21 A22  . . .  A2,p-1 A2p ]         [ 0  ]
C              [  0  A32  . . .  A3,p-1 A3p ]         [ 0  ]
C      Acont = [  .   .   . . .    .     .  ],   Bc = [ .  ],        (2)
C              [  .   .     . .    .     .  ]         [ .  ]
C              [  .   .       .    .     .  ]         [ .  ]
C              [  0   0   . . .  Ap,p-1 App ]         [ 0  ]
C
C   where the blocks  B1, A21, ..., Ap,p-1  have full row ranks and
C   p is the controllability index of the pair.  The size of the
C   block Auncont is equal to the dimension of the uncontrollable
C   subspace of the pair (A, B).  The first stage of the reduction,
C   the "forward" stage, accomplishes the reduction to the orthogonal
C   canonical form. The matrix V is an identity matrix for this stage.
C   The blocks B1, A21, ..., Ap,p-1 are further reduced in a second,
C   "backward" stage to upper triangular form using RQ factorization.
C
C   If IStair = -1, the forward step only is performed, and the matrices
C   Ac and Bc in the staircase form (1) and (2), with identity V, are
C   returned in Ao and Bo, respectively.
C
C   If IStair = -2 or IStair = 1, the matrices A and B are assumed to be
C   given in the staircase form (1) and (2), with identity V, and the
C   backward step only is performed.
C
C   If IStair < 0, the deadbeat feedback matrix is not computed.
C
C   If IStair = 2, the matrices A and B are assumed to be given in the
C   staircase form (1) and (2) with upper triangular stairs, and the
C   deadbeat feedback matrix is directly computed.
C
C   The deadbeat feedback matrix F can only be computed if the submatrix
C   Auncont in (1) is either empty or zero. The controllable subsystem
C   (Acont,Bcont) is procesed as described below, where (Acont,Bcont) is
C   redefined as (A,B).
C
C   Starting from the (A,B)-pair in "staircase form" with "triangular"
C   stairs, dimensions kstair(i+1) x kstair(i):
C
C                  [ B1 | A11   A12 * . .  *  ]
C                  [    | A21   A22   .    .  ]
C                  [    |    .      .   .  .  ]
C    [ B | A ]  =  [    |      .      .       ]
C                  [    |        .      .  *  ]
C                  [ 0  |   0                 ]
C                  [    |          Ap,p-1 App ]
C
C   where the i-th diagonal block of A has dimension kstair(i), for
C   i = 1,2,...,p, the feedback matrix F is constructed recursively in
C   p steps. In each step a unitary state-space transformation U2 and a
C   part of Fo (in the reduced system coordinates) are updated in order
C   to achieve the final form:
C
C                             [ 0   A12    *   . . .  *    ]
C                             [                .      .    ]
C                             [     0    A23     .    .    ]
C                             [         .     .            ]
C   [ U2'*A*U2 + U2'*B*Fo ] = [           .     .     *    ] .
C                             [             .              ]
C                             [               .     Ap-1,p ]
C                             [                            ]
C                             [                       0    ]
C
C   If IStair = -1 or IStair = 0, and bal = 0, the matrices A and B are
C   preliminarily scaled, A <- inv(D)*A*D, B <- inv(D)*B. For the whole
C   procedure (IStair = 0), U = U1*U2, and F = V*Fo*U'*inv(D).
C
C Comments
C   The (A,B)-pair must be controllable, unless the uncontrollable
C   part submatrix Auncont is zero.
C   The eigenvalues of the closed-loop matrix, A+B*F, could be far from
C   zero (especially for large values of n - m), since the deadbeat
C   problem could be ill-conditioned.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Sep. 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2005,
C   Apr. 2009, Dec. 2010, July 2012, Aug. 2012, Oct. 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HUND
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, HUND = 1.0D2 )
C
C .. Mex-file interface parameters ..
      mwPointer         PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      mwPointer         mxCalloc, mxCreateDoubleMatrix,
     $                  mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         JOBU, JOBV, STAGES
      INTEGER           INDCON, INFO, INFO1, LDA, LDB, LDF, LDU, LDV,
     $                  LDWORK, M, N, NCONT
      DOUBLE PRECISION  MAXRED, TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:), KSTAIR(:)
      DOUBLE PRECISION, ALLOCATABLE :: DK(:)
      mwPointer         A, B, DW, DWORK, F, FS, SCAL, U, V
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           BAL, I, IP, IPS, ISTAIR, NV, WITHU, WITHV
      DOUBLE PRECISION  TEMP
C
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE
C
C .. External subroutines ..
      EXTERNAL          AB01OD, DCOPY, DGEMM, DLACPY, DLASET, MA02KD,
     $                  MB01SD, SB06ND, TB01ID
C
C .. Intrinsic functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'DEADBEAT requires at least 2 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'DEADBEAT requires at least 1 output argument' )
      END IF
C
C   A(nxn), B(nxm), IStair, kstair(indcon), tol, WithU, U1(nxn), WithV,
C   bal, scale(n).
C
      N = mxGetM( PRHS(1) )
      M = mxGetN( PRHS(2) )
C
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF ( mxGetN( PRHS(1) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF ( mxGetM( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same number of rows as A' )
      END IF
C
      IP = 3
      IF ( NRHS.GE.IP ) THEN
C
C   IStair
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ISTAIR must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ISTAIR must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ISTAIR = TEMP
         IF ( ABS( ISTAIR ).GT.2 ) THEN
            CALL mexErrMsgTxt
     $      ( 'ISTAIR has -2, -1, 0, 1 or 2 the only admissible values')
         END IF
         IP = IP + 1
      ELSE
         ISTAIR = 0
      END IF
C
      IF ( ISTAIR.LT.0 .AND. NLHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'DEADBEAT requires at least 2 output arguments' )
      END IF
C
      INDCON = N
      IF ( NRHS.GE.IP ) THEN
         IF ( ISTAIR.EQ.-1 .OR. ISTAIR.EQ.0 ) THEN
C
C   tol
C
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'TOL must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'TOL must be a real scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TOL, 1 )
         ELSE
C
C   kstair
C
            INDCON = mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) )
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'KSTAIR must be an integer vector' )
            END IF
            IF ( INDCON.GT.N ) THEN
               WRITE( TEXT, '('' KSTAIR must have at most '',I7,
     $                        '' entries'')' ) N
               CALL mexErrMsgTxt( TEXT )
            END IF
         END IF
         IP = IP + 1
      ELSE
         TOL = ZERO
      END IF
C
      IF ( NRHS.GE.IP ) THEN
C
C   WithU
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'WITHU must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'WITHU must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         WITHU = TEMP
         IF ( WITHU.LT.0 .OR. WITHU.GT.1 ) THEN
            CALL mexErrMsgTxt
     $        ( 'WITHU has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
      ELSE
         WITHU = 1
      END IF
C
      IF ( WITHU.EQ.0 ) THEN
         JOBU = 'N'
      ELSE
         JOBU = 'I'
      END IF
C
      IF ( NRHS.GE.IP .AND. ( ISTAIR.EQ.-2 .OR. ISTAIR.GE.1 ) .AND.
     $     WITHU.EQ.1 ) THEN
C
C   U1
C
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'U1 must be a real matrix' )
         END IF
         IF ( mxGetM( PRHS(IP) ).NE.N .OR.
     $        mxGetN( PRHS(IP) ).NE.N ) THEN
            CALL mexErrMsgTxt( 'U1 must have the same order as A' )
         END IF
         IP = IP + 1
      END IF
C
      IF ( ISTAIR.NE.-1 .AND. ISTAIR.NE.2 ) THEN
         IF ( NRHS.GE.IP ) THEN
C
C   WithV
C
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'WITHV must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'WITHV must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            WITHV = TEMP
            IF ( WITHV.LT.0 .OR. WITHV.GT.1 ) THEN
               CALL mexErrMsgTxt
     $            ( 'WITHV has 0 or 1 the only admissible values' )
            END IF
            IP = IP + 1
         ELSE
            WITHV = 1
         END IF
      ELSE
         WITHV = 0
      END IF
C
      IF ( WITHV.EQ.0 ) THEN
         JOBV = 'N'
      ELSE
         JOBV = 'I'
      END IF
C
      IF ( NRHS.GE.IP .AND. ISTAIR.GE.-1 ) THEN
C
C   bal
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'BAL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'BAL must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         BAL = TEMP
         IF ( BAL.LT.0 .OR. BAL.GT.1 ) THEN
            CALL mexErrMsgTxt
     $         ( 'BAL has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
      ELSE
         BAL = 0
      END IF
C
      IF ( NRHS.GE.IP .AND. ISTAIR.GT.0 .AND. BAL.EQ.0 ) THEN
C
C   scale
C
         IF ( mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) ).NE.N ) THEN
            WRITE( TEXT, '('' SCALE must be a vector with '',I7,
     $                     '' elements'')' ) N
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'SCALE must be a real vector' )
         END IF
         IPS = IP
      ELSE
         IPS = 10
      END IF
C
C Determine the lenghts of working arrays.
C For ISTAIR <= 1, use a larger value for LDWORK for enabling calls of
C block algorithms in AB01OD.
C
      LDA = MAX( 1, N )
      LDB = LDA
      LDU = LDA
      LDF = MAX( 1, M )
      IF ( WITHV.EQ.1 ) THEN
         LDV = LDF
         NV  = M
      ELSE
         LDV = 1
         NV  = 0
      END IF
      IF ( ISTAIR.EQ.-2 ) THEN
         LDWORK = MAX( 1, M + MAX( N, M ) )
      ELSE IF ( ISTAIR.EQ.1 ) THEN
         LDWORK = MAX( 1, M + MAX( N, M ), 2*N )
      ELSE IF ( ISTAIR.LE.0 ) THEN
         LDWORK = MAX( 1, N + MAX( N, 3*M ) )
      ELSE
         LDWORK = 2*N
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( KSTAIR( INDCON ) )
      A     = mxCalloc(  LDA*N, 8 )
      B     = mxCalloc(  LDB*M, 8 )
      DWORK = mxCalloc( LDWORK, 8 )
      U     = mxCalloc(  LDU*N, 8 )
      IF ( ISTAIR.GE.-1 .AND. BAL.EQ.0 )
     $   SCAL = mxCalloc( N, 8 )
      IF ( ISTAIR.GE.0 )
     $   F = mxCalloc( LDF*N, 8 )
      IF ( ISTAIR.NE.2 )
     $   V = mxCalloc( LDV*NV, 8 )
      IF ( ISTAIR.EQ.-1 .OR. ISTAIR.EQ.0 ) THEN
         ALLOCATE ( IWORK( M ) )
      ELSE IF ( ISTAIR.LE.1 ) THEN
         ALLOCATE ( IWORK( 0 ) )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), %VAL( A ), N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( B ), N*M )
      IF ( ISTAIR.EQ.-2 .OR. ISTAIR.GE.1 ) THEN
         ALLOCATE ( DK( INDCON ) )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), DK, INDCON )
         NCONT = 0
C
         DO 10 I = 1, INDCON
            KSTAIR(I) = DK(I)
            NCONT = NCONT + KSTAIR(I)
   10    CONTINUE
C
         DEALLOCATE ( DK )
         IF ( NCONT.GT.N ) THEN
            WRITE( TEXT, '('' SUM( KSTAIR ) must not exceed '', I7)' ) N
            CALL mexErrMsgTxt( TEXT )
         END IF
C
         IF ( NRHS.GE.6 .AND. WITHU.EQ.1 ) THEN
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), %VAL( U ),
     $                             N*N )
         ELSE IF ( WITHU.EQ.1 .OR. ISTAIR.EQ.2 ) THEN
            CALL DLASET( 'Full', N, N, ZERO, ONE, %VAL( U ), LDU )
         END IF
      END IF
      IF ( WITHU.EQ.0 .AND. ( ISTAIR.EQ.0 .OR. ISTAIR.EQ.1 ) )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, %VAL( U ), LDU )
C
      IF ( NRHS.GE.IPS .AND. ISTAIR.GT.0 .AND. BAL.EQ.0 )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(IPS) ),
     $                          %VAL( SCAL ), N )
C
C Do the actual computations.
C
      IF ( ( ISTAIR.EQ.-1 .OR. ISTAIR.EQ.0 ) .AND. BAL.EQ.0 ) THEN
         MAXRED = ZERO
         CALL TB01ID( 'B and A', N, M, 0, MAXRED, %VAL( A ), LDA,
     $                %VAL( B ), LDB, %VAL( DWORK ), 1, %VAL( SCAL ),
     $                INFO )
      END IF
C
      IF ( ISTAIR.LE.1 ) THEN
         IF ( ISTAIR.EQ.0 ) THEN
            STAGES = 'A'
         ELSE IF ( ISTAIR.EQ.-1 ) THEN
            STAGES = 'F'
         ELSE
            STAGES = 'B'
         END IF
         CALL AB01OD( STAGES, JOBU, JOBV, N, M, %VAL( A ), LDA,
     $                %VAL( B ), LDB, %VAL( U ), LDU, %VAL( V ), LDV,
     $                NCONT, INDCON, KSTAIR, TOL, IWORK, %VAL( DWORK ),
     $                LDWORK, INFO1 )
         IF ( ISTAIR.GE.0 .AND. NCONT.NE.N ) THEN
            DW = mxCalloc( (N-NCONT)*(N-NCONT), 8 )
            CALL MA02KD( 'Full', N-NCONT, N-NCONT, NCONT+1, NCONT+1, 1,
     $                   1, %VAL( A ), LDA, %VAL( DW ), N-NCONT )
            IF( DLANGE( '1-norm', N-NCONT, N-NCONT, %VAL( DW ), N-NCONT,
     $                  %VAL( DW ) ).GT.HUND*DBLE( N )*DLAMCH( 'Epsilon'
     $                  ) ) THEN
               CALL mexErrMsgTxt
     $            ( 'The deadbeat feedback matrix cannot be computed (e.
     $g., the system is not controllable).' )
            END IF
            CALL mxFree( DW  )
         END IF
      ELSE
         INFO1 = 0
      END IF
C
      IF ( ISTAIR.GE.0 .AND. INFO1.EQ.0 ) THEN
         CALL SB06ND( N, M, INDCON, %VAL( A ), LDA, %VAL( B ), LDB,
     $                KSTAIR, %VAL( U ), LDU, %VAL( F ), LDF,
     $                %VAL( DWORK ), INFO )
         IF ( MIN( M, NCONT ).GT.0 .AND. MAX( WITHU, WITHV ).EQ.1 ) THEN
            FS = mxCalloc( LDF*NCONT, 8 )
            IF ( WITHV.EQ.1 ) THEN
               CALL DLACPY( 'Full', M, NCONT, %VAL( F ), LDF,
     $                      %VAL( FS ), LDF )
               CALL DGEMM(  'NoTranspose', 'NoTranspose', M, NCONT, M,
     $                      ONE, %VAL( V ), LDV, %VAL( FS ), LDF, ZERO,
     $                      %VAL( F ), LDF )
            END IF
            IF ( WITHU.EQ.1 ) THEN
               CALL DLACPY( 'Full', M, NCONT, %VAL( F ), LDF,
     $                      %VAL( FS ), LDF )
               CALL DGEMM(  'NoTranspose', 'Transpose', M, N, NCONT,
     $                      ONE, %VAL( FS ), LDF, %VAL( U ), LDU, ZERO,
     $                      %VAL( F ), LDF )
            END IF
            CALL mxFree( FS )
         END IF
C
         IF ( BAL.EQ.0 ) THEN
            ALLOCATE ( DK( N ) )
            CALL DCOPY( N, %VAL( SCAL ), 1, DK, 1 )
C
            DO 20 I = 1, N
               DK(I) = ONE/DK(I)
   20       CONTINUE
C
            CALL MB01SD( 'Column', M, N, %VAL( F ), LDF, DK, DK )
            DEALLOCATE ( DK )
         END IF
      ELSE
         INFO = 0
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO1.EQ.0 .AND. INFO.EQ.0 ) THEN
         IP = 1
         IF ( ISTAIR.GE.0 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( M, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( F ), mxGetPr( PLHS(IP) ),
     $                             M*N )
            IP = IP + 1
         END IF
         IF ( NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS(IP) ),
     $                             N*N )
            IP = IP + 1
         END IF
         IF ( NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( B ), mxGetPr( PLHS(IP) ),
     $                             N*M )
            IP = IP + 1
         END IF
         IF ( NLHS.GE.IP .AND. ( ISTAIR.EQ.-1 .OR. ISTAIR.EQ.0 ) ) THEN
            ALLOCATE ( DK( INDCON ) )
C
            DO 30 I = 1, INDCON
               DK(I) = KSTAIR(I)
   30       CONTINUE
C
            PLHS(IP) = mxCreateDoubleMatrix( INDCON,
     $                                       MIN( INDCON, 1 ), 0 )
            CALL mxCopyReal8ToPtr( DK, mxGetPr( PLHS(IP) ),
     $                             INDCON*MIN( INDCON, 1 ) )
            DEALLOCATE ( DK )
            IP = IP + 1
         END IF
         IF ( NLHS.GE.IP .AND. WITHU.EQ.1 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( U ), mxGetPr( PLHS(IP) ),
     $                             N*N )
            IP = IP + 1
         END IF
         IF ( NLHS.GE.IP .AND. WITHV.EQ.1 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( V ), mxGetPr( PLHS(IP) ),
     $                             M*M )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. ( ISTAIR.EQ.-1 .OR. ISTAIR.EQ.0 )
     $                   .AND. BAL.EQ.0 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, MIN( N, 1 ), 0 )
            CALL mxCopyReal8ToPtr( %VAL( SCAL ),
     $                             mxGetPr( PLHS(IP) ),
     $                             N*MIN( N, 1 ) )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE ( KSTAIR )
      CALL mxFree( A     )
      CALL mxFree( B     )
      CALL mxFree( DWORK )
      CALL mxFree( U     )
      IF ( ISTAIR.GE.-1 .AND. BAL.EQ.0 )
     $   CALL mxFree( SCAL )
      IF ( ISTAIR.GE.0 )
     $   CALL mxFree( F )
      IF ( ISTAIR.NE.2 )
     $   CALL mxFree( V )
      IF ( ISTAIR.LE.1 )
     $   DEALLOCATE ( IWORK )
C
C Error and warning handling.
C
      IF ( INFO1.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB01OD'')'
     $        ) INFO1
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB06ND'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of DEADBEAT ***
      END
