#include "fintrf.h"
C SYSCOM.F - Gateway function for computing controllability,
C            observability, and minimality forms of a given system,
C            based on staircase form reductions using the
C            SLICOT routines AB07MD, TB01PD, TB01UD, and TB01ZD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [A,B,C(,N,U,sizes)] = syscom(task,A,B,C,tol(,bal))
C
C   task = 1 :  [A,B,C,Nc,U,sizes] = syscom(1,A,B,C,tol)
C   task = 2 :  [A,B,C,No,U,sizes] = syscom(2,A,B,C,tol)
C   task = 3 :          [Ar,Br,Cr] = syscom(3,A,B,C,tol,bal)
C
C Purpose:
C   To transform the matrix triple (A, B, C) to various staircase forms
C   showing controllability, observability or minimality.
C   The staircase forms are as follows:
C
C   a) Controllability (block) form
C
C                                |-  Nc -|
C                           ---[ * . . . * * . . * | * ]
C                           |  [ *       . .     . | 0 ]
C                           Nc [   .     . .     . | . ]
C                           |  [     .   . .     . | . ]
C       [ U'AU | U'B ]      ---[       * * * . . * | . ]
C       [------|-----] =       [           * . . * | . ];          (1)
C       [  CU  |  0  ]         [           .     . | . ]
C                              [           * . . * | 0 ]
C                              [-------------------|---]
C                              [ * . . . * * . . * | 0 ]
C
C   b) Observability (block) form
C
C                                |-  No -|
C                           ---[ * *               | * ]
C                           |  [ .   .             | * ]
C                           No [ .     .           | . ]
C                           |  [ .       *         | . ]
C       [ U'AU | U'B ]      ---[ * . . . *         | . ]
C       [------|-----] =       [ * . . . * * . . * | . ];          (2)
C       [  CU  |  0  ]         [ .       . .     . | . ]
C                              [ * . . . * * . . * | * ]
C                              [-------------------|---]
C                              [ * 0 . . . . . . 0 | 0 ]
C
C   c) Minimality (block) form
C
C                            |--Nm --|
C                        --[ * . . . *         * . . . * | * ]
C                        | [ *       .         .       . | . ]
C                        Nm[   .  Ar .         .       . | Br]
C                        | [     .   .         .       . | . ]
C                        --[       * *         * . . . * | * ]
C                          [ * . . . * * . . * * . . . * | * ]
C          [ Ar | Br ]     [ .       . .     . .       . | . ]
C          [----|----] =   [ .       . .     . .       . | . ].    (3)
C          [ Cr | 0  ]     [ * . . . * * . . * * . . . * | * ]
C                          [                   * . . . * | 0 ]
C                          [                   .       . | . ]
C                          [                   * . . . * | 0 ]
C                          [-----------------------------|---]
C                          [* . Cr . * 0 . . 0 * . . . * | 0 ]
C
C Input parameters:
C   task  - integer option to determine which form is computed:
C           = 1 : compute the controllability form (1);
C           = 2 : compute the observabilty form (2);
C           = 3 : compute the minimal realization subsystem (Ar,Br,Cr)
C                 from (3).
C   A     - real N-by-N state matrix.
C   B     - real N-by-M input matrix.
C   C     - real P-by-N output matrix.
C   tol   - (optional) real tolerance value used for rank decisions,
C           as a lower bound on the reciprocal condition numbers, or
C           for checking controllability or observability for
C           single-input or single-output systems, respectively.
C           Default:  tol = N*N*machine_epsilon, for rank decisions;
C                     tol = N*machine_epsilon*max(norm(A),norm(BC)), for
C                           checking controllability or observability,
C                           with BC = B or BC = C, respectively.
C   bal   - (for task = 3 only: optional) integer indicating whether
C           or not to balance the system triple before computing the
C           minimal realization subsystem.
C           bal = 0 : use balancing;
C                     otherwise, do not use balancing.
C           Default:  bal = 0.
C
C Output parameters:
C
C   For task = 1, 2,
C
C   A     - transformed N-by-N state matrix.
C   B     - transformed N-by-M input matrix.
C   C     - transformed P-by-N output matrix.
C   Nc    - (optional) order of the controllable subsystem when
C           task = 1.
C   No    - (optional) order of the observable subsystem when
C           task = 2.
C   U     - (optional) N-by-N orthogonal matrix carrying out the
C           transformations.
C   sizes - (optional) integer vector indicating the sizes of the
C           diagonal blocks of A corresponding to controllable or
C           observable parts, respectively. The length of sizes
C           indicates the associated controllable or observable index.
C
C   For task = 3,
C
C   Ar    - state matrix of the minimal realization subsystem.
C   Br    - input matrix of the minimal realization subsystem.
C   Cr    - output matrix of the minimal realization subsystem.
C
C Contributor:
C   H. Xu, TU Chemnitz, FR Germany, Dec. 1998.
C
C Revisions:
C   V. Sima, Katholieke Univ. Leuven, Belgium, May 1999.
C   V. Sima, Research Institute for Informatics, Bucharest, Jan. 2002,
C   Oct. 2004, Apr. 2009, Dec. 2012, Sep. 2015, Jul. 2016.
C
C **********************************************************************
C
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
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
      CHARACTER         JOBU, EQUIL
      INTEGER           INFO, LDA, LDB, LDC, LDU, LDWORK, M, MP, N,
     $                  NIND, NR, P, PM
      DOUBLE PRECISION  TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:), NBLK(:)
      DOUBLE PRECISION, ALLOCATABLE :: RNBLK(:)
      mwPointer         A, B, C, DWORK, TAU, U
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           BAL, ISIZE, K, KL, MXMP, TASK
      DOUBLE PRECISION  ANORM, BNORM, CNORM, NU, TEMP
C
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE
C
C .. External subroutines ..
      EXTERNAL          AB07MD, DLACPY, TB01PD, TB01UD, TB01ZD
C
C .. Intrinsic functions ..
      INTRINSIC         DBLE, MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSCOM requires at least 4 input arguments' )
      ELSE IF ( NLHS.GT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSCOM requires at most 6 output arguments' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   task
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 )
     $      THEN
         CALL mexErrMsgTxt( 'TASK must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
      CALL mexErrMsgTxt( 'TASK must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      TASK = TEMP
      IF ( TASK.LT.1 .OR. TASK.GT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'TASK has 1, 2, or 3 the only admissible values' )
      END IF
C
C   A(NxN), B(NxM), C(PxN)
C
      N = mxGetM( PRHS(2) )
      M = mxGetN( PRHS(3) )
      P = mxGetM( PRHS(4) )
C
      IF ( mxGetN( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same row dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
C
      IF ( mxGetN( PRHS(4) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'C must have the same column dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
C
      TOL = ZERO
      BAL = 0
      IF ( NRHS.GE.5 ) THEN
C
C   tol
C
         ISIZE = mxGetM( PRHS(5) ) * mxGetN( PRHS(5) )
         IF ( ISIZE.GT.1 )
     $      CALL mexErrMsgTxt( 'TOL must be a scalar' )
         IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(5) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), TEMP, ISIZE )
         IF ( ISIZE.GT.0 ) TOL = TEMP
C
         IF ( TASK.EQ.3 .AND. NRHS.GE.6 ) THEN
C
C   bal
C
            ISIZE = mxGetM( PRHS(6) ) * mxGetN( PRHS(6) )
            IF ( ISIZE.GT.1 )
     $         CALL mexErrMsgTxt( 'BAL must be a scalar' )
            IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(6) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'BAL must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), TEMP, ISIZE )
            IF ( ISIZE.GT.0 ) BAL = TEMP
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms:
C     LDWORK = MAX( 1, N*NB, 3*M, P*NB ), for TASK = 1,
C     LDWORK = MAX( 1, N*NB, M*NB, 3*P ), for TASK = 2,
C     LDWORK = MAX( 1, N + MAX( N, 3*M, 3*P )*NB ), for TASK = 3,
C where NB is the maximal block size for DGEHRD, DORGQR, DORMHR, and
C DORMQR.
C
      MXMP = MAX( M, P )
      IF ( TASK.EQ.1 ) THEN
         LDWORK = MAX( 1, N, 3*M, P )
         MP = M
         PM = P
      ELSE IF ( TASK.EQ.2 ) THEN
         LDWORK = MAX( 1, N, M, 3*P )
         MP = P
         PM = M
      ELSE
         LDWORK = MAX( 1, N + MAX( N, 3*MXMP) )
      END IF
C
      LDA = MAX( 1, N )
      LDB = LDA
      IF ( TASK.EQ.1 ) THEN
         LDC = MAX( 1, P )
      ELSE
         LDC = MAX( 1, MXMP )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      A = mxCalloc( LDA*N, 8 )
      C = mxCalloc( LDC*N, 8 )
      IF ( TASK.EQ.1 ) THEN
         B = mxCalloc( LDB*M, 8 )
      ELSE
         B = mxCalloc( LDB*MXMP, 8 )
      END IF
      DWORK = mxCalloc( LDWORK, 8 )
C
      IF ( TASK.NE.3 ) THEN
         IF ( NLHS.GT.4 ) THEN
            JOBU = 'I'
            LDU  = LDA
            U = mxCalloc( LDU*N, 8 )
         ELSE
            JOBU = 'N'
            LDU  = 1
            U = mxCalloc( 1, 8 )
         END IF
         TAU = mxCalloc( N, 8 )
         ALLOCATE( IWORK( MP ), NBLK( N ), RNBLK( N ) )
      ELSE
         ALLOCATE( IWORK( N + MXMP ) )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( A ), N*N )
      IF ( TASK.EQ.1 .OR. P.GE.M ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL( C ), P*N )
      ELSE
C
C        P < M. Use B temporarily for loading C.
C
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL( B ),
     $                          P*N )
         CALL DLACPY( 'Full', P, N, %VAL( B ), MAX( 1, P ), %VAL( C ),
     $                LDC )
C
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL( B ), N*M )
C
C Do the actual computations.
C
      IF ( TASK.EQ.2 ) THEN
C
         CALL AB07MD( 'Z', N, M, P, %VAL( A ), LDA, %VAL( B ), LDB,
     $                %VAL( C ), LDC, TEMP, 1, INFO )
C
         IF ( INFO.NE.0 )  THEN
            WRITE( TEXT, '('' INFO = '', I4, '' ON EXIT FROM AB07MD '')'
     $           ) INFO
            GO TO 40
         END IF
      END IF
C
      IF ( TASK.NE.3 ) THEN
         IF ( MP.EQ.1 ) THEN
            CALL TB01ZD( JOBU, N, PM, %VAL( A ), LDA, %VAL( B ),
     $                   %VAL( C ), LDC, NR, %VAL( U ), LDU,
     $                   %VAL( TAU ), TOL, %VAL( DWORK ), LDWORK, INFO )
         ELSE
            CALL TB01UD( JOBU, N, MP, PM, %VAL( A ), LDA, %VAL( B ),
     $                   LDB, %VAL( C ), LDC, NR, NIND, NBLK, %VAL( U ),
     $                   LDU, %VAL( TAU ), TOL, IWORK, %VAL( DWORK ),
     $                   LDWORK, INFO )
         END IF
C
         IF ( INFO.EQ.0 .AND. TASK.EQ.2 )
     $      CALL AB07MD( 'Z', N, P, M, %VAL( A ), LDA, %VAL( B ), LDB,
     $                   %VAL( C ), LDC, TEMP, 1, INFO )
C
      ELSE
C
         BNORM = DLANGE( 'F', N, M, %VAL( B ), LDB, %VAL( DWORK ) )
         CNORM = DLANGE( 'F', P, N, %VAL( C ), LDC, %VAL( DWORK ) )
C
         IF ( MIN( M, P ).EQ.1 .AND. TOL.EQ.ZERO ) THEN
C
C Compute the default tolerance for minimal realization of single-input
C or single-output systems.
C
            ANORM = DLANGE( 'F', N, N, %VAL( A ), LDA, %VAL( DWORK ) )
            IF ( M.EQ.1 )
     $         TOL = DBLE( N )*MAX( ANORM, BNORM )
            IF ( P.EQ.1 )
     $         TOL = MAX( TOL, DBLE( N )*MAX( ANORM, CNORM ) )
            TOL = TOL*DLAMCH( 'EPSILON' )
         END IF
C
         NR = N
         IF ( TOL.GE.ZERO ) THEN
C
C Check for nonminimal single-input or single-output systems.
C
            IF ( ( M.EQ.1 .AND. BNORM.LE.TOL ) .OR.
     $           ( P.EQ.1 .AND. CNORM.LE.TOL ) )
     $         NR = 0
         END IF
C
         IF ( NR.EQ.N ) THEN
            IF ( BAL.EQ.0 ) THEN
               EQUIL = 'S'
            ELSE
               EQUIL = 'N'
            END IF
            CALL TB01PD( 'M', EQUIL, N, M, P, %VAL( A ), LDA, %VAL( B ),
     $                   LDB, %VAL( C ), LDC, NR, TOL, IWORK,
     $                   %VAL( DWORK ), LDWORK, INFO )
         END IF
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         IF ( TASK.EQ.3 ) THEN
C
C           Redefine N.
C
            N = NR
            CALL DLACPY( 'Full', N, N, %VAL( A ), LDA, %VAL( A ),
     $                   MAX( 1, N ) )
            CALL DLACPY( 'Full', N, M, %VAL( B ), LDB, %VAL( B ),
     $                   MAX( 1, N ) )
         END IF
         IF ( NLHS.GE.1 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS(1) ),
     $                             N*N )
         END IF
         IF ( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( N, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( B ), mxGetPr( PLHS(2) ),
     $                             N*M )
         END IF
         IF ( NLHS.GE.3 ) THEN
            IF ( TASK.EQ.2 .OR. TASK.EQ.3 )
     $         CALL DLACPY( 'Full', P, N, %VAL( C ), LDC, %VAL( C ),
     $                      MAX( 1, P ) )
            PLHS(3) = mxCreateDoubleMatrix( P, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( C ), mxGetPr( PLHS(3) ),
     $                             P*N )
            IF ( TASK.NE.3 ) THEN
               IF ( NLHS.GE.4 ) THEN
                  PLHS(4) = mxCreateDoubleMatrix( 1, 1, 0 )
                  NU = NR
                  CALL mxCopyReal8ToPtr( NU, mxGetPr( PLHS(4) ),
     $                                   1 )
               END IF
               IF ( NLHS.GE.5 ) THEN
                  PLHS(5) = mxCreateDoubleMatrix( N, N, 0 )
                  CALL mxCopyReal8ToPtr( %VAL( U ),
     $                                   mxGetPr( PLHS(5) ), N*N )
               END IF
               IF ( NLHS.GE.6 ) THEN
                  IF ( MP.GT.1 ) THEN
                     KL = NIND
                     DO 20 K = 1, KL
                        RNBLK( K ) = NBLK( K )
   20                CONTINUE
                  ELSE
                     KL = NR
                     DO 30 K = 1, KL
                        RNBLK( K ) = ONE
   30                CONTINUE
                  END IF
                  PLHS(6) = mxCreateDoubleMatrix( KL, 1, 0 )
                  CALL mxCopyReal8ToPtr( RNBLK ,
     $                                   mxGetPr( PLHS(6) ), KL )
               END IF
            END IF
         END IF
      END IF
C
   40 CONTINUE
C
C Deallocate local arrays.
C !Fortran 90/95
C
      CALL mxFree( A     )
      CALL mxFree( B     )
      CALL mxFree( C     )
      CALL mxFree( DWORK )
      DEALLOCATE ( IWORK )
      IF ( TASK.NE.3 ) THEN
         CALL mxFree( TAU )
         CALL mxFree( U   )
         DEALLOCATE( NBLK, RNBLK )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( TASK.NE.3 ) THEN
            IF ( MP.EQ.1 ) THEN
               WRITE( TEXT,'('' INFO = '',I4,'' ON EXIT FROM TB01ZD'')'
     $              ) INFO
            ELSE
               WRITE( TEXT,'('' INFO = '',I4,'' ON EXIT FROM TB01UD'')'
     $              ) INFO
            END IF
         ELSE
            WRITE( TEXT,'('' INFO = '',I4,'' ON EXIT FROM TB01PD'')' )
     $         INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SYSCOM ***
      END
