C CONHIN.F - Gateway function for H_infinity or H_2 design of
C            continuous-time systems using SLICOT routines
C            SB10FD and SB10HD.
C
C RELEASE 2.0 of SLICOT Robust Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [AK,BK,CK,DK,(RCOND)] = conhin(task,A,B,C,D,ncon,nmeas,(gamma))
C
C   task = 1 :  [AK,BK,CK,DK,(RCOND)] = conhin(1,A,B,C,D,ncon,nmeas,
C                                              gamma)
C   task = 2 :  [AK,BK,CK,DK,(RCOND)] = conhin(2,A,B,C,D,ncon,nmeas)
C
C Purpose:
C
C   task = 1:
C
C     To compute the matrices of an H-infinity (sub)optimal n-state
C     controller
C
C              | AK | BK |
C          K = |----|----|,
C              | CK | DK |
C
C     for the continuous-time system
C
C              | A  | B1  B2  |   | A | B |
C          P = |----|---------| = |---|---|
C              | C1 | D11 D12 |   | C | D |
C              | C2 | D21 D22 |
C
C     and for a given value of gamma, where B2 has column size of the
C     number of control inputs (ncon) and C2 has row size of the number
C     of measurements (nmeas) being provided to the controller.
C
C     It is assumed that
C
C     (A1) (A,B2) is stabilizable and (C2,A) is detectable,
C
C     (A2) D12 is full column rank and D21 is full row rank,
C
C     (A3) | A-j*omega*I  B2  | has full column rank for all omega,
C          |    C1        D12 |
C
C     (A4) | A-j*omega*I  B1  |  has full row rank for all omega.
C          |    C2        D21 |
C
C   task = 2:
C
C     To compute the matrices of the H2 optimal n-state controller
C
C              | AK | BK |
C          K = |----|----|
C              | CK | DK |
C
C     for the continuous-time system
C
C              | A  | B1  B2  |   | A | B |
C          P = |----|---------| = |---|---|
C              | C1 |  0  D12 |   | C | D |
C              | C2 | D21 D22 |
C
C     where B2 has column size of the number of control inputs (ncon)
C     and C2 has row size of the number of measurements (nmeas) being
C     provided to the controller.
C
C     It is assumed that
C
C     (A1) (A,B2) is stabilizable and (C2,A) is detectable,
C
C     (A2) The block D11 of D is zero,
C
C     (A3) D12 is full column rank and D21 is full row rank.
C
C Input parameters:
C   task  - integer option to determine the type of the design:
C           = 1 : H_infinity design;
C           = 2 : H_2 design;
C   A     - the n-by-n system state matrix A.
C   B     - the n-by-m system input matrix B.
C   C     - the p-by-n system output matrix C.
C   D     - the p-by-m system matrix D.
C   ncon  - the number of control inputs. m >= ncon >= 0,
C           p-nmeas >= ncon.
C   nmeas - the number of measurements. p >= nmeas >= 0,
C           m-ncon >= nmeas.
C   gamma - (task 1 only) the parameter gamma used in H_infinity design.
C           It is assumed that gamma is sufficiently large so that the
C           controller is admissible. gamma >= 0.
C
C Output parameters:
C   AK    - the n-by-n controller state matrix AK.
C   BK    - the n-by-nmeas controller input matrix BK.
C   CK    - the ncon-by-n controller output matrix CK.
C   DK    - the ncon-by-nmeas controller matrix DK.
C   RCOND - (optional) a vector containing estimates of the reciprocal
C           condition numbers of the matrices which are to be inverted
C           and estimates of the reciprocal condition numbers of the
C           Riccati equations which have to be solved during the
C           computation of the controller. (See the description of
C           the algorithm in [1].)
C           RCOND(1) contains the reciprocal condition number of the
C                    control transformation matrix TU,
C           RCOND(2) contains the reciprocal condition number of the
C                    measurement transformation matrix TY,
C           RCOND(3) contains an estimate of the reciprocal condition
C                    number of the X-Riccati equation,
C           RCOND(4) contains an estimate of the reciprocal condition
C                    number of the Y-Riccati equation.
C
C References
C   [1] P.Hr. Petkov, D.W. Gu and M.M. Konstantinov. Fortran 77 routines
C       for Hinf and H2 design of continuous-time linear control systems.
C       Report98-14, Department of Engineering, Leicester University,
C       August 1998.
C
C Contributor:
C   P.Hr. Petkov, TU Sofia, Bulgaria, Sep. 1999.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009,
C   Dec. 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C
C .. Mex-file interface parameters ..
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDAK, LDBK, LDCK,
     $                  LDDK, LDWORK, M, N, NCON, NMEAS, P
      DOUBLE PRECISION  GAMMA, TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      LOGICAL,          ALLOCATABLE :: BWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), D(:,:),
     $                  AK(:,:), BK(:,:), CK(:,:), DK(:,:), DWORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           LBWORK, LIWORK, Q, TASK
      DOUBLE PRECISION  TEMP
      DOUBLE PRECISION  RCOND(4)
C
C .. External subroutines ..
      EXTERNAL          SB10FD, SB10HD
C
C ..Intrinsic functions..
      INTRINSIC         MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.7 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CONHIN requires at least 7 input arguments' )
      ELSE IF ( NLHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CONHIN requires at least 4 output arguments' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   task
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'TASK must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'TASK must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      TASK = TEMP
      IF ( TASK.LT.1 .OR. TASK.GT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'TASK has 1 or 2 the only admissible values' )
      END IF
C
      IF ( TASK.EQ.1 ) THEN
         IF ( NRHS.NE.8 ) THEN
            CALL mexErrMsgTxt
     $          ( 'CONHIN requires 8 input arguments' )
         END IF
      END IF
C
C   A(nxn), B(nxm), C(pxn), D(pxm)
C
      N = mxGetM( PRHS(2) )
      M = mxGetN( PRHS(3) )
      P = mxGetM( PRHS(4) )
C
      IF ( M.LE.0 ) THEN
         CALL mexErrMsgTxt( 'The system has no inputs' )
      END IF
      IF ( P.LE.0 ) THEN
         CALL mexErrMsgTxt( 'The system has no outputs' )
      END IF
C
      IF ( mxGetN( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF ( mxGetM( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same row dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF ( mxGetN( PRHS(4) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $       ('C must have the same column dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
      IF ( mxGetM( PRHS(5) ).NE.P ) THEN
         CALL mexErrMsgTxt( 'D must have the same row dimension as C' )
      END IF
      IF ( mxGetN( PRHS(5) ).NE.M ) THEN
         CALL mexErrMsgTxt
     $       ( 'D must have the same column dimension as B' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a real matrix' )
      END IF
C
C     ncon
      IF ( mxGetM( PRHS(6) ).NE.1 .OR. mxGetN( PRHS(6) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'NCON must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(6) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'NCON must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), TEMP, 1 )
      NCON = TEMP
      IF ( NCON.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'NCON must be a positive integer' )
      END IF
C
C     nmeas
      IF ( mxGetM( PRHS(7) ).NE.1 .OR. mxGetN( PRHS(7) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'NMEAS must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(7) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'NMEAS must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), TEMP, 1 )
      NMEAS = TEMP
      IF ( NMEAS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'NMEAS must be a positive integer' )
      END IF
C
C     gamma
      IF ( TASK.EQ.1 ) THEN
         IF ( mxGetM( PRHS(8) ).NE.1 .OR. mxGetN( PRHS(8) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'GAMMA must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(8) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'GAMMA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(8) ), GAMMA, 1 )
         IF ( GAMMA.LE.ZERO ) THEN
            CALL mexErrMsgTxt( 'GAMMA must be a positive scalar' )
         END IF
      END IF
C
      TOL = -ONE
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms
C in DGEES.
C
      LDA  = MAX( 1, N )
      LDB  = MAX( 1, N )
      LDC  = MAX( 1, P )
      LDD  = MAX( 1, P )
      LDAK = MAX( 1, N )
      LDBK = MAX( 1, N )
      LDCK = MAX( 1, NCON )
      LDDK = MAX( 1, NCON )
      LBWORK = 2*N
      Q = MAX( M - NCON, NCON, P - NMEAS, NMEAS )
      IF ( TASK.EQ.1 ) THEN
         LIWORK = MAX ( 2*MAX( N, M-NCON, P-NMEAS, NCON ), N*N )
         LDWORK = 2*Q*( 3*Q + 2*N ) +
     $            MAX( 1, ( N + Q )*( N + Q + 6 ),
     $                 Q*( Q + MAX( N, Q, 5 ) + 1 ), 2*N*( N + 2*Q ) +
     $                 MAX( 1, 4*Q*Q +
     $                      MAX( 2*Q, 3*N*N +
     $                           MAX( 2*N*Q, 10*N*N + 12*N + 5 ) ),
     $                           Q*( 3*N + 3*Q +
     $                               MAX( 2*N, 4*Q + max( N, Q ) ) ) ) )
      ELSE
         LIWORK = MAX ( 2*N, N*N )
         LDWORK = 2*Q*( 3*Q + 2*N ) +
     $            MAX( 1, Q*( Q + MAX( N, 5 ) + 1 ),
     $                 N*( 14*N + 12 + 2*Q ) + 5 )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( A( LDA, N ), B( LDB, M ), C( LDC, N ), D( LDD, M ),
     $           AK( LDAK, N ), BK( LDBK, NMEAS ), CK( LDCK, N ),
     $           DK( LDDK, NMEAS ) )
      ALLOCATE ( BWORK( LBWORK ), IWORK( LIWORK ), DWORK( LDWORK ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), C, P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), D, P*M )
C
C Do the actual computations.
C
      IF ( TASK.EQ.1 ) THEN
         CALL SB10FD( N, M, P, NCON, NMEAS, GAMMA, A, LDA, B, LDB,
     $                C, LDC, D, LDD, AK, LDAK, BK, LDBK, CK, LDCK,
     $                DK, LDDK, RCOND, TOL, IWORK, DWORK, LDWORK,
     $                BWORK, INFO )
      ELSE
         CALL SB10HD( N, M, P, NCON, NMEAS, A, LDA, B, LDB, C, LDC,
     $                D, LDD, AK, LDAK, BK, LDBK, CK, LDCK, DK,
     $                LDDK, RCOND, TOL, IWORK, DWORK, LDWORK, BWORK,
     $                INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
      CALL mxCopyReal8ToPtr( AK, mxGetPr( PLHS(1) ), N*N )
      PLHS(2) = mxCreateDoubleMatrix( N, NMEAS, 0 )
      CALL mxCopyReal8ToPtr( BK, mxGetPr( PLHS(2) ), N*NMEAS )
      PLHS(3) = mxCreateDoubleMatrix( NCON, N, 0 )
      CALL mxCopyReal8ToPtr( CK, mxGetPr( PLHS(3) ), NCON*N )
      PLHS(4) = mxCreateDoubleMatrix( NCON, NMEAS, 0 )
      CALL mxCopyReal8ToPtr( DK, mxGetPr( PLHS(4) ), NCON*NMEAS )
      IF ( NLHS.GT.4 ) THEN
         PLHS(5) = mxCreateDoubleMatrix( 4, 1, 0 )
         CALL mxCopyReal8ToPtr( RCOND, mxGetPr( PLHS(5) ), 4 )
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( A, B, C, D, AK, BK, CK, DK, BWORK, IWORK, DWORK )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( TASK.EQ.1 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB10FD'')'
     $           ) INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB10HD'')'
     $           ) INFO
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of CONHIN ***
      END
