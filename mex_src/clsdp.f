C CLSDP.F - Gateway function for Loop Shaping Design of
C           continuous-time systems using SLICOT routine
C           SB10ID.
C
C RELEASE 2.0 of SLICOT Robust Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [AK,BK,CK,DK,(RCOND)] = clsdp(A,B,C,D,factor)
C
C Purpose:
C     To compute the matrices of the positive feedback controller
C
C              | Ak | Bk |
C          K = |----|----|
C              | Ck | Dk |
C
C     for the shaped plant
C
C              | A | B |
C          G = |---|---|
C              | C | D |
C
C     in the McFarlane/Glover Loop Shaping Design Procedure.
C
C Input parameters:
C   A      - the n-by-n system state matrix A.
C   B      - the n-by-m system input matrix B.
C   C      - the p-by-n system output matrix C.
C   D      - the p-by-m system matrix D.
C   factor - = 1 implies that an optimal controller is required
C            > 1 implies that a suboptimal controller is required
C                achieving a performance FACTOR less than optimal.
C
C Output parameters:
C   AK    - the nk-by-nk controller state matrix Ak.
C   BK    - the nk-by-np controller input matrix Bk.
C   CK    - the m-by-nk controller output matrix Ck.
C   DK    - the m-by-np controller matrix Dk.
C   RCOND - (optional) a vector containing estimates of the reciprocal
C           condition numbers of the Riccati equations which have to be
C           solved during the computation of the controller.
C           RCOND(1) contains an estimate of the reciprocal condition
C                    number of the X-Riccati equation,
C           RCOND(2) contains an estimate of the reciprocal condition
C                    number of the Z-Riccati equation.
C
C References
C   [1] D.W. Gu, P.Hr. Petkov, D.W. Gu and M.M. Konstantinov.
C       Hinf Loop Shaping Design procedure routines in SLICOT.
C       NICONET Report 1999-15, November 1999.
C
C Contributor:
C   P.Hr. Petkov, TU Sofia, Bulgaria, Oct. 2000.
C
C Revisions:
C   V. Sima, February 2001, Apr. 2009, Dec. 2012.
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
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDAK, LDBK, LDCK,
     $                  LDDK, LDWORK, M, N, NK, P
      DOUBLE PRECISION  FACTOR
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      LOGICAL,          ALLOCATABLE :: BWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), D(:,:),
     $                  AK(:,:), BK(:,:), CK(:,:), DK(:,:),
     $                  AKA(:,:), BKA(:,:), CKA(:,:), DKA(:,:), DWORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           LBWORK, LIWORK
      DOUBLE PRECISION  RCOND(2)
C
C .. External subroutines ..
      EXTERNAL          SB10ID
C
C ..Intrinsic functions..
      INTRINSIC         MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.NE.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CLSDP requires 5 input arguments' )
      ELSE IF ( NLHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CLSDP requires at least 4 output arguments' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   A(nxn), B(nxm), C(pxn), D(pxm)
C
      N = mxGetM( PRHS(1) )
      M = mxGetN( PRHS(2) )
      P = mxGetM( PRHS(3) )
C
      IF ( mxGetN( PRHS(1) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a numeric matrix' )
      END IF
      IF ( mxGetM( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same row dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a numeric matrix' )
      END IF
      IF ( mxGetN( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $       ('C must have the same column dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a numeric matrix' )
      END IF
      IF ( mxGetM( PRHS(4) ).NE.P ) THEN
         CALL mexErrMsgTxt( 'D must have the same row dimension as C' )
      END IF
      IF ( mxGetN( PRHS(4) ).NE.M ) THEN
         CALL mexErrMsgTxt
     $       ( 'D must have the same column dimension as B' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a numeric matrix' )
      END IF
      IF ( M.LE.0 ) THEN
         CALL mexErrMsgTxt( 'The system has no inputs' )
      END IF
      IF ( P.LE.0 ) THEN
         CALL mexErrMsgTxt( 'The system has no outputs' )
      END IF
C
C     factor
      IF ( mxGetM( PRHS(5) ).NE.1 .OR. mxGetN( PRHS(5) ).NE.1 ) THEN
          CALL mexErrMsgTxt( 'FACTOR must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
           CALL mexErrMsgTxt( 'FACTOR must be a numeric scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), FACTOR, 1 )
      IF ( FACTOR.LT.ONE ) THEN
          CALL mexErrMsgTxt( 'FACTOR must be equal or greater than 1' )
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms
C in DGEES.
C
      LDA = MAX( 1, N )
      LDB = MAX( 1, N )
      LDC = MAX( 1, P )
      LDD = MAX( 1, P )
      LDAK = MAX( 1, N )
      LDBK = MAX( 1, N )
      LDCK = MAX( 1, M )
      LDDK = MAX( 1, M )
      LBWORK = 2*N
      LIWORK = MAX ( 2*N, N*N, M, P )
      LDWORK = 4*N*N + M*M + P*P + 2*M*N + N*P + 4*N +
     $         MAX( 6*N*N + 5 + MAX( 1, 4*N*N + 8*N ), N*P + 2*N )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( A( LDA, N ), B( LDB, M ), C( LDC, N ), D( LDD, M ),
     $           AK( LDAK, N ), BK( LDBK, P ), CK( LDCK, N ),
     $           DK( LDDK, P ), AKA( LDAK, N ), BKA( LDBK, P ),
     $           CKA( LDCK, N ), DKA( LDDK, P ) )
      ALLOCATE ( BWORK( LBWORK ), IWORK( LIWORK ), DWORK( LDWORK ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), C, P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), D, P*M )
C
C Do the actual computations.
C
      CALL SB10ID( N, M, P, A, LDA, B, LDB, C, LDC, D, LDD, FACTOR,
     $             NK, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK, RCOND,
     $             IWORK, DWORK, LDWORK, BWORK, INFO )
C
C Copy output which is stored in local array to matrix output
C
      CALL DLACPY( 'F', NK, NK, AK, LDAK, AKA, NK )
      CALL DLACPY( 'F', NK, P, BK, LDBK, BKA, NK )
      CALL DLACPY( 'F', M, NK, CK, LDCK, CKA, M )
      CALL DLACPY( 'F', M, P, DK, LDDK, DKA, M )
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( NK, NK, 0 )
      CALL mxCopyReal8ToPtr( AKA, mxGetPr( PLHS(1) ), NK*NK )
      PLHS(2) = mxCreateDoubleMatrix( NK, P, 0 )
      CALL mxCopyReal8ToPtr( BKA, mxGetPr( PLHS(2) ), NK*P )
      PLHS(3) = mxCreateDoubleMatrix( M, NK, 0 )
      CALL mxCopyReal8ToPtr( CKA, mxGetPr( PLHS(3) ), M*NK )
      PLHS(4) = mxCreateDoubleMatrix( M, P, 0 )
      CALL mxCopyReal8ToPtr( DKA, mxGetPr( PLHS(4) ), M*P )
      IF ( NLHS.GT.4 ) THEN
         PLHS(5) = mxCreateDoubleMatrix( 2, 1, 0 )
         CALL mxCopyReal8ToPtr( RCOND, mxGetPr( PLHS(5) ), 2 )
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( A, B, C, D, AK, BK, CK, DK, AKA, BKA, CKA, DKA, BWORK,
     $            IWORK, DWORK )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
          WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB10ID'')'
     $         ) INFO
      END IF
C
      IF ( INFO.NE.0 ) THEN
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of CLSDP ***
      END
