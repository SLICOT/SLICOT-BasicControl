C DLSDZ.F - Gateway function for Loop Shaping Design of
C           discrete-time systems using SLICOT routine SB10ZD.
C
C RELEASE 2.0 of SLICOT Robust Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [AK,BK,CK,DK,(RCOND)] = dlsdz(A,B,C,D,factor(,tol))
C
C Purpose:
C   To compute the matrices of the positive feedback controller
C
C              | Ak | Bk |
C          K = |----|----|
C              | Ck | Dk |
C
C   for the discrete-time shaped plant
C
C              | A | B |
C          G = |---|---|
C              | C | D |
C
C   in the McFarlane/Glover Loop Shaping Design Procedure.
C
C Input parameters:
C   A      - the n-by-n system state matrix A.
C   B      - the n-by-m system input matrix B.
C   C      - the p-by-n system output matrix C.
C   D      - the p-by-m system matrix D.
C   factor - = 1 implies that an optimal controller is required
C                (not reccomended);
C            > 1 implies that a suboptimal controller is required
C                achieving a performance FACTOR less than optimal.
C   tol    - (optional) tolerance used for checking the nonsingularity
C            of the matrices to be inverted.
C            Default:  tol = sqrt(epsilon_machine), where
C            epsilon_machine is the relative machine precision.
C
C Output parameters:
C   AK    - the n-by-n controller state matrix Ak.
C   BK    - the n-by-p controller input matrix Bk.
C   CK    - the m-by-n controller output matrix Ck.
C   DK    - the m-by-p controller matrix Dk.
C   RCOND - (optional) a vector containing estimates of the reciprocal
C           condition numbers of the matrices which have to be inverted
C           during the computation of the controller.
C           RCOND(1) contains an estimate of the reciprocal condition
C                    number of the linear system of equations from
C                    which the solution of the P-Riccati equation is
C                    obtained;
C           RCOND(2) contains an estimate of the reciprocal condition
C                    number of the linear system of equations from
C                    which the solution of the Q-Riccati equation is
C                    obtained;
C           RCOND(3) contains an estimate of the reciprocal condition
C                    number of the matrix (gamma^2-1)*In - P*Q;
C           RCOND(4) contains an estimate of the reciprocal condition
C                    number of the matrix Rx + Bx'*X*Bx;
C           RCOND(5) contains an estimate of the reciprocal condition
C                                                ^
C                    number of the matrix Ip + D*Dk;
C           RCOND(6) contains an estimate of the reciprocal condition
C                                              ^
C                    number of the matrix Im + Dk*D.
C
C References
C   [1] Gu, D.-W., Petkov, P.H., and Konstantinov, M.M.
C       On discrete H-infinity loop shaping design procedure routines.
C       Technical Report 00-6, Dept. of Engineering, Univ. of
C       Leicester, UK, 2000.
C
C Contributor:
C   P.Hr. Petkov, TU Sofia, Bulgaria, May 2001.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, July 2001,
C   Apr. 2009, Dec. 2012.
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
      INTEGER           INFO, LDA, LDAK, LDB, LDBK, LDC, LDCK, LDD,
     $                  LDDK, LDWORK, M, N, P
      DOUBLE PRECISION  FACTOR, TOL
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
      INTEGER           LBWORK, LIWORK
      DOUBLE PRECISION  RCOND(6)
C
C .. External subroutines ..
      EXTERNAL          SB10ZD
C
C ..Intrinsic functions..
      INTRINSIC         MAX
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'DLSDZ requires at least 5 input arguments' )
      ELSE IF ( NLHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'DLSDZ requires at least 4 output arguments' )
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
C
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
C     tol
C
      IF ( NRHS.GT.5 ) THEN
         IF ( mxGetM( PRHS(6) ).NE.1 .OR. mxGetN( PRHS(6) ).NE.1 ) THEN
             CALL mexErrMsgTxt( 'TOL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(6) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'TOL must be a numeric scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), TOL, 1 )
         IF ( TOL.GE.ONE ) THEN
             CALL mexErrMsgTxt( 'TOL must be less than 1' )
         END IF
      ELSE
         TOL = ZERO
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms
C in DGEES.
C
      LDA  = MAX( 1, N )
      LDB  = LDA
      LDC  = MAX( 1, P )
      LDD  = LDC
      LDAK = LDA
      LDBK = LDA
      LDCK = MAX( 1, M )
      LDDK = LDCK
C
      LBWORK = 2*N
      LIWORK = 2*MAX( N, P + M )
      LDWORK = 16*N*N + 5*M*M + 7*P*P + 6*M*N + 7*M*P + 7*N*P + 6*N +
     $         2*( M + P ) + MAX( 14*N + 23, 16*N, 2*M - 1, 2*P - 1 )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( A( LDA, N ), B( LDB, M ), C( LDC, N ), D( LDD, M ),
     $           AK( LDAK, N ), BK( LDBK, P ), CK( LDCK, N ),
     $           DK( LDDK, P ) )
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
      CALL SB10ZD( N, M, P, A, LDA, B, LDB, C, LDC, D, LDD, FACTOR, AK,
     $             LDAK, BK, LDBK, CK, LDCK, DK, LDDK, RCOND, TOL,
     $             IWORK, DWORK, LDWORK, BWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
      CALL mxCopyReal8ToPtr( AK, mxGetPr( PLHS(1) ), N*N )
      PLHS(2) = mxCreateDoubleMatrix( N, P, 0 )
      CALL mxCopyReal8ToPtr( BK, mxGetPr( PLHS(2) ), N*P )
      PLHS(3) = mxCreateDoubleMatrix( M, N, 0 )
      CALL mxCopyReal8ToPtr( CK, mxGetPr( PLHS(3) ), M*N )
      PLHS(4) = mxCreateDoubleMatrix( M, P, 0 )
      CALL mxCopyReal8ToPtr( DK, mxGetPr( PLHS(4) ), M*P )
      IF ( NLHS.GT.4 ) THEN
         PLHS(5) = mxCreateDoubleMatrix( 6, 1, 0 )
         CALL mxCopyReal8ToPtr( RCOND, mxGetPr( PLHS(5) ), 6 )
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
          WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB10ZD'')'
     $         ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of DLSDZ ***
      END
