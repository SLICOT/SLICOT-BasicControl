C LDSIMT.F - Gateway function for computing the output response of a
C            linear discrete-time system using SLICOT routines TF01MD
C            and TF01ND.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [Y(,x)] = ldsimt(A,B,C,D,U(,x,HessA))
C
C Purpose:
C   To compute the output vector sequence y(1), y(2),..., y(t)
C
C        x(k+1) = A x(k) + B u(k)
C        y(k)   = C x(k) + D u(k),
C
C   given an initial state vector x(1), and the input vector sequence
C   u(1), u(2),..., u(t), where y(k) and u(k) are vectors of length p
C   and m, respectively. The input trajectories are given as
C
C        U = ( u(1)  u(2) ...  u(t) )
C
C   and the output trajectories result in similarly. Optionally, the
C   matrix A can be given in upper or lower Hessenberg form.
C
C Input parameters:
C   A      - the n-by-n state matrix A. If HessA = 1 or 2 only the upper
C            or lower Hessenberg part, respectively, must be defined.
C   B      - the n-by-m input matrix B.
C   C      - the p-by-n output matrix C.
C   D      - the p-by-m input-output matrix D.
C   U      - the m-by-t matrix U.
C   x      - (optional) the initial state x(1).
C            Default: x = 0.
C   HessA  - (optional) scalar indicating whether the matrix A is
C            general or in an upper/lower Hessenberg form:
C            = 0 :  general matrix;
C            = 1 :  upper Hessenberg matrix;
C            = 2 :  lower Hessenberg matrix.
C            Default: HessA = 0.
C
C Output parameters:
C   Y      - the p-by-t output matrix Y.
C   x      - (optional) the final state x(t+1).
C
C Comments
C   This gateway function stores the input and output trajectories in
C   a transposed form comparing to the LDSIM function. Moreover, A could
C   be given in a Hessenberg form.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003.
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
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
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
      CHARACTER         UPLO
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDU, LDY, M, N, NY, P
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), D(:,:),
     $                                 DWORK(:), U(:,:), X(:), Y(:,:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           HESSA, LDWORK
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          DCOPY, TF01MD, TF01ND
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'LDSIMT requires at least 5 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'LDSIMT requires at least 1 output argument' )
      END IF
C
C   A(nxn), B(nxm), C(pxn), D(pxm), and U(mxt).
C
      N  = mxGetM( PRHS(1) )
      M  = mxGetN( PRHS(2) )
      P  = mxGetM( PRHS(3) )
      NY = mxGetN( PRHS(5) )
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
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
      IF ( mxGetN( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $           ( 'C must have the same number of columns as A' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a real matrix' )
      END IF
      IF ( mxGetM( PRHS(4) ).NE.P ) THEN
         CALL mexErrMsgTxt( 'D must have the same number of rows as C' )
      END IF
      IF ( mxGetN( PRHS(4) ).NE.M ) THEN
         CALL mexErrMsgTxt
     $           ( 'D must have the same number of columns as B' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'U must be a real matrix' )
      END IF
      IF ( mxGetM( PRHS(5) ).NE.M ) THEN
         CALL mexErrMsgTxt( 'The number of rows of U must equal the numb
     $er of columns of B' )
      END IF
C
      IF ( NRHS.GE.6 ) THEN
C
C   x(1)
C
         IF ( mxGetM( PRHS(6) )*mxGetN( PRHS(6) ).LT.N ) THEN
            WRITE( TEXT, '(''X must have '',I7,'' entries'')' ) N
            CALL mexErrMsgTxt( TEXT )
         END IF
      END IF
C
      IF ( NRHS.GE.7 ) THEN
C
C   HessA
C
         IF ( mxGetM( PRHS(7) ).NE.1 .OR.
     $        mxGetN( PRHS(7) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'HESSA must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(7) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'HESSA must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), TEMP, 1 )
         HESSA = TEMP
         IF ( HESSA.LT.0 .OR. HESSA.GT.2 ) THEN
            CALL mexErrMsgTxt
     $        ( 'HESSA has 0, 1, or 2 the only admissible values' )
         END IF
      ELSE
         HESSA = 0
      END IF
C
      IF ( HESSA.EQ.1 ) THEN
         UPLO = 'U'
      ELSEIF ( HESSA.EQ.2 ) THEN
         UPLO = 'L'
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA = MAX( 1, N )
      LDB = LDA
      LDC = MAX( 1, P )
      LDD = LDC
      LDU = MAX( 1, M )
      LDY = LDC
      LDWORK = N
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( A( LDA, N ), B( LDB, M ), C( LDC, N ), D( LDD, M ),
     $           DWORK( LDWORK ), U( LDU, NY ), X( N ), Y( LDY, NY ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), C, P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), D, P*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), U, M*NY )
      IF ( NRHS.GE.6 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), X, N )
      ELSEIF ( N.GT.0 ) THEN
         X(1) = ZERO
         CALL DCOPY( N, X(1), 0, X(1), 1 )
      END IF
C
C Do the actual computations.
C
      IF ( HESSA.EQ.0 ) THEN
         CALL TF01MD( N, M, P, NY, A, LDA, B, LDB, C, LDC, D, LDD,
     $                U, LDU, X, Y, LDY, DWORK, INFO )
      ELSE
         CALL TF01ND( UPLO, N, M, P, NY, A, LDA, B, LDB, C, LDC, D,
     $                LDD, U, LDU, X, Y, LDY, DWORK, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( P, NY, 0 )
         CALL mxCopyReal8ToPtr( Y, mxGetPr( PLHS(1) ), P*NY )
         IF ( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( N, MIN( N, 1 ), 0 )
            CALL mxCopyReal8ToPtr( X, mxGetPr( PLHS(2) ), N )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( A, B, C, D, DWORK, U, X, Y )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( HESSA.EQ.0 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM TF01MD'')'
     $           ) INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM TF01ND'')'
     $           ) INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of LDSIMT ***
      END
