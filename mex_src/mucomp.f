C MUCOMP.F - Gateway function for computing the structured
C            singular value using the SLICOT routine AB13MD.
C
C RELEASE 2.0 of SLICOT Robust Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [BOUND, D, G] = mucomp(Z, K, T)
C      [BOUND, D] = mucomp(Z, K, T)
C           BOUND = mucomp(Z, K, T)
C
C Purpose:
C   To compute an upper bound on the structured singular value for a
C   given square complex matrix and given block structure of the
C   uncertainty.
C
C Input parameters:
C   Z      - the complex n-by-n matrix for which the structured
C            singular value is to be computed.
C   K      - the vector of length m containing the block structure
C            of the uncertainty; K(i) is the size of block i, i = 1:m.
C   T      - the vector of length m indicating the type of each block.
C            T(i) = 1 if the corresponding block is real;
C            T(i) = 2 if the corresponding block is complex.
C
C Output parameters:
C   BOUND  - the upper bound on the structured singular value.
C   D, G   - vectors of length n containing the diagonal entries
C            of the diagonal matrices D and G, respectively, such that
C            the matrix Z'*D^2*Z + sqrt(-1)*(G*Z-Z'*G) - bound^2*D^2
C            is negative semidefinite.
C
C Comments:
C   Currently, there are two limitations:
C   1. The size of a real block should be 1 (K(i) = 1 if T(i) = 1).
C   2. The sum of block sizes must be equal to n.
C
C References
C   [1] Fan, M.K.H., Tits, A.L., and Doyle, J.C.
C       Robustness in the presence of mixed parametric uncertainty and
C       unmodeled dynamics.
C       IEEE Trans. Automatic Control, vol. AC-36, 1991, pp. 25-38.
C
C Contributor:
C   P.Hr. Petkov, TU Sofia, Bulgaria, Oct. 2000.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2005,
C   Apr. 2009, Dec. 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C     .. Mex-file interface parameters ..
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C     .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mxGetPi, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      INTEGER           INFO, LDWORK, LDZ, LZWORK, M, N
      DOUBLE PRECISION  BOUND
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: ITYPE(:), IWORK(:), NBLOCK(:)
      DOUBLE PRECISION, ALLOCATABLE :: D(:), DWORK(:), G(:), X(:)
      COMPLEX*16,       ALLOCATABLE :: Z(:,:), ZWORK(:)
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           I, L, LIWORK
C
C     .. External subroutines ..
      EXTERNAL          AB13MD
C
C     ..Intrinsic functions..
      INTRINSIC         INT, MAX
C
C Check for proper number of arguments.
C
      IF( NRHS.NE.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'MUCOMP requires 3 input arguments' )
      ELSE IF ( NLHS.GT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'MUCOMP requires at most 3 output arguments' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   Z(nxn), K(m), T(m)
C
      N = mxGetM( PRHS(1) )
C
      IF( mxGetN( PRHS(1) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'Z must be a square matrix' )
      END IF
      IF( mxIsNumeric( PRHS(1) ).EQ.0 ) THEN
         CALL mexErrMsgTxt( 'Z must be a numeric matrix' )
      END IF
C
      M = MAX( mxGetM( PRHS(2) ), mxGetN( PRHS(2) ) )
      IF( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'K must be a real vector' )
      END IF
C
      L = MAX( mxGetM( PRHS(3) ), mxGetN( PRHS(3) ) )
      IF( L.NE.M ) THEN
         CALL mexErrMsgTxt
     $       ( 'T must have the same length as K' )
      END IF
      IF( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'T must be a real vector' )
      END IF
C
C Determine the lenghts of working arrays.
C
      LDZ    = MAX( 1, N )
      LIWORK = MAX( 4*M - 2, N )
      LDWORK = 2*N*N*M - N*N + 9*M*M + N*M + 11*N + 33*M - 11
      LZWORK = 6*N*N*M + 12*N*N + 6*M + 6*N - 3
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( D( N ), DWORK( LDWORK ), G( N ), ITYPE( M ),
     $           IWORK( LIWORK ), NBLOCK( M ), X( 2*M ), Z( LDZ, N ),
     $           ZWORK( LZWORK ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToComplex16( mxGetPr( PRHS(1) ), mxgetPi( PRHS(1) ),
     $                           Z, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), X, M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), X( M+1 ), M )
C
      DO 10 I = 1, M
         NBLOCK( I ) = INT( X( I ) )
         ITYPE( I )  = INT( X( M+I ) )
   10 CONTINUE
C
C Do the actual computations.
C
      CALL AB13MD( 'No information', N, Z, LDZ, M, NBLOCK, ITYPE, X,
     $             BOUND, D, G, IWORK, DWORK, LDWORK, ZWORK, LZWORK,
     $             INFO )
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( 1, 1, 0 )
      CALL mxCopyReal8ToPtr( BOUND, mxGetPr( PLHS(1) ), 1 )
      IF( NLHS.GT.1 ) THEN
         PLHS(2) = mxCreateDoubleMatrix( N, 1, 0 )
         CALL mxCopyReal8ToPtr( D, mxGetPr( PLHS(2) ), N )
         IF( NLHS.GT.2 ) THEN
            PLHS(3) = mxCreateDoubleMatrix( N, 1, 0 )
            CALL mxCopyReal8ToPtr( G, mxGetPr( PLHS(3) ), N )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( D, DWORK, G, ITYPE, IWORK, NBLOCK, X, Z, ZWORK )
C
C Error and warning handling.
C
      IF( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '', I4, '' ON EXIT FROM AB13MD'')'
     $          ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C
C *** Last line of MUCOMP ***
      END
