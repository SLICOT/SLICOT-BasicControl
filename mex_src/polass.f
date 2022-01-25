C POLASS.F - Gateway function to perform (partial) pole assignment,
C            using SLICOT routine SB01BD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [F(,split,WRo,WIo,Z,Ao)] = polass(A,B,WR,WI(,tol,discr,alpha))
C
C Purpose:
C  To determine the state feedback matrix F for a given system (A,B)
C  such that the closed-loop state matrix A+B*F has specified
C  eigenvalues (closed-loop system poles).
C
C Input parameters:
C   A     - the n-by-n system state matrix A.
C   B     - the n-by-m system input matrix B.
C   WR    - the np-vector of real parts of the desired system poles,
C           np <= n.
C   WI    - the np-vector of imaginary parts of the desired system
C           poles. Complex conjugate pairs must appear consecutively.
C   tol   - (optional) absolute tolerance level below which the elements
C           of A or B are considered zero (used for controllability
C           tests). If tol <= 0, then a default tolerance is used.
C           Default:    n*epsilon_machine*max(norm(A),norm(B)), where
C                       epsilon_machine is the relative machine
C                       precision, and norm denotes the 1-norm.
C   discr - (optional) scalar indicating the type of system:
C           = 0 : continuous-time (default);
C           = 1 : discrete-time.
C   alpha - (optional) scalar specifying the maximum admissible value,
C           either for real parts, if discr = 0, or for moduli,
C           if discr = 1, of the eigenvalues of A which will not be
C           modified by the eigenvalue assignment algorithm
C           (alpha >= 0 if discr = 1).
C           Default:    -sqrt(epsilon_machine)  for continuous-time;
C                    1.0-sqrt(epsilon_machine)  for discrete-time.
C
C Output parameters:
C   F     - the m-by-n state feedback matrix, which assigns nap
C           closed-loop poles and keeps unaltered n-nap open-loop poles.
C   split - (optional) a 3-vector containing the pole splitting details;
C           split = [nfp; nap; nup], where
C                    nfp - the number of fixed poles, not modified by
C                          the assignment algorithm (poles having real
C                          parts, if discr = 0, or moduli, if discr = 1,
C                          less than alpha);
C                    nap - the number of assigned poles,
C                          nap = n - nfp - nup;
C                    nup - the number of uncontrollable poles detected
C                          by the algorithm.
C   WRo,  - (optional) the np-vector of real and imaginary parts,
C   WIo     respectively, of the closed-loop poles of interest;
C           the leading nap elements contain the real/imaginary parts
C           of the assigned poles, and the trailing np-nap elements
C           contain the unassigned poles.
C   Z     - (optional) the n-by-n orthogonal matrix which reduces the
C           closed-loop system state matrix A+B*F to real Schur form.
C   Ao    - (optional) the n-by-n matrix Z'*(A+B*F)*Z in a real Schur
C           form. The leading nfp-by-nfp diagonal block of Ao
C           corresponds to the fixed poles. The trailing nup-by-nup
C           diagonal block of A corresponds to the uncontrollable poles.
C
C Comments
C   1. Not all uncontrollable poles of the pair (A,B) are necessarily
C      detected by the eigenvalue assignment algorithm. Undetected
C      uncontrollable poles may exist if nfp > 0 and/or np < n-nfp.
C   2. To assign all controllable poles, set alpha = -Inf, for
C      discr = 0, and alpha = 0, for discr = 1.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, June 2002.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, July 2002,
C   Apr. 2009, Dec. 2012.
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
      CHARACTER         DICO
      INTEGER           INFO, IWARN, LDA, LDB, LDF, LDWORK, LDZ, M, N,
     $                  NAP, NFP, NP, NUP
      DOUBLE PRECISION  ALPHA, TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), DWORK(:), F(:,:),
     $                                 WI(:), WR(:), Z(:,:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           IP, ITMP, DISCR, NCW, NRW
      DOUBLE PRECISION  TEMP, TMP(3)
C
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C
C .. External subroutines ..
      EXTERNAL          SB01BD
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN, SQRT
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'POLASS requires at least 4 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'POLASS requires 1 output argument' )
      END IF
C
C   A(nxn), B(nxm), WR(np), WI(np) (, tol, discr, alpha).
C
      N = mxGetM( PRHS(1) )
      M = mxGetN( PRHS(2) )
C
      IF ( mxGetN( PRHS(1) ).NE.N ) THEN
         WRITE( TEXT, '(''A must have '',I5,'' rows and columns'')' ) N
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(2) ).NE.N ) THEN
         WRITE( TEXT, '(''B must have '',I5,'' rows'')' ) N
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
C
      NRW = mxGetM( PRHS(3) )
      NCW = mxGetN( PRHS(3) )
      NP  = NRW*NCW
      IF ( MIN( NRW, NCW ).NE.1 .AND. NP.GT.N ) THEN
         WRITE( TEXT, '(''WR must be a vector with '',I5,'' elements'',
     $                  '' at most'')' ) N
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'WR must be a real vector' )
      END IF
C
      NRW = mxGetM( PRHS(4) )
      NCW = mxGetN( PRHS(4) )
      IF ( MIN( NRW, NCW ).NE.1 .AND. NRW*NCW.NE.NP ) THEN
         WRITE( TEXT, '(''WI must be a vector with '',I5,'' elements''
     $                 )' ) NP
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'WI must be a real vector' )
      END IF
C
C   tol
C
      TOL = -ONE
      IP  = 5
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TOL, 1 )
         IP = IP + 1
      END IF
C
C   discr
C
      DISCR = 0
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'DISCR must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'DISCR must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ITMP = TEMP
         IF ( ITMP.LT.0 .OR. ITMP.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'DISCR has 0 or 1 the only admissible values' )
         END IF
         DISCR = ITMP
         IP  = IP + 1
      END IF
C
      ALPHA = -SQRT( DLAMCH( 'Epsilon' ) )
      IF ( DISCR.EQ.0 ) THEN
         DICO  = 'C'
      ELSE
         DICO  = 'D'
         ALPHA = ONE + ALPHA
      END IF
C
C   alpha
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IF ( DISCR.EQ.1 .AND. TEMP.LT.ZERO ) THEN
            CALL mexErrMsgTxt
     $           ( 'ALPHA must be positive in the discrete-time case' )
         END IF
         ALPHA = TEMP
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA = MAX( 1, N )
      LDF = MAX( 1, M )
      LDB = LDA
      LDZ = LDA
C
C   ldwork
C
      LDWORK = MAX( 1, 5*M, 5*N, 2*N + 4*M )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( A( LDA, N ), B( LDB, M ), DWORK( LDWORK ),
     $           F( LDF, N ), WI( NP ), WR( NP ), Z( LDZ, N ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), WR, NP )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), WI, NP )
C
C Do the actual computations.
C
      CALL SB01BD( DICO, N, M, NP, ALPHA, A, LDA, B, LDB, WR, WI, NFP,
     $             NAP, NUP, F, LDF, Z, LDZ, TOL, DWORK, LDWORK, IWARN,
     $             INFO )
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( M, N, 0 )
      CALL mxCopyReal8ToPtr( F, mxGetPr( PLHS(1) ), M*N )
      IF ( NLHS.GE.2 ) THEN
         TMP(1)  = NFP
         TMP(2)  = NAP
         TMP(3)  = NUP
         PLHS(2) = mxCreateDoubleMatrix( 3, 1, 0 )
         CALL mxCopyReal8ToPtr( TMP, mxGetPr( PLHS(2) ), 3 )
         IF ( NLHS.GE.3 ) THEN
            PLHS(3) = mxCreateDoubleMatrix( NP, 1, 0 )
            CALL mxCopyReal8ToPtr( WR, mxGetPr( PLHS(3) ), NP )
            IF ( NLHS.GE.4 ) THEN
               PLHS(4) = mxCreateDoubleMatrix( NP, 1, 0 )
               CALL mxCopyReal8ToPtr( WI, mxGetPr( PLHS(4) ), NP )
               IF ( NLHS.GE.5 ) THEN
                  PLHS(5) = mxCreateDoubleMatrix( N, N, 0 )
                  CALL mxCopyReal8ToPtr( Z, mxGetPr( PLHS(5) ), N*N )
                  IF ( NLHS.GE.6 ) THEN
                     PLHS(6) = mxCreateDoubleMatrix( N, N, 0 )
                     CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(6) ), N*N )
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( A, B, DWORK, F, WI, WR, Z )
C
C Error and warning handling.
C
      IF ( IWARN.NE.0 ) THEN
         WRITE( TEXT, '(''  IWARN = '',I4,'' ON EXIT FROM SB01BD'')'
     $        ) IWARN
         CALL mexPrintf( TEXT )
      END IF
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB01BD'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of POLASS ***
      END
