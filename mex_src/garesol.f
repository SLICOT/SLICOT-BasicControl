C GARESOL.F - Gateway function for solving descriptor algebraic Riccati
C             equations using SLICOT routine SG02AD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C       [X,F,ev,db,rcond1] = garesol(A,E,Q,R,B,L,flag)
C         [X,ev,db,rcond1] = garesol(A,E,Q,R,B,L,flag)
C         [X,ev,db,rcond1] = garesol(A,E,Q,G,flag)
C   [Z,scale,ev,db,rcond1] = garesol(A,E,Q,R,B,L,flag)
C   [Z,scale,ev,db,rcond1] = garesol(A,E,Q,G,flag)
C
C Purpose:
C   To solve the continuous-time algebraic Riccati equations
C                                              -1
C      0 = Q + A'*X*E + E'*X*A - (L + E'*X*B)*R  *(L + E'*X*B)' ,  (1a)
C
C   or the discrete-time algebraic Riccati equations
C                                                 -1
C      E'*X*E = A'*X*A - (L + A'*X*B)*(R + B'*X*B)  (L + A'*X*B)' + Q;
C                                                                  (1b)
C   or, in alternate forms,
C
C                 -1
C      0 = Q - L*R  *L' + A1'*X*E + E'*X*A1 - E'*X*G*X*E,          (2a)
C
C   or
C                      -1                     -1
C      E'*X*E = Q - L*R  *L' + A1'*X*(I + G*X)  *A1,               (2b)
C
C                -1                 -1
C   where G = B*R  *B', A1 = A - B*R  *L'. The optimal feedback gain
C   matrix is
C
C           -1
C      F = R  (L+E'XB)' ,        for (1a),
C
C     and
C                  -1
C      F = (R+B'XB)  (L+A'XB)' , for (2a).
C
C   The function can also be used to compute an orthogonal basis of the
C   subspace which generates the solution.
C
C   The generalized Schur method is used on the extended pencils
C
C             [ A   0   B ]      [ E   0   0 ]
C           a [ Q   A'  L ] -  b [ 0  -E'  0 ] ,                   (3a)
C             [ L'  B'  R ]      [ 0   0   0 ]
C   or
C
C             [ A   0   B ]      [ E   0   0 ]
C           a [ Q  -E'  L ] -  b [ 0  -A'  0 ] ,                   (3b)
C             [ L'  0   R ]      [ 0  -B'  0 ]
C
C   to solve the equations (1a) or (1b).
C
C Input parameters:
C   A      - real N-by-N system state matrix.
C   E      - real N-by-N descriptor system matrix.
C   Q      - normally, real symmetric N-by-N state weighting matrix.
C            If flag(5) <> 0, array Q stores the P-by-N factor C of Q,
C            Q = C'*C.
C   R      - normally, real symmetric M-by-M input weighting matrix.
C            If R is an input parameter and flag(6) <> 0, array R stores
C            the P-by-M factor D of R, R = D'*D.
C   B      - real N-by-M input matrix.
C   L      - real N-by-M coupling matrix.
C   G      - real symmetric N-by-N matrix.
C   flag   - (optional) vector of length 10 containing options:
C               flag(1)  = 0 : solve the continuous-time equation (1a);
C                              otherwise, solve the discrete-time
C                              equation (1b).
C               flag(2)  = 0 : compute the stabilizing solution;
C                              otherwise, compute the anti-stabilizing
C                              solution.
C               flag(3)  = 0 : compute both the solution X and the
C                              feedback gain matrix F;
C                        > 0 : compute the solution only;
C                        < 0 : compute the related subspace.
C               flag(4)      : tolerance to check the singularity of
C                              the matrix pencil. If tol <= 0, tol will
C                              be replaced by EPS, the machine epsilon.
C               flag(5)  = 0 : matrix Q is not factored; otherwise Q is
C                              assumed factored as Q = C'*C, where C is
C                              P-by-N.
C               flag(6)  = 0 : matrix R is not factored; otherwise R is
C                              assumed factored as R = D'*D, where D is
C                              P-by-M.
C               flag(7)  = 0 : the upper triangles of matrices Q and G,
C                              or Q and R, are stored; otherwise, lower
C                              triangles are stored.
C               flag(8)  = 0 : matrix L is zero; otherwise, L is given.
C               flag(9)  = 0 : use a scaling strategy (for given R, B);
C                              otherwise, do not use a scaling strategy.
C               flag(10) = 0 : iterative refinement should not be used
C                              for solving the system of algebraic
C                              equations giving the solution matrix X;
C                              otherwise, use iterative refinement.
C            Default:    flag(1:4) = [0,0,1,0,0,0,0,1,0,0].
C
C Output parameters:
C   X      - N-by-N real symmetric solution of Riccati equations (1)
C            or (2).
C   F      - M-by-N real state feedback matrix, returned only when R
C            and B are input arguments, the number of output parameters
C            is at least 2, and X is computed.
C   Z      - 2N-by-N real matrix with orthonormal columns, the basis of
C            the subspace related to the solution X.
C   scale  - the scaling factor used internally, which should multiply
C            the submatrix Z(N+1:2N,:) to recover X from Z.
C   ev     - complex vector of length N.
C            the "a" part of the associated eigenvalues in (3).
C   db     - real vector of length N containing the "b" part of
C            the associated eigenvalues in (3), i.e., the associated
C            eigenvalues with respect to X are ev(k)/db(k), k=1,...,N.
C   rcond1 - estimate of the reciprocal of the 1-norm condition number
C            of the N-th order system of algebraic equations from which
C            the solution matrix X is obtained.
C
C Contributor:
C   V. Sima, Katholieke Univ. Leuven, Belgium, March 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Jan. 2005,
C   Apr. 2009, Dec. 2012, Dec. 2014, Jan. 2015, Apr. 2017.
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
      CHARACTER         ACC, DICO, FACT, FACTR, JOBB, JOBL, SCAL, SORT,
     $                  UPLO
      INTEGER           INFO, IWARN, LDA, LDB, LDE, LDF, LDL, LDQ, LDR,
     $                  LDS, LDT, LDU, LDWORK, LDX, M, N, P
      DOUBLE PRECISION  RCOND1, RNORM, TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: BWORK(:)
      INTEGER,          ALLOCATABLE :: IPIV(:), IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), BETA(:),
     $                                 DWORK(:), E(:,:), F(:,:), L(:,:),
     $                                 Q(:,:), R(:,:), S(:,:), T(:,:),
     $                                 U(:,:), WI(:), WR(:), X(:,:)
      COMPLEX*16,       ALLOCATABLE :: EV(:)
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           WITHF, WITHL, WITHR
      CHARACTER*120     TEXT
      INTEGER           FLAG(10), FLAGF, I, IP, ISIZE, LIWORK, LOUT(2),
     $                  NCB, NN, NO, NP
      DOUBLE PRECISION  FLAGR(10), SCALE
C
C     .. External Functions ..
      LOGICAL           LSAME, MA02HD
      EXTERNAL          LSAME, MA02HD
C
C .. External subroutines ..
      EXTERNAL          SB02ND, SG02AD, SG02ND
C
C ..Intrinsic functions..
      INTRINSIC         DCMPLX, MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'GARESOL requires at least 4 input arguments' )
      ELSE IF ( NLHS.GT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'GARESOL requires at most 5 output arguments' )
      END IF
C
C   A(NxN), E(NxN), Q(NxN), (R(MxM), B(NxM), L(NxM),) (G(NxN))
C
      N  = mxGetM( PRHS(1) )
      P  = mxGetM( PRHS(3) )
      NP = mxGetM( PRHS(4) )
      M  = mxGetN( PRHS(4) )
C
      IF ( mxGetN( PRHS(1) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(2) ).NE.N .OR. mxGetN( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'E must have the same dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'E must be a real matrix' )
      END IF
C
      IF ( mxGetN( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $           ( 'Q must have the same column dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'Q must be a real matrix' )
      END IF
C
      IF ( NRHS.GE.6 ) THEN
         JOBB  = 'B'
         WITHR = .TRUE.
         IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(4) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'R must be a real matrix' )
         END IF
C
         IF ( mxGetM( PRHS(5) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $           ( 'B must have the same row dimension as A' )
         END IF
         IF ( mxGetN( PRHS(5) ).NE.M ) THEN
            CALL mexErrMsgTxt
     $           ( 'B must have the same column dimension as R' )
         END IF
         IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(5) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'B must be a real matrix' )
         END IF
C
         IF ( mxGetM( PRHS(6) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $           ( 'L must have the same row dimension as A' )
         END IF
         IF ( mxGetN( PRHS(6) ).NE.M ) THEN
            CALL mexErrMsgTxt
     $           ( 'L must have the same column dimension as R' )
         END IF
         IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'L must be a real matrix' )
         END IF
C
         IP = 7
      ELSE
         JOBB  = 'G'
         WITHR = .FALSE.
         IF ( NP.NE.N .OR. M.NE.N ) THEN
            CALL mexErrMsgTxt
     $           ( 'G must have the same dimension as A' )
         END IF
         IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(4) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'G must be a real matrix' )
         END IF
         IP = 5
      END IF
C
      DO 10 I = 1, 10
         FLAG(I) = 0
   10 CONTINUE
C
      FLAG(3) = 1
      FLAG(8) = 1
C
      IF ( NRHS.GE.IP ) THEN
C
C   flag
C
         ISIZE = mxGetM( PRHS(IP) ) * mxGetN( PRHS(IP) )
         IF ( ISIZE.GT.10 )
     $      CALL mexErrMsgTxt
     $           ( 'FLAG must be a vector with at most 10 elements' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'FLAG must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), FLAGR, ISIZE )
C
         DO 20 I = 1, ISIZE
            FLAG(I) = FLAGR(I)
   20    CONTINUE
C
      END IF
C
      IF ( FLAG(1).EQ.0 ) THEN
         DICO = 'C'
      ELSE
         DICO = 'D'
      END IF
C
      IF ( FLAG(2).EQ.0 ) THEN
         SORT = 'S'
      ELSE
         SORT = 'U'
      END IF
      FLAGF = FLAG(3)
      WITHF = FLAGF.EQ.0 .AND. WITHR .AND. NLHS.GE.2
      TOL   = FLAGR(4)
C
      IF ( FLAG(5).EQ.0 .AND. P.NE.N )
     $   CALL mexErrMsgTxt( 'Q must have the same size as A' )
      IF ( WITHR ) THEN
         IF ( FLAG(6).EQ.0 ) THEN
            FACTR = 'N'
            IF ( NP.NE.M )
     $         CALL mexErrMsgTxt( 'R must be a square matrix' )
         ELSE
            FACTR = 'D'
            IF ( FLAG(5).NE.0 .AND. NP.NE.P ) THEN
               CALL mexErrMsgTxt
     $            ( 'D must must have the same row dimension as C' )
            END IF
         END IF
      ELSE
         FLAG(6) = 0
         FLAG(8) = 0
      END IF
C
      IF ( FLAG(5).EQ.0 ) THEN
         IF ( FLAG(6).EQ.0 ) THEN
            FACT = 'N'
         ELSE
            FACT = 'D'
         END IF
      ELSE
         IF ( FLAG(6).EQ.0 ) THEN
            FACT = 'C'
         ELSE
            FACT = 'B'
         END IF
      END IF
C
      IF ( FLAG(7).EQ.0 ) THEN
         UPLO = 'U'
      ELSE
         UPLO = 'L'
      END IF
C
      IF ( FLAG(8).EQ.0 ) THEN
         JOBL = 'Z'
         WITHL = .FALSE.
      ELSE
         JOBL = 'N'
         WITHL = .TRUE.
      END IF
C
      IF ( FLAG(9).EQ.0 .AND. WITHR ) THEN
         SCAL = 'G'
      ELSE
         SCAL = 'N'
      END IF
C
      IF ( FLAG(10).EQ.0 ) THEN
         ACC = 'N'
      ELSE
         ACC = 'R'
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms
C in DGEES or DGGES, and also in DGEQRF/DGEQLF (M*NB), DGETRI (N*NB),
C DSYTRF (M*NB), DORMQL (2*N*NB).
C
      NN  = 2*N
      LDA = MAX( 1, N )
      LDE = LDA
      LDB = LDA
      LDU = MAX( 1, NN )
      LDL = 1
C
      LDWORK = MAX( 7*(2*N + 1) + 16, 16*N )
      IF ( WITHR ) THEN
         NCB    = M
         LDS    = MAX( 1, NN + M )
         LIWORK = MAX( 1, M, NN )
         LDWORK = MAX( LDWORK, 2*N + M, 3*M )
         IF ( WITHL )
     $      LDL = LDA
      ELSE
         NCB    = N
         LDS    = LDU
         LIWORK = LDU
      END IF
      LDQ = MAX( 1, P )
      LDR = MAX( 1, NP )
      LDT = LDS
      LDX = LDA
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE( A(LDA,N), E(LDE,N), B(LDB,NCB), Q(LDQ,N), U(LDU,NN),
     $          EV(N), WR(NN), WI(NN), IWORK(LIWORK), DWORK(LDWORK),
     $          BWORK(NN), T(LDT,NN), BETA(NN), S(LDS,LDS), X(LDX,N) )
      IF ( WITHR ) THEN
         ALLOCATE( R(LDR,M) )
         IF ( WITHL )
     $      ALLOCATE( L(LDL,M) )
         IF ( WITHF ) THEN
            LDF = LDR
            ALLOCATE( F(LDF,N), IPIV(M) )
         END IF
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), E, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), Q, P*N )
C
      IF ( WITHR ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), R, NP*M )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), B, N*M )
         IF ( WITHL )
     $      CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), L, N*M )
C
      ELSE
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), B, N*N )
      END IF
C
C Do the actual computations.
C
      CALL SG02AD( DICO, JOBB, FACT, UPLO, JOBL, SCAL, SORT, ACC, N,
     $             M, P, A, LDA, E, LDE, B, LDB, Q, LDQ, R, LDR, L, LDL,
     $             RCOND1, X, LDX, WR, WI, BETA, S, LDS, T, LDT, U, LDU,
     $             TOL, IWORK, DWORK, LDWORK, BWORK, IWARN, INFO )
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '', I4, '' ON EXIT FROM SG02AD'')' )
     $          INFO
         GO TO 40
      END IF
C
      SCALE = DWORK(4)
      IF ( WITHF ) THEN
         IF ( FLAG(1).EQ.0 ) THEN
            IF ( N.EQ.0 .OR. MA02HD( 'Full', N, N, ONE, E, LDE ) ) THEN
               CALL SB02ND( DICO, FACTR, UPLO, JOBL, N, M, P, A, LDA, B,
     $                      LDB, R, LDR, IPIV, L, LDL, X, LDX, RNORM, F,
     $                      LDF, LOUT, IWORK, DWORK, LDWORK, INFO )
            ELSE
               CALL SG02ND( DICO, 'GenE', 'Kalman', 'No XE', FACTR,
     $                      UPLO, JOBL, 'NoTrans', N, M, P, A, LDA, E,
     $                      LDE, B, LDB, R, LDR, IPIV, L, LDL, X, LDX,
     $                      RNORM, F, LDF, A, 1, A, 1, LOUT, IWORK,
     $                      DWORK, LDWORK, INFO )
            END IF
         ELSE
            CALL SB02ND( DICO, FACTR, UPLO, JOBL, N, M, P, A, LDA, B,
     $                   LDB, R, LDR, IPIV, L, LDL, X, LDX, RNORM, F,
     $                   LDF, LOUT, IWORK, DWORK, LDWORK, INFO )
         END IF
C
         IF ( INFO.NE.0 ) THEN
            WRITE( TEXT, '('' INFO = '', I4, '' ON EXIT FROM SB02ND'')')
     $             INFO
            GO TO 40
         END IF
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         IF ( FLAGF.GE.0 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( X, mxGetPr( PLHS(1) ), N*N )
            NO = 2
            IF ( WITHF ) THEN
               PLHS(2) = mxCreateDoubleMatrix( M, N, 0 )
               CALL mxCopyReal8ToPtr( F, mxGetPr( PLHS(2) ), M*N )
               NO = NO + 1
            END IF
         ELSE
            PLHS(1) = mxCreateDoubleMatrix( NN, N, 0 )
            CALL mxCopyReal8ToPtr( U, mxGetPr( PLHS(1) ), NN*N )
            PLHS(2) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( SCALE, mxGetPr( PLHS(2) ), 1 )
            NO = 3
         END IF
C
         IF ( NLHS.GE.NO ) THEN
C
            DO 30 I = 1, N
               EV(I) = DCMPLX( WR(I), WI(I) )
   30       CONTINUE
C
            NP = MIN( 1, N )
            PLHS(NO) = mxCreateDoubleMatrix( N, NP, 1 )
            CALL mxCopyComplex16ToPtr( EV, mxGetPr( PLHS(NO) ),
     $                                 mxGetPi( PLHS(NO) ), N*NP )
            NO = NO + 1
            PLHS(NO) = mxCreateDoubleMatrix( N, NP, 0 )
            CALL mxCopyReal8ToPtr( BETA, mxGetPr( PLHS(NO) ), N*NP )
         END IF
C
         NO = NO + 1
         IF ( NLHS.GE.NO ) THEN
            PLHS(NO) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( RCOND1, mxGetPr( PLHS(NO) ), 1 )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
   40 CONTINUE
C
      DEALLOCATE( A, E, B, Q, S, U, EV, WR, WI, IWORK, DWORK, BWORK,
     $            T, X, BETA )
      IF ( WITHR ) THEN
         DEALLOCATE( R )
         IF ( WITHL )
     $      DEALLOCATE( L )
         IF ( WITHF )
     $      DEALLOCATE( F, IPIV )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         CALL mexErrMsgTxt( TEXT )
      ELSE IF ( IWARN.NE.0 ) THEN
         WRITE( TEXT, '(''  IWARN = '',I4,'' ON EXIT FROM SG02AD'')'
     $        ) IWARN
         CALL mexPrintf( TEXT )
      END IF
C
      RETURN
C
C *** Last line of GARESOL ***
      END
