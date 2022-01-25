C ARESOLC.F - Gateway function for solving algebraic Riccati equations
C             using SLICOT routines SB02RD, SB02ND, SB02MT and SB02OD.
C             SB02RD is a more accurate version of SB02MD, also
C             computing condition numbers for Riccati equations.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [X(,F),ev(,db),rcond1(,acc)] = aresolc(method,A,Q,R,B,L,flag)
C       [X,ev(,db),rcond1(,acc)] = aresolc(method,A,Q,G,flag)
C       [Z,scale,ev(,db),rcond1] = aresolc(method,A,Q,R,B,L,flag)
C       [Z,scale,ev(,db),rcond1] = aresolc(method,A,Q,G,flag)
C
C   method = 1:        [X,F,ev,rcond1] = aresolc(1,A,Q,R,B,L,flag)
C                    [X,ev,rcond1,acc] = aresolc(1,A,Q,R,B,L,flag)
C                    [X,ev,rcond1,acc] = aresolc(1,A,Q,G,flag)
C                  [Z,scale,ev,rcond1] = aresolc(1,A,Q,R,B,L,flag)
C                  [Z,scale,ev,rcond1] = aresolc(1,A,Q,G,flag)
C   method = 2:     [X,F,ev,db,rcond1] = aresolc(2,A,Q,R,B,L,flag)
C                     [X,ev,db,rcond1] = aresolc(2,A,Q,R,B,L,flag)
C                     [X,ev,db,rcond1] = aresolc(2,A,Q,G,flag)
C               [Z,scale,ev,db,rcond1] = aresolc(2,A,Q,R,B,L,flag)
C               [Z,scale,ev,db,rcond1] = aresolc(2,A,Q,G,flag)
C   method = 3:     [X,F,ev,db,rcond1] = aresolc(3,A,Q,R,B,L,flag)
C                     [X,ev,db,rcond1] = aresolc(3,A,Q,R,B,L,flag)
C               [Z,scale,ev,db,rcond1] = aresolc(3,A,Q,R,B,L,flag)
C
C Purpose:
C   To solve the continuous-time algebraic Riccati equations
C                                      -1
C      0 = Q + A'*X + X*A - (L + X*B)*R  *(L + X*B)' ,             (1a)
C
C   or the discrete-time algebraic Riccati equations
C                                            -1
C      X = A'*X*A - (L + A'*X*B)*(R + B'*X*B)  (A'*X*B + L)' + Q;  (1b)
C
C   or, in alternate forms,
C
C      0 = Q + A'*X + X*A - X*G*X,                                 (2a)
C
C   or
C                            -1
C      X = Q + A'*X*(I + G*X)  *A.                                 (2b)
C
C   The function can also be used to compute an orthogonal basis of the
C   subspace which generates the solution.
C
C   When method = 1, the Schur method is used on Hamiltonian matrices
C                _      _
C             [  A     -G  ]         [ A  -G ]
C             [  _      _  ]   or    [       ]                     (3a)
C             [ -Q     -A' ]         [-Q  -A']
C   or symplectic matrices
C
C        [  _ -1       _ -1  _    ]       [   -1         -1     ]
C        [  A          A   * G    ]       [  A          A  *G   ]
C        [ _  _-1   _    _ _ -1 _ ]   or  [    -1          -1   ]  (3b)
C        [ Q* A     A' + Q*A   *G ]       [ Q*A    A' + Q*A  *G ]
C
C         _          -1     _          -1     _      -1
C   where A = A - B*R  *L', Q = Q - L*R  *L', G = B*R  *B'.
C
C   When method = 2, the generalized Schur method is used for discrete-
C   time problems, on the pencil
C               _                  _
C             [ A   0 ]      [ I   G  ]
C           a [ _     ]  - b [     _  ],                           (4a)
C             [-Q   I ]      [ 0   A' ]
C         _  _     _
C   where A, Q and G are defined above, or on the pencil
C
C             [ A   0 ]      [ I   G  ]
C           a [       ]  - b [        ],                           (4b)
C             [-Q   I ]      [ 0   A' ]
C
C   for (1b) or (2b).
C
C   When method = 3, the generalized Schur method is used on the
C   extended pencils
C
C             [ A  0  B ]      [ I   0  0 ]
C           a [ Q  A' L ] -  b [ 0  -I  0 ],                       (5a)
C             [ L' B' R ]      [ 0   0  0 ]
C   or
C
C             [ A  0  B ]      [ I  0  0 ]
C           a [ Q -I  L ] -  b [ 0 -A' 0 ],                        (5b)
C             [ L' 0  R ]      [ 0 -B' 0 ]
C
C   to solve the equations (1a) or (1b).
C
C Input parameters:
C   method - integer option to indicate the method to be used:
C            = 1 : use ordinary Schur method on matrix (3a) or (3b).
C            = 2 : use generalized Schur method on pencil (4a) or (4b).
C            = 3 : use generalized Schur method on pencil (5a) or (5b).
C   A      - real N-by-N system state matrix.
C   Q      - real symmetric N-by-N state weighting matrix.
C   R      - real symmetric M-by-M input weighting matrix.
C            When method = 1, or 2, R must be nonsingular.
C   B      - real N-by-M input matrix.
C   L      - real N-by-M coupling matrix.
C   G      - real symmetric N-by-N matrix.
C   flag   - (optional) vector containing options and parameters.
C            method = 1 : flag has length 4
C                flag(1) = 0 : solve the continuous-time equation
C                              (1a) or (2a); otherwise, solve the
C                              discrete-time equation (1b) or (2b).
C                flag(2) = 0 : compute the stabilizing solution;
C                              otherwise, compute the anti-stabilizing
C                              solution.
C                flag(3) = 0 : use a scaling strategy;
C                              otherwise, do not use a scaling strategy.
C                flag(4) = 0 : compute both the solution X and the
C                              feedback gain matrix F;
C                        = 1 : compute the solution X only;
C                        = 2 : compute the solution X and the accuracy
C                              estimates (separation, reciprocal
C                              condition number and forward error bound);
C                        < 0 : compute the invariant subspace.
C            Default:    flag(1:4) = [0,0,0,1].
C
C            method = 2 : flag has length 3
C                flag(1) = 0 : compute the stabilizing solution;
C                              otherwise, compute the anti-stabilizing
C                              solution.
C                flag(2) = 0 : compute both the solution X and the
C                              feedback gain matrix F;
C                        > 0 : compute the solution X only;
C                        < 0 : compute the deflating subspace.
C                flag(3)     : tolerance to check the singularity of
C                              the matrix pencil. If tol <= 0, tol will
C                              be replaced by EPS, the machine epsilon.
C            Default:    flag(1:3) = [0,1,0].
C
C            method = 3 : flag has length 4
C                flag(1) = 0 : solve the continuous-time equation (1a);
C                              otherwise, solve the discrete-time
C                              equation (1b).
C                flag(2) = 0 : compute the stabilizing solution;
C                              otherwise, compute the anti-stabilizing
C                              solution.
C                flag(3) = 0 : compute both the solution X and the
C                              feedback gain matrix F;
C                        > 0 : compute the solution only;
C                        < 0 : compute the related subspace.
C                flag(4)     : tolerance to check the singularity of
C                              the matrix pencil. If tol <= 0, tol will
C                              be replaced by EPS, the machine epsilon.
C            Default:    flag(1:4) = [0,0,1,0].
C
C Output parameters:
C   X      - N-by-N real symmetric solution of Riccati equations (1)
C            or (2).
C   F      - M-by-N real state feedback matrix, returned only when R
C            and B are input arguments, NLHS >= 2 and X is computed.
C   Z      - 2N-by-N real matrix with orthonormal columns, the basis of
C            the subspace related to the solution X.
C   scale  - the scaling factor used (internally), which should multiply
C            the submatrix Z(N+1:2N,:) to recover X from Z.
C   ev     - complex vector of length N.
C            method = 1   : the associated eigenvalues.
C            method = 2,3 : the "a" part of the associated eigenvalues
C                           (see (4) or (5)).
C   db     - real vector of length N containing the "b" part of
C            the associated eigenvalues (see (4) or (5)), i.e., the
C            associated eigenvalues with respect to X are ev(k)/db(k),
C            k=1,...,N.
C   rcond1 - estimate of the reciprocal of the 1-norm condition number
C            of the N-th order system of algebraic equations from which
C            the solution matrix X is obtained.
C   acc    - real vector of length 3 containing estimates of the
C            separation, reciprocal condition number and forward error
C            bound, repectively. These results can be obtained only if
C            meth = 1 and flag(4) = 2.
C
C Further Comments:
C   For METH = 2 or 3, this function behaves as aresol.
C
C Contributor:
C   V. Sima, Katholieke Univ. Leuven, Belgium, August 1999.
C
C Revisions:
C   S. Steer and V. Sima, February 2001.
C   V. Sima, Katholieke Univ. Leuven, Belgium, September 2001.
C   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2002.
C   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2005,
C   Apr. 2009, Dec. 2012.
C
C **********************************************************************
C
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
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
      CHARACTER         DICO, HINV, JOB, SCAL, SORT, UPLO
      INTEGER           INFO, LDA, LDA1, LDB, LDB1, LDF, LDG, LDL, LDL1,
     $                  LDQ, LDR, LDR1, LDS, LDT, LDU, LDWORK, LDX,
     $                  M, N
      DOUBLE PRECISION  FERR, RCOND, RCOND1, RNORM, SEP, TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: BWORK(:)
      INTEGER,          ALLOCATABLE :: IPIV(:), IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), A1(:,:), B(:,:), B1(:,:),
     $                                 BETA(:), DWORK(:), F(:,:),
     $                                 G(:,:), L(:,:), L1(:,:), Q(:,:),
     $                                 R(:,:), R1(:,:), S(:,:), T(:,:),
     $                                 U(:,:), WI(:), WR(:), X(:,:)
      COMPLEX*16,       ALLOCATABLE :: EV(:)
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           IWARN, JOBX
      CHARACTER*120     TEXT
      INTEGER           FLAG(4), FLAGF, IP, ISIZE, LIWORK, LOUT(2),
     $                  METH, NN, NO
      DOUBLE PRECISION  ACC(3), DUM(1), FLAGR(4), SCALE, TEMP
C
C .. External functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C
C .. External subroutines ..
      EXTERNAL          DLACPY, SB02MT, SB02ND, SB02OD, SB02RD
C
C ..Intrinsic functions..
      INTRINSIC         DCMPLX, MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ARESOLC requires at least 4 input arguments' )
      ELSE IF ( NLHS.GT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ARESOLC requires at most 5 output arguments' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   meth
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'METHOD must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'METHOD must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      METH = TEMP
      IF ( METH.LT.1 .OR. METH.GT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'METHOD has 1, 2, or 3 the only admissible values' )
      END IF
      IF ( METH.EQ.3 .AND. NRHS.LT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ARESOLC requires at least 6 input arguments' )
      END IF
C
C   A(NxN), Q(NxN), (R(MxM), B(NxM), L(NxM),) (G(NxN))
C
      N = mxGetM( PRHS(2) )
C
      IF ( mxGetN( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(3) ).NE.N .OR. mxGetN( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'Q must have the same size as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'Q must be a real matrix' )
      END IF
C
      IP = 5
      IF ( NRHS.GE.6 ) THEN
         M  = mxGetM( PRHS(4) )
         IF ( mxGetN( PRHS(4) ).NE.M  ) THEN
            CALL mexErrMsgTxt( 'R must be a square matrix' )
         END IF
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
      ELSE IF ( METH.LE.2 ) THEN
         IF ( mxGetM( PRHS(4) ).NE.N .OR.
     $        mxGetN( PRHS(4) ).NE.N ) THEN
            CALL mexErrMsgTxt( 'G must have the same size as A' )
         END IF
         IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(4) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'G must be a real matrix' )
         END IF
      END IF
C
      FLAG(1) = 0
      FLAG(2) = 0
      FLAG(3) = 0
      FLAG(4) = 0
      IF ( METH.EQ.1 ) THEN
         FLAG(4) = 1
      ELSE IF ( METH.EQ.2 ) THEN
         FLAG(2) = 1
      ELSE
         FLAG(3) = 1
      END IF
      IF ( NRHS.GE.IP ) THEN
C
C   flag
C
         ISIZE = mxGetM( PRHS(IP) ) * mxGetN( PRHS(IP) )
         IF ( METH.NE.2 ) THEN
            IF ( ISIZE.GT.4 )
     $         CALL mexErrMsgTxt
     $              ( 'FLAG must be a vector with at most 4 elements' )
         ELSE
            IF ( ISIZE.GT.3 )
     $         CALL mexErrMsgTxt
     $              ( 'FLAG must be a vector with at most 3 elements' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'FLAG must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), FLAGR, ISIZE )
         IF ( ISIZE.GT.0 ) FLAG(1) = FLAGR(1)
         IF ( ISIZE.GT.1 ) FLAG(2) = FLAGR(2)
         IF ( ISIZE.GT.2 ) FLAG(3) = FLAGR(3)
         IF ( ISIZE.GT.3 ) FLAG(4) = FLAGR(4)
      END IF
C
      IF ( METH.EQ.1 ) THEN
         FLAGF = FLAG(4)
         IF ( NLHS.GT.4 ) THEN
            CALL mexErrMsgTxt
     $         ( 'ARESOLC requires at most 4 output arguments' )
         END IF
      ELSE IF ( METH.EQ.2 ) THEN
         FLAGF = FLAG(2)
      ELSE
         FLAGF = FLAG(3)
      END IF
C
      IF ( METH.GE.2 ) THEN
         IF ( FLAG(METH).GT.0 .AND. NLHS.GT.4 ) THEN
            CALL mexErrMsgTxt
     $         ( 'ARESOLC requires at most 4 output arguments' )
         END IF
      END IF
C
      IF ( METH.NE.2 ) THEN
         IF ( FLAG(1).EQ.0 ) THEN
            DICO = 'C'
         ELSE
            DICO = 'D'
         END IF
      ELSE
         DICO = 'D'
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms
C in DGEES or DGGES, and also in DGEQRF/DGEQLF (M*NB), DGETRI (N*NB),
C DSYTRF (M*NB), DORMQL (2*N*NB).
C
      NN  = 2*N
      LDA = MAX( 1, N )
      LDQ = LDA
      LDX = LDA
      IF ( METH.LE.2 ) THEN
         LDG = LDA
         LDS = MAX( 1, NN )
         IF ( METH.EQ.2 )
     $      LDU = LDS
         LIWORK = LDS
         IF ( FLAGF.EQ.0 .OR. NRHS.GE.6 )
     $      LIWORK = MAX( LIWORK, M )
         IF ( METH.EQ.1 .AND. FLAGF.EQ.2 )
     $      LIWORK = MAX( LIWORK, N*N )
      ELSE
         LDS = MAX( 1, NN+M )
         LDU = MAX( 1, NN )
         LIWORK = MAX( 1, M, NN )
      END IF
      IF ( NRHS.GE.6 ) THEN
         LDR = MAX( 1, M )
         LDB = LDA
         LDL = LDA
C
         IF ( METH.LE.2 .AND. FLAGF.EQ.0 .AND. NLHS.GE.2 ) THEN
            LDA1 = LDA
            LDR1 = LDR
            LDB1 = LDA
            LDL1 = LDA
         END IF
      END IF
C
      IF ( METH.EQ.1 ) THEN
         LDT = LDA
         LDU = LDA
         LDWORK = 5 + MAX( 1, 4*N*N + 8*N )
         IF ( NRHS.GE.6 )
     $      LDWORK = MAX( LDWORK, 3*M, N*M )
      ELSE
         LDT = LDS
         LDWORK = MAX( 7*(2*N + 1) + 16, 16*N )
         IF ( METH.EQ.2 ) THEN
            IF ( NRHS.GE.6 )
     $         LDWORK = MAX( LDWORK, 3*M, N*M )
         ELSE
            LDWORK = MAX( LDWORK, 2*N + M, 3*M )
         END IF
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE( A(LDA,N), Q(LDQ,N), X(LDX,N), EV(N), WR(NN), WI(NN),
     $          IWORK(LIWORK), DWORK(LDWORK), BWORK(NN) )
      IF ( METH.LE.2 )
     $   ALLOCATE( G(LDG,N) )
C
      IF ( NRHS.GE.6 ) THEN
         ALLOCATE( R(LDR,M), B(LDB,M), L(LDL,M), IPIV(M) )
         IF ( FLAGF.EQ.0 .AND. NLHS.GE.2 ) THEN
            IF ( METH.LE.2 ) THEN
               ALLOCATE( R1(LDR1,M), B1(LDB1,M), L1(LDL1,M) )
               IF ( LSAME( DICO, 'D' ) ) THEN
                  ALLOCATE( A1(LDA1,N) )
               ELSE
                  ALLOCATE( A1(1,1) )
               END IF
            END IF
            LDF = LDR
            ALLOCATE( F(LDF,N) )
         END IF
      END IF
C
      IF ( METH.EQ.1 ) THEN
         ALLOCATE( S(LDS,NN), T(LDT,N), U(LDU,N) )
      ELSE
         ALLOCATE( T(LDT,NN), U(LDU,NN), BETA(NN) )
         IF ( METH.EQ.2 ) THEN
            ALLOCATE( S(LDS,NN) )
         ELSE
            ALLOCATE( S(LDS,NN+M) )
         END IF
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), Q, N*N )
C
      IF ( NRHS.GE.6 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), R, M*M )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), B, N*M )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), L, N*M )
C
      ELSE
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), G, N*N )
      END IF
C
C Do the actual computations.
C
      UPLO = 'U'
      IF ( METH.EQ.1 ) THEN
         IF ( FLAGF.EQ.2 .AND. NLHS.EQ.4 ) THEN
            JOB  = 'A'
            JOBX = .FALSE.
         ELSE
            JOB  = 'X'
            JOBX = .TRUE.
         END IF
         SORT = 'U'
         HINV = 'D'
         IF ( FLAG(2).EQ.0 ) THEN
            IF ( LSAME( DICO, 'C' ) ) THEN
               SORT = 'S'
            END IF
         ELSE IF ( LSAME( DICO, 'D' ) ) THEN
            HINV = 'I'
         END IF
C
         IF ( FLAG(3).NE.0 ) THEN
            SCAL = 'N'
         ELSE
            SCAL = 'G'
         END IF
C
      ELSE IF ( METH.EQ.2 ) THEN
         TOL = FLAGR(3)
         IF ( FLAG(1).EQ.0 ) THEN
            SORT = 'S'
         ELSE
            SORT = 'U'
         END IF
C
      ELSE
         TOL = FLAGR(4)
         IF ( FLAG(2).EQ.0 ) THEN
            SORT = 'S'
         ELSE
            SORT = 'U'
         END IF
      END IF
C
      IF ( METH.LE.2 .AND. NRHS.GE.6 ) THEN
         IF ( FLAGF.EQ.0 .AND. NLHS.GE.2 ) THEN
            IF ( LSAME( DICO, 'D' ) )
     $          CALL DLACPY( 'Full', N, N, A, LDA, A1, LDA1 )
            CALL DLACPY( 'Full', M, M, R, LDR, R1, LDR1 )
            CALL DLACPY( 'Full', N, M, B, LDB, B1, LDB1 )
            CALL DLACPY( 'Full', N, M, L, LDL, L1, LDL1 )
         END IF
C
         CALL SB02MT( 'G', 'N', 'N', UPLO, N, M, A, LDA, B, LDB,
     $                Q, LDQ, R, LDR, L, LDL, IPIV, LOUT, G, LDG,
     $                IWORK, DWORK, LDWORK, INFO )
C
         IF ( INFO.NE.0 ) THEN
            WRITE( TEXT, '(''INFO = '', I4, '' ON EXIT FROM SB02MT'')'
     $           ) INFO
            GO TO 20
         END IF
      END IF
C
      IWARN = .FALSE.
      IF ( METH.EQ.1 ) THEN
         CALL SB02RD( JOB, DICO, HINV, 'N', UPLO, SCAL, SORT, 'N',
     $                'R', N, A, LDA, T, LDT, U, LDU, G, LDG, Q,
     $                LDQ, X, LDX, SEP, RCOND, FERR, WR, WI, S, LDS,
     $                IWORK, DWORK, LDWORK, BWORK, INFO )
C
         IWARN = INFO.EQ.7 .OR.
     $           ( FLAGF.LT.0 .AND. JOBX .AND. INFO.EQ.5 )
         IF ( IWARN ) THEN
            WRITE( TEXT,
     $          '(''Warning: INFO = '', I4, '' ON EXIT FROM SB02RD'')'
     $           ) INFO
         ELSE IF ( INFO.NE.0 ) THEN
            WRITE( TEXT, '(''INFO = '', I4, '' ON EXIT FROM SB02RD'')' )
     $             INFO
            GO TO 20
         END IF
         RCOND1 = DWORK(2)
         IF ( N.EQ.0 )
     $      SEP = ONE
         IF ( .NOT.JOBX ) THEN
            ACC(1) = SEP
            ACC(2) = RCOND
            ACC(3) = FERR
         END IF
C
      ELSE IF ( METH.EQ.2 ) THEN
C
         CALL SB02OD( 'D', 'G', 'N', UPLO, 'Z', SORT, N, M, M, A,
     $                LDA, G, LDG, Q, LDQ, DUM, 1, DUM, 1, RCOND1, X,
     $                LDX, WR, WI, BETA, S, LDS, T, LDT, U, LDU,
     $                TOL, IWORK, DWORK, LDWORK, BWORK, INFO )
C
      ELSE
C
         CALL SB02OD( DICO, 'B', 'N', UPLO, 'N', SORT, N, M, M, A,
     $                LDA, B, LDB, Q, LDQ, R, LDR, L, LDL, RCOND1, X,
     $                LDX, WR, WI, BETA, S, LDS, T, LDT, U, LDU,
     $                TOL, IWORK, DWORK, LDWORK, BWORK, INFO )
C
      END IF
C
      IF ( METH.GE.2 ) THEN
         IWARN = FLAGF.LT.0 .AND. INFO.EQ.6
         IF ( IWARN ) THEN
            WRITE( TEXT,
     $          '(''Warning: INFO = '', I4, '' ON EXIT FROM SB02OD'')'
     $           ) INFO
         ELSE IF ( INFO.NE.0 ) THEN
            WRITE( TEXT, '(''INFO = '', I4, '' ON EXIT FROM SB02OD'')' )
     $             INFO
            GO TO 20
         END IF
         SCALE = DWORK(3)
      END IF
      IF ( FLAGF.EQ.0 .AND. NRHS.GE.6 .AND. NLHS.GE.2 ) THEN
         IF ( METH.LE.2 ) THEN
            CALL SB02ND( DICO, 'N', UPLO, 'N', N, M, M, A1, LDA1, B1,
     $                   LDB1, R1, LDR1, IPIV, L1, LDL1, X, LDX, RNORM,
     $                   F, LDF, LOUT, IWORK, DWORK, LDWORK, INFO )
         ELSE
            CALL SB02ND( DICO, 'N', UPLO, 'N', N, M, M, A, LDA, B,
     $                   LDB, R, LDR, IPIV, L, LDL, X, LDX, RNORM,
     $                   F, LDF, LOUT, IWORK, DWORK, LDWORK, INFO )
         END IF
C
         IF ( INFO.NE.0 ) THEN
            WRITE( TEXT, '(''INFO = '', I4, '' ON EXIT FROM SB02ND'')' )
     $             INFO
            GO TO 20
         END IF
      END IF
C
C Copy output to MATLAB workspace.
C
      NO = 2
      IF ( INFO.EQ.0 .OR. IWARN ) THEN
         IF ( FLAGF.GE.0 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( X, mxGetPr( PLHS(1) ), N*N )
            IF ( FLAGF.EQ.0 .AND. NRHS.GE.6 .AND. NLHS.GE.NO ) THEN
               PLHS(NO) = mxCreateDoubleMatrix( M, N, 0 )
               CALL mxCopyReal8ToPtr( F, mxGetPr( PLHS(NO) ), M*N )
               NO = NO + 1
            END IF
         ELSE
            IF ( METH.EQ.1 ) THEN
               PLHS(1) = mxCreateDoubleMatrix( NN, N, 0 )
               CALL mxCopyReal8ToPtr( DWORK(6),
     $                      mxGetPr( PLHS(1) ), NN*N )
               IF ( JOBX .AND. NLHS.GE.NO ) THEN
                  PLHS(NO) = mxCreateDoubleMatrix( 1, 1, 0 )
                  CALL mxCopyReal8ToPtr( SEP, mxGetPr( PLHS(NO) ), 1 )
                  NO = NO + 1
               END IF
            ELSE
               PLHS(1) = mxCreateDoubleMatrix( NN, N, 0 )
               CALL mxCopyReal8ToPtr( U, mxGetPr( PLHS(1) ), NN*N )
               IF ( NLHS.GE.NO ) THEN
                  PLHS(NO) = mxCreateDoubleMatrix( 1, 1, 0 )
                  CALL mxCopyReal8ToPtr( SCALE, mxGetPr( PLHS(NO) ), 1 )
                  NO = NO + 1
               END IF
            END IF
         END IF
C
         IF ( NLHS.GE.NO ) THEN
            DO 10 I = 1, N
               EV(I) = DCMPLX( WR(I), WI(I) )
   10       CONTINUE
            PLHS(NO) = mxCreateDoubleMatrix( N, MIN( 1, N ), 1 )
            CALL mxCopyComplex16ToPtr( EV, mxGetPr( PLHS(NO) ),
     $                                 mxGetPi( PLHS(NO) ), N )
            IF ( METH.GE.2 ) THEN
               NO = NO + 1
               PLHS(NO) = mxCreateDoubleMatrix( N, MIN( 1, N ), 0 )
               CALL mxCopyReal8ToPtr( BETA, mxGetPr( PLHS(NO) ), N )
            END IF
         END IF
C
         NO = NO + 1
         IF ( NLHS.GE.NO ) THEN
            PLHS(NO) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( RCOND1, mxGetPr( PLHS(NO) ), 1 )
            NO = NO + 1
         END IF
C
         IF ( FLAGF.EQ.2 .AND. NLHS.GE.NO ) THEN
            PLHS(NO) = mxCreateDoubleMatrix( 3, 1, 0 )
            CALL mxCopyReal8ToPtr( ACC, mxGetPr( PLHS(NO) ), 3 )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
   20 CONTINUE
C
      DEALLOCATE( A, Q, S, T, U, X, EV, WR, WI, IWORK, DWORK, BWORK )
      IF ( METH.LE.2 )
     $   DEALLOCATE( G )
      IF ( NRHS.GE.6 ) THEN
         DEALLOCATE( R, B, L, IPIV )
         IF ( FLAGF.EQ.0 .AND. NLHS.GE.2 ) THEN
            IF ( METH.LE.2 )
     $         DEALLOCATE( A1, R1, B1, L1 )
            DEALLOCATE( F )
         END IF
      END IF
      IF ( METH.GE.2 )
     $   DEALLOCATE ( BETA )
C
C Error and warning handling.
C
      IF ( IWARN ) THEN
         CALL mexPrintf( TEXT )
      ELSE IF ( INFO.NE.0 ) THEN
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C
C *** Last line of ARESOLC ***
      END
