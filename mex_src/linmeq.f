C LINMEQ.F - Gateway function for solving Sylvester and Lyapunov matrix
C            equations using SLICOT routines SB04MD, SB04ND, SB04PD,
C            SB04QD, SB04RD, SB03MD, and SB03OD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [X(,sep)] = linmeq(task,A(,B),C,flag,trans(,schur))
C
C   task = 1 :      [X] = linmeq(1,A,B,C,flag,trans,schur)
C   task = 2 :  [X,sep] = linmeq(2,A,C,flag,trans)
C                   [X] = linmeq(2,A,C,flag,trans)
C   task = 3 :      [X] = linmeq(3,A,C,flag,trans)
C
C Purpose:
C   To solve the Sylvester and Lyapunov linear matrix equations
C
C   task = 1:
C
C         op(A)*X + X*op(B) = C,                          (1a)
C
C         op(A)*X*op(B) + X = C,                          (1b)
C
C   task = 2:
C
C         op(A)'*X + X*op(A) = C,                         (2a)
C
C         op(A)'*X*op(A) - X = C,                         (2b)
C
C   task = 3:
C
C         op(A)'*(op(X)'*op(X)) + (op(X)'*op(X))*op(A) =
C                               -  op(C)'*op(C),          (3a)
C
C         op(A)'*(op(X)'*op(X))*op(A) - op(X)'*op(X) =
C                                     - op(C)'*op(C),     (3b)
C
C   where op(M) = M, if trans = 0, and op(M) = M', if trans = 1.
C
C Input parameters:
C   task  - integer option to determine the equation type:
C           = 1 : solve the Sylvester equation (1a) or (1b);
C           = 2 : solve the Lyapunov equation (2a) or (2b);
C           = 3 : solve for the Cholesky factor op(X) the Lyapunov
C                 equation (3a) or (3b).
C   A     - real coefficient N-by-N matrix.
C           When task = 3, matrix A must be stable.
C   B     - another real coefficient M-by-M matrix for
C           equations (1a) or (1b).
C   C     - right hand side matrix.
C           task = 1 : C is N-by-M;
C           task = 2 : C is N-by-N symmetric;
C           task = 3 : op(C) is P-by-N.
C   flag  - (optional) integer vector of length 3 or 2 containing
C           options.
C           task = 1 : flag has length 3
C                flag(1) = 0 : solve the continuous-time equation (1a);
C                              otherwise, solve the discrete-time
C                              equation (1b).
C                flag(2) = 1 : A is (quasi) upper triangular;
C                          2 : A is upper Hessenberg;
C                              otherwise, A is in general form.
C                flag(3) = 1 : B is (quasi) upper triangular;
C                          2 : B is upper Hessenberg;
C                              otherwise, B is in general form.
C           task = 2 : flag has length 2
C                flag(1) = 0 : solve continuous-time equation (2a);
C                              otherwise, solve discrete-time
C                              equation (2b).
C                flag(2) = 1 : A is (quasi) upper triangular;
C                              otherwise, A is in general form.
C           task = 3 : flag has length 2
C                flag(1) = 0 : solve continuous-time equation (3a);
C                              otherwise, solve discrete-time
C                              equation (3b).
C                flag(2) = 1 : A is (quasi) upper triangular;
C                              otherwise, A is in general form.
C           Default:    flag(1) = 0, flag(2) = 0 (, flag(3) = 0).
C   trans - (optional) integer specifying a transposition option.
C           trans = 0 : solve the equations (1) - (3) with op(M) = M.
C           trans = 1 : solve the equations (1) - (3) with op(M) = M'.
C           trans = 2 : solve the equations (1) with op(A) = A',
C                                                    op(B) = B.
C           trans = 3 : solve the equations (1) with op(A) = A,
C                                                    op(B) = B'.
C           Default:    trans = 0.
C   schur - (optional) integer specifying whether the Hessenberg-Schur
C           or Schur method should be used.
C           Available for task = 1.
C           schur = 1 : Hessenberg-Schur method (one matrix is reduced
C                       to Schur form).
C           schur = 2 : Schur method (two matrices are reduced to Schur
C                       form).
C           Default:    schur = 1.
C
C Output parameters:
C   X     - solution of the equation (or its Cholesky factor for (3)).
C   sep   - (optional) estimator of Sep(op(A),-op(A)') for (2.a) or
C           Sepd(A,A') for (2.b).
C
C Comments:
C   1. For equation (1a) or (1b), when schur = 1, the Hessenberg-Schur
C      method is used, reducing one matrix to Hessenberg form and the
C      other one to a real Schur form.
C      Otherwise, both matrices are reduced to real Schur forms.
C      If one or both matrices are already reduced to Schur/Hessenberg
C      forms, this could be specified by flag(2) and flag(3).
C      For general matrices, the Hessenberg-Schur method could be
C      significantly more efficient than the Schur method.
C   2. For equation (3a) or (3b), the computed matrix X is the Cholesky
C      factor of the solution, i.e., the real solution is op(X)'*op(X),
C      where X is an N-by-N upper triangular matrix.
C
C Contributor:
C   H. Xu, TU Chemnitz, FR Germany, Dec. 1998.
C
C Revisions:
C   V. Sima, Katholieke Univ. Leuven, Belgium, May 1999.
C   V. Sima, Katholieke Univ. Leuven, Belgium, May 2000.
C   D. Sima, University of Bucharest, Romania, May 2000.
C   V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 2000.
C   V. Sima, Research Institute for Informatics, Bucharest, March 2006.
C   V. Sima, Research Institute for Informatics, Bucharest, Feb. 2009,
C   Apr. 2009, July 2011, Dec. 2012.
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
      CHARACTER         DICO, FACT, FACTA, FACTB, JOB, SCHU, TRANA,
     $                  TRANB, ULA, ULB
      INTEGER           INFO, ISGN, LDA, LDB, LDC, LDU, LDV, LDWORK, M,
     $                  N, P
      DOUBLE PRECISION  FERR, SCALE, SEP, TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), DWORK(:),
     $                                 U(:,:), V(:,:), WR(:), WI(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           BWORK(1), PERTRB
      INTEGER           FLAG(3), IB, IP, ISIZE, J, LDWO, LDW1, LDW2,
     $                  LIWORK, MC, MXMN, NB1, NB2, NC, NM, NSCHUR,
     $                  TASK, TRANS
      DOUBLE PRECISION  DUM(1), FLAGR(3), TEMP
C
C .. External functions ..
      LOGICAL           LSAME, SELECT
      INTEGER           ILAENV
      EXTERNAL          LSAME, ILAENV, SELECT
C
C .. External subroutines ..
      EXTERNAL          DLACPY, DGEES, DLASET, DSWAP, DTRSYL, SB03MD,
     $                  SB03OD, SB04MD, SB04ND, SB04PD, SB04PY, SB04QD,
     $                  SB04RD
C
C ..Intrinsic functions..
      INTRINSIC         MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'LINMEQ requires at least 3 input arguments' )
      ELSE IF ( NLHS.GT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'LINMEQ requires at most 2 output arguments' )
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
      IF ( TASK.LT.1 .OR. TASK.GT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'TASK has 1, 2, or 3 the only admissible values' )
      END IF
C
      IF ( TASK.EQ.1 ) THEN
         IF ( NRHS.LT.4 ) THEN
            CALL mexErrMsgTxt
     $          ( 'LINMEQ requires at least 4 input arguments' )
         END IF
         IP = 6
      ELSE
         IP = 5
      END IF
C
C   trans
C
      TRANS = 0
      IF ( NRHS.GE.IP ) THEN
         ISIZE = mxGetM( PRHS(IP) ) * mxGetN( PRHS(IP) )
         IF ( ISIZE.GT.1 ) THEN
            CALL mexErrMsgTxt( 'TRANS must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TRANS must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, ISIZE )
         IF ( ISIZE.GT.0 ) TRANS = TEMP
         IF ( TASK.EQ.1 .AND. ( TRANS.LT.0 .OR. TRANS.GT.3 ) ) THEN
            CALL mexErrMsgTxt
     $         ( 'TRANS has 0, 1, 2, or 3 the only admissible values' )
         ELSE IF ( TASK.NE.1 .AND. ( TRANS.LT.0 .OR. TRANS.GT.1 ) ) THEN
            CALL mexErrMsgTxt
     $         ( 'TRANS has 0, or 1 the only admissible values' )
         END IF
      END IF
C
C   schur
C
      IF ( TASK.EQ.1 ) THEN
         NSCHUR = 1
         IF ( NRHS.EQ.IP+1 ) THEN
            ISIZE = mxGetM( PRHS(IP+1) ) * mxGetN( PRHS(IP+1) )
            IF ( ISIZE.GT.1 ) THEN
               CALL mexErrMsgTxt( 'SCHUR must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP+1) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP+1) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'SCHUR must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP+1) ), TEMP, ISIZE )
            IF ( ISIZE.GT.0 ) NSCHUR = TEMP
            IF ( NSCHUR.LT.1 .OR. NSCHUR.GT.2 ) THEN
               CALL mexErrMsgTxt
     $              ( 'SCHUR has 1, or 2 the only admissible values' )
            END IF
         END IF
      ELSE
         NSCHUR = 0
      END IF
C
C   A(NxN), (B(MxM),) C(NxM), or C(NxN), or op(C)(PxN)
C
      N = mxGetM( PRHS(2) )
      IF ( TASK.EQ.1 ) THEN
         M = mxGetM( PRHS(3) )
      ELSE IF ( TASK.EQ.3 ) THEN
         IF ( TRANS.EQ.0 ) THEN
            P = mxGetM( PRHS(3) )
         ELSE
            P = mxGetN( PRHS(3) )
         END IF
      END IF
C
      IF ( mxGetN( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
      IP = 3
      IF ( TASK.EQ.1 ) THEN
C
         IF ( mxGetN( PRHS(IP) ).NE.M ) THEN
            CALL mexErrMsgTxt( 'B must be a square matrix' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'B must be a real matrix' )
         END IF
         IP = IP + 1
      END IF
C
      MC = mxGetM( PRHS(IP) )
      NC = mxGetN( PRHS(IP) )
C
      IF ( TASK.NE.3 ) THEN
C
         IF ( MC.NE.N ) THEN
            CALL mexErrMsgTxt
     $           ( 'C must have the same row dimension as A' )
         END IF
         IF ( TASK.EQ.1 ) THEN
            IF ( NC.NE.M )
     $         CALL mexErrMsgTxt
     $              ( 'C must have the same column dimension as B' )
         ELSE
            IF ( NC.NE.N )
     $         CALL mexErrMsgTxt
     $              ( 'C must have the same column dimension as A' )
         END IF
C
      ELSE
C
         IF ( TRANS.EQ.0 ) THEN
            IF ( NC.NE.N )
     $         CALL mexErrMsgTxt
     $              ( 'C must have the same column dimension as A' )
         ELSE
            IF ( MC.NE.N )
     $         CALL mexErrMsgTxt
     $              ( 'C must have the same row dimension as A' )
         END IF
C
      END IF
C
      IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
      IP = IP + 1
C
      FLAG(1) = 0
      FLAG(2) = 0
      FLAG(3) = 0
      IF ( NRHS.GE.IP ) THEN
C
C   flag
C
         ISIZE = mxGetM( PRHS(IP) ) * mxGetN( PRHS(IP) )
         IF ( TASK.EQ.1 ) THEN
            IF ( ISIZE.GT.3 )
     $         CALL mexErrMsgTxt
     $              ( 'FLAG must be a vector with at most 3 elements' )
         ELSE
            IF ( ISIZE.GT.2 )
     $         CALL mexErrMsgTxt
     $              ( 'FLAG must be a vector with at most 2 elements' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'FLAG must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), FLAGR, ISIZE )
         IF ( ISIZE.GT.0 ) FLAG(1) = FLAGR(1)
         IF ( ISIZE.GT.1 ) FLAG(2) = FLAGR(2)
         IF ( ISIZE.GT.2 ) FLAG(3) = FLAGR(3)
      END IF
      IF ( FLAG(2).NE.0 .AND. FLAG(2).NE.1 ) THEN
         IF ( TASK.NE.1 .OR. ( TASK.EQ.1 .AND. FLAG(2).NE.2 ) )
     $      FLAG(2) = 0
      END IF
      IF ( TASK.EQ.1 .AND. ( FLAG(3).LT.0 .OR. FLAG(3).GT.2 ) )
     $   FLAG(3) = 0
C
C Determine the lenghts of working arrays.
C Use a value for LDWORK enabling calls of block algorithms
C in DGEES, and possibly in DGEHRD, DGEQRF, DGERQF, SB04PD.
C
      LDA = MAX( 1, N )
      IF ( NSCHUR.NE.1 ) THEN
         IP = N
      ELSE
         IP = M
      END IF
      IF ( FLAG(2).NE.1 ) THEN
         IF ( TASK.NE.1 .OR. ( TASK.EQ.1 .AND. NSCHUR.EQ.2 ) ) THEN
            CALL DGEES( 'Vectors', 'No sort', SELECT, IP, DUM,
     $                  MAX( 1, IP ), IP, DUM, DUM, DUM, MAX( 1, IP ),
     $                  DUM, -1, BWORK, INFO )
            LDWO = DUM(1)
         ELSE
            LDWO = 0
         END IF
      ELSE
         LDWO = 0
      END IF
      IF ( TASK.EQ.1 ) THEN
         LDB = MAX( 1, M )
         IF ( NSCHUR.EQ.2 ) THEN
C
C           Need    LDWORK = MAX( 1, LDW1 + LDW2 ) with LDWO = 0.
C           Prefer  larger.
C
            IF ( FLAG(2).EQ.1 ) THEN
               LDW1 = 0
               LDW2 = 0
            ELSE
               LDW1 = 3*N + 1
               LDW2 = 2*N + LDWO
            END IF
            IB = 0
            IF ( FLAG(3).NE.1 ) THEN
               IB = 2*M
               IF ( FLAG(2).EQ.1 )
     $            IB = IB + 1
               IF ( FLAG(3).NE.2 ) THEN
                  CALL DGEES( 'Vectors', 'No sort', SELECT, M, DUM, LDB,
     $                        IP, DUM, DUM, DUM, LDB, DUM, -1, BWORK,
     $                        INFO )
                  LDWO = DUM(1)
               END IF
               LDW2 = MAX( LDW2, MAX( IB + M, LDWO ) + 2*M )
            END IF
            LDW2   = MAX( LDW2, IB + 2*N )
            LDWORK = MAX( 1, LDW1 + LDW2 )
            IF ( FLAG(2)*FLAG(3).NE.1 ) THEN
     $         LDWORK = MAX( LDWORK, LDW1 - N + IB + M*N )
         ELSE
            IF ( FLAG(2)*FLAG(3).EQ.1 ) THEN
               LIWORK = 0
               IF ( FLAG(1).NE.0 ) THEN
                  LDWORK = 2*N
               ELSE
                  LDWORK = 0
               END IF
            ELSE IF ( FLAG(2)*FLAG(3).EQ.2 ) THEN
               MXMN   = MAX( M, N )
               LIWORK = 2*MXMN
               LDWORK = 2*MXMN*( 4 + 2*MXMN )
            ELSE
               LIWORK = 4*N
               LDWORK = MAX( 1, 5*M, N + M )
               IF ( FLAG(1).EQ.0 ) THEN
                  LDWORK = MAX( LDWORK, 2*N*N + 8*N )
               ELSE
                  LDWORK = MAX( LDWORK, 2*N*N + 9*N )
               END IF
C
C              Need    LDWORK above.
C              Prefer  LDWORK below.
C
               LDWORK = MAX( LDWORK, 2*M + LDWO, N + N*M, N +
     $                       N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 ),
     $                       N +
     $                       M*ILAENV( 1, 'DORMQR', 'LT', N, M, N, -1 ),
     $                       N +
     $                       M*ILAENV( 1, 'DORMQR', 'LN', N, M, N, -1 )
     $                     )
            END IF
         END IF
         NM = M
C
      ELSE IF ( TASK.EQ.2 ) THEN
C
C        Need    LDWORK = MAX( 1, N*N, 3*N )
C        Prefer  larger.
C
         LDWORK = MAX( 1, N*N, 3*N, LDWO )
         IF ( NLHS.EQ.2 ) THEN
            LDWORK = MAX( LDWORK, 2*N*N )
            IF ( FLAG(1).NE.0 )
     $         LDWORK = MAX( LDWORK, 2*( N*N + N ) )
         END IF
         NM = N
      END IF
      IF ( TASK.NE.3 ) THEN
         LDC = LDA
      ELSE
         IF ( TRANS.EQ.0 ) THEN
            LDC = MAX( 1, N, P )
         ELSE
            LDC = LDA
         END IF
         MXMN = MIN( P, N )
C
C        Need    LDWORK = MAX( 1, 4*N + MXMN )
C        Prefer  larger.
C
         IF ( TRANS.EQ.0 ) THEN
            NB1 = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
            NB2 = ILAENV( 1, 'DGEQRF', ' ', N, N, -1, -1 )
         ELSE
            NB1 = ILAENV( 1, 'DGERQF', ' ', N, M, -1, -1 )
            NB2 = ILAENV( 1, 'DGERQF', ' ', N, N, -1, -1 )
         END IF
         LDWORK = MAX( 1, MAX( LDWO, MXMN + MAX( N*NB1, 4*N ), N*MXMN,
     $                         N*N, N + N*NB2 ) )
         NM = N
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      IF ( TASK.EQ.1 ) THEN
         ALLOCATE ( A( LDA, N ), B( LDB, M ), C( LDC, M ) )
         ALLOCATE ( DWORK( LDWORK ) )
         IF ( NSCHUR.EQ.2 ) THEN
            IF ( FLAG(2).EQ.1 ) THEN
               FACTA = 'S'
               LDU   = 1
               ALLOCATE ( U( LDU, 1 ) )
            ELSE
               FACTA = 'N'
               LDU   = LDA
               ALLOCATE ( U( LDU, N ) )
            END IF
            IF ( FLAG(3).EQ.1 ) THEN
               FACTB = 'S'
               LDV   = 1
               ALLOCATE ( V( LDV, 1 ) )
            ELSE
               FACTB = 'N'
               LDV   = LDB
               ALLOCATE ( V( LDV, M ) )
            END IF
         ELSE
            SCHU = 'N'
            ALLOCATE ( IWORK( LIWORK ) )
            IF ( FLAG(2).EQ.1 ) THEN
               IF ( FLAG(3).EQ.1 ) THEN
                  SCHU = 'S'
               ELSE IF ( FLAG(3).EQ.2 ) THEN
                  SCHU = 'A'
               END IF
            ELSE IF ( FLAG(2).EQ.2 .AND. FLAG(3).EQ.1 ) THEN
               SCHU = 'B'
            END IF
C
            IF ( LSAME( SCHU, 'N' ) ) THEN
               LDU = LDB
               ALLOCATE ( U( LDU, M ) )
            END IF
         END IF
C
      ELSE IF ( TASK.EQ.2 ) THEN
         LDU = LDA
         ALLOCATE ( A( LDA, N ), C( LDC, N ), U( LDU, N ) )
         ALLOCATE ( DWORK( LDWORK ), WI( N ), WR( N ) )
         IF ( NLHS.EQ.1 ) THEN
            ALLOCATE ( IWORK( 1 ) )
         ELSE
            ALLOCATE ( IWORK( N*N ) )
         END IF
      ELSE
         IF ( TRANS.EQ.0 ) THEN
            ALLOCATE ( A( LDA, N ), C( LDC, N ) )
         ELSE
            ALLOCATE ( A( LDA, N ), C( LDC, MAX( P, N ) ) )
         END IF
         LDU = LDA
         ALLOCATE ( U( LDU, N ), DWORK( LDWORK ), WI( N ), WR( N ) )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      IF ( TASK.EQ.3 ) THEN
         IF ( P.GE.N .OR. TRANS.NE.0 ) THEN
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), C, P*N )
         ELSE
C
C           P < N and trans = 0. Use A temporarily for loading C.
C
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), A, P*N )
            CALL DLACPY( 'Full', P, N, A, MAX( 1, P ), C, LDC )
         END IF
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, N*N )
      IF ( TASK.EQ.1 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), B, M*M )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), C, N*M )
      ELSE IF ( TASK.EQ.2 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), C, N*N )
      END IF
C
C Do the actual computations.
C
      IF ( TASK.EQ.1 ) THEN
C
         IF ( NSCHUR.EQ.2 ) THEN
            IF ( TRANS.EQ.0 ) THEN
               TRANA = 'N'
               TRANB = 'N'
            ELSE IF ( TRANS.EQ.1 ) THEN
               TRANA = 'T'
               TRANB = 'T'
            ELSE IF ( TRANS.EQ.2 ) THEN
               TRANA = 'T'
               TRANB = 'N'
            ELSE IF ( TRANS.EQ.3 ) THEN
               TRANA = 'N'
               TRANB = 'T'
            END IF
            IF ( FLAG(1).NE.0 ) THEN
               DICO = 'D'
            ELSE
               DICO = 'C'
            END IF
            ISGN = 1
            CALL SB04PD( DICO, FACTA, FACTB, TRANA, TRANB, ISGN, N, M,
     $                   A, LDA, U, LDU, B, LDB, V, LDV, C, LDC, SCALE,
     $                   DWORK, LDWORK, INFO )
         ELSE
            IF ( TRANS.EQ.0 ) THEN
               IF ( LSAME( SCHU, 'S' ) ) THEN
                  TRANA = 'N'
                  TRANB = 'N'
               ELSE
                  ULA = 'U'
                  ULB = 'U'
               END IF
            ELSE IF ( TRANS.EQ.1 ) THEN
               IF ( LSAME( SCHU, 'S' ) ) THEN
                  TRANA = 'T'
                  TRANB = 'T'
               ELSE
                  ULA = 'L'
                  ULB = 'L'
                  DO 10 J = 2, N
                     CALL DSWAP( J-1, A( 1, J ), 1, A( J, 1 ), LDA )
   10             CONTINUE
                  DO 20 J = 2, M
                     CALL DSWAP( J-1, B( 1, J ), 1, B( J, 1 ), LDB )
   20             CONTINUE
               END IF
            ELSE IF ( TRANS.EQ.2 ) THEN
               IF ( LSAME( SCHU, 'S' ) ) THEN
                  TRANA = 'T'
                  TRANB = 'N'
               ELSE
                  ULA = 'L'
                  ULB = 'U'
                  DO 30 J = 2, N
                     CALL DSWAP( J-1, A( 1, J ), 1, A( J, 1 ), LDA )
   30             CONTINUE
               END IF
            ELSE IF ( TRANS.EQ.3 ) THEN
               IF ( LSAME( SCHU, 'S' ) ) THEN
                  TRANA = 'N'
                  TRANB = 'T'
               ELSE
                  ULA = 'U'
                  ULB = 'L'
                  DO 40 J = 2, M
                     CALL DSWAP( J-1, B( 1, J ), 1, B( J, 1 ), LDB )
   40             CONTINUE
               END IF
            END IF
C
            IF ( LSAME( SCHU, 'N' ) ) THEN
C
               SCALE = ONE
               IF ( FLAG(1).EQ.0 ) THEN
                  CALL SB04MD( N, M, A, LDA, B, LDB, C, LDC, U, LDU,
     $                         IWORK, DWORK, LDWORK, INFO )
               ELSE
                  CALL SB04QD( N, M, A, LDA, B, LDB, C, LDC, U, LDU,
     $                         IWORK, DWORK, LDWORK, INFO )
               END IF
C
            ELSE IF ( LSAME( SCHU, 'S' ) ) THEN
C
               IF ( FLAG(1).EQ.0 ) THEN
                  CALL DTRSYL( TRANA, TRANB, 1, N, M, A, LDA, B, LDB, C,
     $                         LDC, SCALE, INFO )
                  IF ( MIN( N, M ).EQ.0 )
     $               SCALE = ONE
               ELSE
                  CALL SB04PY( TRANA, TRANB, 1, N, M, A, LDA, B, LDB, C,
     $                         LDC, SCALE, DWORK, INFO )
               END IF
            ELSE
               SCALE = ONE
               TOL   = ZERO
C
C              Default tolerance (epsilon_machine) is used.
C
               IF ( FLAG(1).EQ.0 ) THEN
                  CALL SB04ND( SCHU, ULA, ULB, N, M, A, LDA, B, LDB, C,
     $                         LDC, TOL, IWORK, DWORK, LDWORK, INFO )
               ELSE
                  CALL SB04RD( SCHU, ULA, ULB, N, M, A, LDA, B, LDB, C,
     $                         LDC, TOL, IWORK, DWORK, LDWORK, INFO )
               END IF
C
            END IF
         END IF
C
      ELSE
C
         IF ( FLAG(1).EQ.0 ) THEN
            DICO = 'C'
         ELSE
            DICO = 'D'
         END IF
         IF ( FLAG(2).NE.1 ) THEN
            FACT = 'N'
         ELSE
            FACT = 'F'
            CALL DLASET( 'Full', N, N, ZERO, ONE, U, LDU )
         END IF
C
         IF ( TRANS.EQ.0 ) THEN
            TRANA = 'N'
         ELSE
            TRANA = 'T'
         END IF
C
         IF ( TASK.EQ.2 ) THEN
            IF ( NLHS.EQ.2 ) THEN
               JOB = 'B'
            ELSE
               JOB = 'X'
            END IF
C
            CALL SB03MD( DICO, JOB, FACT, TRANA, N, A, LDA, U, LDU,
     $                   C, LDC, SCALE, SEP, FERR, WR, WI, IWORK,
     $                   DWORK, LDWORK, INFO )
C
         ELSE
C
            CALL SB03OD( DICO, FACT, TRANA, N, P, A, LDA, U, LDU,
     $                   C, LDC, SCALE, WR, WI, DWORK, LDWORK, INFO )
C
         END IF
      END IF
C
C Copy output to MATLAB workspace.
C
      PERTRB = ( TASK.EQ.1 .AND. ( INFO.EQ.N+M+1 .OR.
     $         ( FLAG(2)*FLAG(3).EQ.1 .AND. INFO.EQ.1 ) ) )
     $    .OR. ( TASK.EQ.2 .AND. INFO.EQ.N+1 )
     $    .OR. ( TASK.EQ.3 .AND. INFO.EQ.1 )
      IF ( INFO.EQ.0 .OR. PERTRB ) THEN
         IF ( NLHS.GE.1 ) THEN
            IF ( TASK.EQ.3 ) THEN
               IF ( TRANS.EQ.0 .AND. P.GT.N )
     $            CALL DLACPY( 'Upper', N, N, C, LDC, C, LDA )
               IF ( N.GT.1 )
     $            CALL DLASET( 'Lower', N-1, N-1, ZERO, ZERO, C(2,1),
     $                         LDA )
            END IF
            PLHS(1) = mxCreateDoubleMatrix( N, NM, 0 )
            CALL mxCopyReal8ToPtr( C, mxGetPr( PLHS(1) ), N*NM )
         END IF
C
         IF ( TASK.EQ.2 ) THEN
            IF ( NLHS.GE.2 ) THEN
               IF ( N.EQ.ZERO )
     $            SEP = ZERO
               PLHS(2) = mxCreateDoubleMatrix( 1, 1, 0 )
               CALL mxCopyReal8ToPtr( SEP, mxGetPr( PLHS(2) ), 1 )
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      IF ( TASK.EQ.1 ) THEN
         DEALLOCATE( A, B, C, DWORK )
         IF ( NSCHUR.EQ.2 ) THEN
            DEALLOCATE ( U, V )
         ELSE
            DEALLOCATE( IWORK )
            IF ( LSAME( SCHU, 'N' ) )
     $         DEALLOCATE ( U )
         END IF
      ELSE IF ( TASK.EQ.2 ) THEN
         DEALLOCATE( A, C, U, WR, WI, IWORK, DWORK )
      ELSE
         DEALLOCATE( A, C, U, WR, WI, DWORK )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( TASK.EQ.1 ) THEN
            IF ( NSCHUR.EQ.2 ) THEN
               WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB04PD'')'
     $              ) INFO
            ELSE
               IF ( LSAME( SCHU, 'N' ) ) THEN
                  IF ( FLAG(1).EQ.0 ) THEN
                     WRITE( TEXT, '('' INFO = '',I4,
     $                      '' ON EXIT FROM SB04MD'')' ) INFO
                  ELSE
                     WRITE( TEXT, '('' INFO = '',I4,
     $                      '' ON EXIT FROM SB04QD'')' ) INFO
                  END IF
               ELSE IF ( LSAME( SCHU, 'S' ) ) THEN
                  IF ( FLAG(1).EQ.0 ) THEN
                     WRITE( TEXT, '('' INFO = '',I4,
     $                      '' ON EXIT FROM DTRSYL'')' ) INFO
                  ELSE
                     WRITE( TEXT, '('' INFO = '',I4,
     $                      '' ON EXIT FROM SB04PY'')' ) INFO
                  END IF
               ELSE
                  IF ( FLAG(1).EQ.0 ) THEN
                     WRITE( TEXT, '('' INFO = '',I4,
     $                      '' ON EXIT FROM SB04ND'')' ) INFO
                  ELSE
                     WRITE( TEXT, '('' INFO = '',I4,
     $                      '' ON EXIT FROM SB04RD'')' ) INFO
                  END IF
               END IF
            END IF
         ELSE IF ( TASK.EQ.2 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB03MD'')'
     $           ) INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB03OD'')' )
     $             INFO
         END IF
      END IF
C
      IF ( ( INFO.EQ.0 .OR. PERTRB ) .AND. SCALE.NE.ONE ) THEN
         IF ( TASK.LE.2 ) THEN
            TEMP = SCALE
         ELSE
            TEMP = SCALE**2
         END IF
         WRITE ( TEXT, '( '' Warning: The right hand sides were '',
     $      ''scaled by'', D13.6, '' to avoid overflow. '' )' ) TEMP
      END IF
C
      IF ( INFO.NE.0 .AND. .NOT.PERTRB ) THEN
         CALL mexErrMsgTxt( TEXT )
      ELSE IF ( SCALE.NE.ONE ) THEN
         CALL mexPrintf( TEXT )
      END IF
      IF ( PERTRB ) THEN
         WRITE ( TEXT, '( '' Warning: The equation is (almost) '',
     $      ''singular; perturbed values have been used. '' )' )
         CALL mexPrintf( TEXT )
      END IF
C
      RETURN
C *** Last line of LINMEQ ***
      END
