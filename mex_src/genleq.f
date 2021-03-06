C GENLEQ.F - Gateway function for solving generalized linear matrix
C            equations using SLICOT routines SB04OD, SG03AD, and SG03BD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [X,(Y,)dif] = genleq(task,A(,D,B),E,C,(F,)flag,trans)
C
C   eq. (1) :  [X,Y,dif] = genleq(1,A,D,B,E,C,F,flag,trans)
C                  [X,Y] = genleq(1,A,D,B,E,C,F,flag,trans)
C   eq. (2) :    [X,dif] = genleq(2,A,D,B,E,C,flag,trans)
C                    [X] = genleq(2,A,D,B,E,C,flag,trans)
C   eq. (3) :    [X,sep] = genleq(3,A,E,C,flag,trans)
C                    [X] = genleq(3,A,E,C,flag,trans)
C   eq. (4) :        [X] = genleq(4,A,E,C,flag,trans)
C
C
C Purpose:
C   To solve the generalized linear matrix equations
C
C   task = 1:
C
C        (  A*X - Y*B = C,
C        <                                                         (1.a)
C        (  D*X - Y*E = F,
C
C        (  A'*X + D'*Y = C,
C        <                                                         (1.b)
C        (  X*B' + Y*E' = -F,
C
C   task = 2:
C
C       A*X*B - D*X*E = C,                                           (2)
C
C   task = 3:
C
C       op(A)'*X*op(E) + op(E)'*X*op(A) = C,                       (3.a)
C
C       op(A)'*X*op(A) - op(E)'*X*op(E) = C,                       (3.b)
C
C   task = 4:
C
C       op(A)'*op(X)'*op(X)*op(E) + op(E)'*op(X)'*op(X)*op(A) =
C                                               -op(C)'*op(C),     (4.a)
C
C       op(A)'*op(X)'*op(X)*op(A) - op(E)'*op(X)'*op(X)*op(E) =
C                                               -op(C)'*op(C),     (4.b)
C
C   where op(M) = M, if trans = 0, and op(M) = M', if trans <> 0, and
C   X in the equations (4) is upper triangular.
C
C Input parameters:
C   task  - integer option to determine the equation type:
C           = 1 : solve the linear matrix equation pairs (1.a) or (1.b);
C           = 2 : solve the linear matrix equation (2);
C           = 3 : solve the linear matrix equation (3.a) or (3.b);
C           = 4 : solve the linear matrix equation (4.a) or (4.b).
C   A, D  - real coefficient N-by-N matrices.
C   B, E  - real coefficient M-by-M matrices, for task = 1 or 2.
C   E     - real coefficient N-by-N matrix,   for task = 3 or 4.
C   C(,F) - right hand side N-by-M matrices (F exists only for (1)),
C           for task = 1 or 2.
C   C     - real coefficient N-by-N symmetric matrix, for task = 3.
C         - real coefficient M-by-N matrix, for task = 4 and trans =  0.
C         - real coefficient N-by-M matrix, for task = 4 and trans <> 0.
C   flag  - (optional) integer vector of length 2 containing options
C           task = 1, 2 :
C                flag(1) = 1 : (A,D) is in generalized Schur form;
C                              otherwise, (A,D) is in general form.
C                flag(2) = 1 : (B,E) is in generalized Schur form;
C                              otherwise, (B,E) is in general form.
C           task = 3, 4 :
C                flag(1) = 0 : solve the continuous-time equation (x.a);
C                              otherwise, solve the discrete-time
C                              equation (x.b), where x is 3 or 4.
C                flag(2) = 1 : (A,E) is in generalized Schur form;
C                              otherwise, (A,E) is in general form.
C           Default:      flag(1) = 0, flag(2) = 0.
C   trans - (optional) integer specifying a transposition option
C           trans = 0 : solve the equations (1.a), (2); or solve (3)
C                       or (4) with op(M) = M.
C                       otherwise, solve the "transposed" equations
C                       (1.b); or solve (3) or (4) with op(M) = M'.
C           Default:      trans = 0.
C
C Output parameters:
C   X(,Y) - solution of the equations (Y appears only for (1)).
C   dif   - (optional) estimator of Dif[(A,D),(B,E)].
C           dif is not computed for (1.b), i.e., for task = 1 with
C           trans <> 0.
C   sep   - (optional) estimator of sep(A,E).
C
C Comments:
C   1. Currently there is no available SLICOT routine for equations (2).
C   2. For task = 4, the pencil (A,E) must be stable, i.e., all
C      eigenvalues must have negative real parts, for (4.a), or moduli
C      strictly less than 1, for (4.b).
C
C Contributor:
C   H. Xu, TU Chemnitz, FR Germany, Dec. 1998.
C   V. Sima, Katholieke Univ. Leuven, Belgium, May 1999.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, April 2006,
C   April 2009, Dec. 2012, Apr. 2014.
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
      INTEGER           mxCreateDoubleMatrix, mxGetM, mxGetN, mxGetPr
      INTEGER*4         mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         DICO, FACT, JOB, SCHU, TRANS, UPLO
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDE, LDF, LDQ, LDU,
     $                  LDV, LDW, LDWORK, LDZ, M, N
      DOUBLE PRECISION  DIF, FERR, SCALE, SEP
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), ALPHAI(:), ALPHAR(:),
     $                                 BETA(:), B(:,:), C(:,:), D(:,:),
     $                                 DWORK(:), E(:,:), F(:,:), Q(:,:),
     $                                 U(:,:), V(:,:), W(:,:), Z(:,:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           FLAG(2), IDIF, IP, ISIZE, LDWO, NB1, NB2, TASK,
     $                  TRANE
      DOUBLE PRECISION  DUM(1), FLAGR(2), TEMP
C
C .. External subroutines ..
      INTEGER           ILAENV
      EXTERNAL          DGEES, DLACPY, DLASET, ILAENV, SB04OD, SG03AD,
     $                  SG03BD
C
C ..Intrinsic functions..
      INTRINSIC         MAX
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'GENLEQ requires at least 4 input arguments' )
      ELSE IF ( NLHS.GT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'GENLEQ requires at most 3 output arguments' )
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
      IF ( TASK.LT.1 .OR. TASK.GT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'TASK has 1 2, 3, or 4 the only admissible values' )
      END IF
C
      IF ( TASK.EQ.1 ) THEN
         IF ( NRHS.LT.7 ) THEN
            CALL mexErrMsgTxt(
     $          'TASK = 1: GENLEQ requires at least 7 input arguments' )
         ELSE
            IP = 9
         END IF
      ELSE IF ( TASK.EQ.2 ) THEN
         IF ( NRHS.LT.6 ) THEN
            CALL mexErrMsgTxt(
     $          'TASK = 2: GENLEQ requires at least 6 input arguments' )
         ELSE
            IP = 8
         END IF
      ELSE
         IP = 6
      END IF
C
C   trans
C
      TRANE = 0
      IF ( NRHS.GE.IP ) THEN
         ISIZE = mxGetM( PRHS(IP) ) * mxGetN( PRHS(IP) )
         IF ( ISIZE.GT.1 )
     $      CALL mexErrMsgTxt( 'TRANS must be a scalar' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TRANS must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IF ( ISIZE.GT.0 ) TRANE = TEMP
      END IF
C
C   A(NxN), D(NxN), B(MxM), E(MxM), C(NxM), F(NxM), for task = 1, 2;
C   A(NxN), E(NxN), C(NxN),                         for task = 3;
C   A(NxN), E(NxN), op(C)(MxN),                     for task = 4.
C
      N = mxGetM( PRHS(2) )
      IF ( TASK.LE.2 ) THEN
         M = mxGetM( PRHS(4) )
      ELSE IF ( TASK.EQ.4 ) THEN
         IF ( TRANE.EQ.0 ) THEN
            M = mxGetM( PRHS(4) )
         ELSE
            M = mxGetN( PRHS(4) )
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
      IF ( TASK.LE.2 ) THEN
         IF ( mxGetM( PRHS(3) ).NE.N .OR. mxGetN( PRHS(3) ).NE.N ) THEN
            CALL mexErrMsgTxt( 'D must have the same size as A' )
         END IF
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(3) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'D must be a real matrix' )
         END IF
C
         IF ( mxGetN( PRHS(4) ).NE.M ) THEN
            CALL mexErrMsgTxt( 'B must be a square matrix' )
         END IF
         IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(4) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'B must be a real matrix' )
         END IF
C
         IP = 5
         IF ( mxGetM( PRHS(IP) ).NE.M .OR. mxGetN( PRHS(IP) ).NE.M )
     $         THEN
            CALL mexErrMsgTxt( 'E must have the same size as B' )
         END IF
C
      ELSE
C
         IP = 3
         IF ( mxGetM( PRHS(IP) ).NE.N .OR. mxGetN( PRHS(IP) ).NE.N )
     $         THEN
            CALL mexErrMsgTxt( 'E must have the same size as A' )
         END IF
      END IF
C
      IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'E must be a real matrix' )
      END IF
      IP = IP + 1
C
      IF ( TASK.LE.2 ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $           ( 'C must have the same row dimension as A' )
         END IF
         IF ( mxGetN( PRHS(IP) ).NE.M ) THEN
            CALL mexErrMsgTxt
     $           ( 'C must have the same column dimension as B' )
         END IF
C
      ELSE IF ( TASK.EQ.3 ) THEN
C
         IF ( mxGetM( PRHS(IP) ).NE.N .OR. mxGetN( PRHS(IP) ).NE.N )
     $         THEN
            CALL mexErrMsgTxt( 'C must have the same size as A' )
         END IF
C
      ELSE IF ( TRANE.EQ.0 ) THEN
         IF ( mxGetN( PRHS(IP) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $           ( 'C must have the same column dimension as A' )
         END IF
      ELSE
         IF ( mxGetM( PRHS(IP) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $           ( 'C must have the same row dimension as A' )
         END IF
      END IF
C
      IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
      IP = IP + 1
C
      IF ( TASK.EQ.1 ) THEN
         IF ( mxGetM( PRHS(7) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $           ( 'F must have the same row dimension as A' )
         END IF
         IF ( mxGetN( PRHS(7) ).NE.M ) THEN
            CALL mexErrMsgTxt
     $           ( 'F must have the same column dimension as B' )
         END IF
         IF ( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(7) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'F must be a real matrix' )
         END IF
         IP = IP + 1
      END IF
C
C   flag(2x1)
C
      FLAG(1) = 0
      FLAG(2) = 0
      IF ( NRHS.GE.IP ) THEN
         ISIZE = mxGetM( PRHS(IP) ) * mxGetN( PRHS(IP) )
         IF ( ISIZE.GT.2 )
     $      CALL mexErrMsgTxt
     $           ( 'FLAG must be a vector with at most 2 elements' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'FLAG must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), FLAGR, ISIZE )
         IF ( ISIZE.GT.0 ) FLAG(1) = FLAGR(1)
         IF ( ISIZE.GT.1 ) FLAG(2) = FLAGR(2)
      END IF
C
C Determine the lenghts of working arrays.
C Use a value for LDWORK enabling calls of block algorithms, when
C possible.
C
      LDA = MAX( 1, N )
      LDQ = LDA
      IF ( TASK.LE.2 ) THEN
         LDB = MAX( 1, M )
         LDC = LDA
         LDD = LDA
         LDE = LDB
         LDF = LDA
         LDU = LDA
         LDV = LDB
         LDW = LDB
      ELSE
         LDE = LDA
         LDZ = LDA
         IF ( TASK.EQ.3 ) THEN
            LDC = LDA
         ELSE IF ( TRANE.EQ.0 ) THEN
            LDC = MAX( LDA, M )
         ELSE
            LDC = LDA
         END IF
      END IF
C
      IF ( TRANE.EQ.0 ) THEN
         TRANS = 'N'
      ELSE
         TRANS = 'T'
      END IF
C
      IF ( TASK.LE.2 ) THEN
         SCHU = 'R'
         IF ( FLAG(1).EQ.1 ) THEN
            IF ( FLAG(2).EQ.1 ) THEN
               SCHU = 'N'
            ELSE
               SCHU = 'B'
            END IF
         ELSE IF ( FLAG(2).EQ.1 ) THEN
            SCHU = 'A'
         END IF
C
         IF ( NLHS.GE.IDIF ) THEN
            JOB = 'D'
         ELSE
            JOB = 'N'
         END IF
      END IF
C
      IDIF = 4 - TASK
      IF ( TASK.EQ.1 ) THEN
         CALL SB04OD( SCHU, TRANS, JOB, N, M, A, LDA, B, LDB, C, LDC, D,
     $                LDD, E, LDE, F, LDF, SCALE, DIF, Q, LDQ, U, LDU,
     $                V, LDV, W, LDW, IWORK, DUM, -1, INFO )
         LDWORK = DUM(1)
      ELSE IF ( TASK.GT.2 .AND. FLAG(2).NE.1 ) THEN
         CALL DGEES( 'Vectors', 'Vectors', N, DUM, LDA, DUM, LDA, DUM,
     $               DUM, DUM, DUM, LDA, DUM, LDA, DUM, -1, INFO )
         LDWO = DUM(1)
      ELSE
         LDWO = 0
      END IF
      IF ( TASK.EQ.3 ) THEN
C
C        Need    LDWORK is computed as follows:
C                LDWORK = MAX( 1, N ),     if flag(2) = 1;
C                LDWORK = MAX( 1, 4*N ),   if flag(2) <> 1;
C                LDWORK = MAX( LDWORK, 2*N*N ), if sep is required.
C        Prefer  larger.
C
         IF ( FLAG(2).EQ.1 ) THEN
            LDWORK = MAX( 1, N )
         ELSE
            LDWORK = MAX( 1, 4*N, LDWO )
         END IF
         IF ( NLHS.NE.1 )
     $      LDWORK = MAX( LDWORK, 2*N*N )
         LDWORK = MAX( LDWORK, N*N )
      ELSE IF ( TASK.EQ.4 ) THEN
C
C        Need    LDWORK is computed as follows:
C                LDWORK = MAX( 1, 2*N, 6*N - 6 ),   if flag(2) = 1;
C                LDWORK = MAX( 1, 4*N, 6*N - 6 ),   if flag(2) <> 1.
C        Prefer  larger.
C
         IF ( FLAG(2).EQ.1 ) THEN
            LDWORK = MAX( 1, 2*N )
         ELSE
            LDWORK = MAX( 1, 4*N, LDWO )
         END IF
         IF ( TRANE.EQ.0 ) THEN
            NB1 = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
            NB2 = ILAENV( 1, 'DGEQRF', ' ', N, N, -1, -1 )
         ELSE
            NB1 = ILAENV( 1, 'DGERQF', ' ', N, M, -1, -1 )
            NB2 = ILAENV( 1, 'DGERQF', ' ', N, N, -1, -1 )
         END IF
         LDWORK = MAX( LDWORK, 6*N - 6, N + N*MAX( M, N, NB1, NB2 ) )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      IF ( TASK.EQ.1 ) THEN
         ALLOCATE( A( LDA, N ), D( LDD, N ), B( LDB, M ), E( LDE, M ),
     $             C( LDC, M ), F( LDF, M ), Q( LDQ, N ), U( LDU, N ),
     $             V( LDV, M ), W( LDW, M ) )
         ALLOCATE ( IWORK( M+N+6 ), DWORK( LDWORK ) )
      ELSE IF ( TASK.EQ.2 ) THEN
         ALLOCATE( A( LDA, N ), D( LDD, N ), B( LDB, M ), E( LDE, M ),
     $             C( LDC, M ), Q( LDQ, N ), U( LDU, N ),
     $             V( LDV, M ), W( LDW, M ) )
         ALLOCATE ( IWORK( M+N+6 ), DWORK( LDWORK ) )
      ELSE
         ALLOCATE( A( LDA, N ), E( LDE, N ), Q( LDQ, N ), Z( LDZ, N ),
     $             ALPHAI( N ), ALPHAR( N ), BETA( N ),
     $             DWORK( LDWORK ) )
         IF ( TASK.EQ.3 ) THEN
            ALLOCATE( C( LDC, N ) )
            IF ( NLHS.GE.2 ) THEN
               ALLOCATE ( IWORK( N*N ) )
            ELSE
               ALLOCATE ( IWORK( 1 ) )
            END IF
         ELSE IF ( TRANE.EQ.0 ) THEN
            ALLOCATE( C( LDC, N ) )
         ELSE
            ALLOCATE( C( LDC, MAX( M, N ) ) )
         END IF
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, N*N )
      IF ( TASK.LE.2 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), D, N*N )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), B, M*M )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), E, M*M )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), C, N*M )
         IF ( TASK.EQ.1 )
     $      CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), F, N*M )
      ELSE IF ( TASK.GE.3 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), E, N*N )
         IF ( TASK.EQ.3 ) THEN
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), C, N*N )
         ELSE
            IF ( TRANE.EQ.0 .AND. M.LT.N ) THEN
               CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), Q, M*N )
               CALL DLACPY( 'Full', M, N, Q, MAX( 1, M ), C, LDC )
            ELSE
               CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), C, M*N )
            END IF
         END IF
      END IF
C
C Do the actual computations.
C
      IF ( TASK.LE.2 ) THEN
         IF ( TASK.EQ.1 ) THEN
            CALL SB04OD( SCHU, TRANS, JOB, N, M, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, E, LDE, F, LDF, SCALE, DIF, Q, LDQ,
     $                   U, LDU, V, LDV, W, LDW, IWORK, DWORK, LDWORK,
     $                   INFO )
         ELSE
            INFO  = 0
            IDIF  = 3
            SCALE = ONE
            CALL mexPrintf( 'TASK = 2: Currently it is not available' )
         END IF
C
      ELSE IF ( TASK.GE.3 ) THEN
         IF ( FLAG(1).EQ.0 ) THEN
            DICO = 'C'
         ELSE
            DICO = 'D'
         END IF
         IF ( FLAG(2).EQ.1 ) THEN
            FACT = 'F'
            CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
            CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
         ELSE
            FACT = 'N'
         END IF
C
         IF ( TASK.EQ.3 ) THEN
            UPLO = 'U'
            IF ( NLHS.EQ.1 ) THEN
               JOB = 'X'
            ELSE
               JOB = 'B'
            END IF
C
            CALL SG03AD( DICO, JOB, FACT, TRANS, UPLO, N, A, LDA,
     $                   E, LDE, Q, LDQ, Z, LDZ, C, LDC, SCALE, SEP,
     $                   FERR, ALPHAR, ALPHAI, BETA, IWORK, DWORK,
     $                   LDWORK, INFO )
            IDIF = 2
         ELSE
            CALL SG03BD( DICO, FACT, TRANS, N, M, A, LDA, E, LDE,
     $                   Q, LDQ, Z, LDZ, C, LDC, SCALE, ALPHAR, ALPHAI,
     $                   BETA, DWORK, LDWORK, INFO )
         END IF
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         IF ( NLHS.GE.1 ) THEN
            IF ( TASK.LE.2 ) THEN
               PLHS(1) = mxCreateDoubleMatrix( N, M, 0 )
               CALL mxCopyReal8ToPtr( C, mxGetPr( PLHS(1) ), N*M )
            ELSE
               PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
               IF ( TASK.EQ.4 .AND. TRANE.EQ.0 .AND. M.GT.N )
     $            CALL DLACPY( 'Full', N, N, C, LDC, C, LDA )
               CALL mxCopyReal8ToPtr( C, mxGetPr( PLHS(1) ), N*N )
            END IF
         END IF
         IF ( NLHS.GE.2 ) THEN
            IF ( TASK.EQ.1 ) THEN
               PLHS(2) = mxCreateDoubleMatrix( N, M, 0 )
               CALL mxCopyReal8ToPtr( F, mxGetPr( PLHS(2) ), N*M )
            END IF
            IF ( NLHS.GE.IDIF ) THEN
               PLHS(IDIF) = mxCreateDoubleMatrix( 1, 1, 0 )
               IF ( TASK.EQ.1 ) THEN
                  CALL mxCopyReal8ToPtr( DIF, mxGetPr( PLHS(IDIF) ), 1 )
               ELSE IF ( TASK.EQ.3 ) THEN
                  CALL mxCopyReal8ToPtr( SEP, mxGetPr( PLHS(IDIF) ), 1 )
               END IF
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      IF ( TASK.EQ.1 ) THEN
         DEALLOCATE( A, B, C, D, E, F, Q, U, V, W, IWORK, DWORK )
      ELSE IF ( TASK.EQ.2 ) THEN
         DEALLOCATE( A, B, C, D, E, Q, U, V, W, IWORK, DWORK )
      ELSE IF ( TASK.EQ.3 ) THEN
         DEALLOCATE( A, C, E, Q, Z, IWORK, DWORK, ALPHAI, ALPHAR, BETA )
      ELSE
         DEALLOCATE( A, C, E, Q, Z, DWORK, ALPHAI, ALPHAR, BETA )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( TASK.EQ.1 ) THEN
            WRITE ( TEXT,'('' INFO = '',I4,'' ON EXIT FROM SB04OD'')' )
     $              INFO
         ELSE IF ( TASK.EQ.3 ) THEN
            WRITE ( TEXT,'('' INFO = '',I4,'' ON EXIT FROM SG03AD'')' )
     $              INFO
         ELSE IF ( TASK.EQ.4 ) THEN
            WRITE ( TEXT,'('' INFO = '',I4,'' ON EXIT FROM SG03BD'')' )
     $              INFO
         END IF
      END IF
C
      IF ( TASK.EQ.1 .OR. TASK.EQ.3 ) THEN
         TEMP = SCALE
      ELSE IF ( TASK.EQ.4 ) THEN
         TEMP = SCALE**2
      END IF
C
      IF ( INFO.EQ.0 .AND. SCALE.NE.ONE ) THEN
         WRITE ( TEXT, '( '' Warning: The right hand side(s) were '',
     $      ''scaled by'', D13.6, '' to avoid overflow'' )' ) TEMP
      END IF
C
      IF ( INFO.NE.0 ) THEN
         CALL mexErrMsgTxt( TEXT )
      ELSE IF ( SCALE.NE.ONE ) THEN
         CALL mexPrintf( TEXT )
      END IF
C
      RETURN
C *** Last line of GENLEQ ***
      END
