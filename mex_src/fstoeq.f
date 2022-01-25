C FSTOEQ.F - Gateway function for the QR factorization of (block)
C            Toeplitz matrices and/or solving associated linear
C            least-squares systems using SLICOT routines MB02HD,
C            MB02ID, MB02JD, MB02JX, MB02KD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C            [(Q,)R(,X,Y)] = fstoeq(task,TC,TR(,B,C))
C                  [(Q,)R] = fstoeq(task,TC,TR,tol)
C                [R(,X,Y)] = fstoeq(task,m,n,TCB,TRB(,B,C))
C                  [X(,Y)] = fstoeq(task,TC,TR,B(,C))
C                      [Y] = fstoeq(task,TC,TR,C)
C
C   task =  1 :        [R] = fstoeq( 1,TC,TR)
C                    [Q,R] = fstoeq( 1,TC,TR)
C   task =  2 :      [R,X] = fstoeq( 2,TC,TR,B)
C                  [Q,R,X] = fstoeq( 2,TC,TR,B)
C   task =  3 :      [R,Y] = fstoeq( 3,TC,TR,C)
C                  [Q,R,Y] = fstoeq( 3,TC,TR,C)
C   task =  4 :    [R,X,Y] = fstoeq( 4,TC,TR,B,C)
C                [Q,R,X,Y] = fstoeq( 4,TC,TR,B,C)
C   task =  5 :      [R,E] = fstoeq( 5,TC,TR,tol)
C                  [Q,R,E] = fstoeq( 5,TC,TR,tol)
C   task =  6 :        [R] = fstoeq( 6,m,n,TCB,TRB)
C   task =  7 :      [R,X] = fstoeq( 7,m,n,TCB,TRB,B)
C   task =  8 :      [R,Y] = fstoeq( 8,m,n,TCB,TRB,C)
C   task =  9 :    [R,X,Y] = fstoeq( 9,m,n,TCB,TRB,B,C)
C   task = 10 :        [X] = fstoeq(10,TC,TR,B)
C   task = 11 :        [Y] = fstoeq(11,TC,TR,C)
C   task = 12 :      [X,Y] = fstoeq(12,TC,TR,B,C)
C   task = 13 :        [X] = fstoeq(13,TC,TR,B)
C
C Purpose:
C   To compute an orthogonal-triangular/trapezoidal decomposition
C   (with partial column pivoting) of a (banded) block Toeplitz matrix.
C   If task <= 5 or task >= 10, the first block column of T is contained
C   in TC and has the dimension M*K-by-L and the first block row of T is
C   contained in TR and has the dimension K-by-N*L, while for task >= 6
C   and task <= 9 the leading nonzero blocks of the first block column
C   of T are contained in TCB which has the dimension (ML+1)*K-by-L
C   and the leading nonzero blocks of the first block row of T are
C   contained in TRB which has the dimesion K-by-(NU+1)*L.
C   If task = 13, the product of a block Toeplitz matrix T with a block
C   column vector C is computed.
C
C   task =  1: Compute R, the Cholesky factor of T'*T, i.e., R is lower
C              triangular, so that R*R' = T'*T. Optionally, a matrix Q
C              is computed so that Q'*Q = I and T = Q*R'.
C
C   task =  2: Compute R and the least-squares solution of T*X = B for a
C              given right hand side matrix B.
C
C   task =  3: Compute R and the minimum norm solution of T'*Y = C for a
C              given right hand side matrix C.
C
C   task =  4: Compute R, X and Y.
C
C   task =  5: Compute R, a low rank Cholesky factor of T'*T, i.e.,
C              R is lower triangular, so that R*R' = (T*E)'*(T*E)
C              for a permutation matrix E. Optionally, a matrix
C              Q is computed so that Q'*Q = I and T*E = Q*R'.
C
C   task =  6: Compute R, the Cholesky factor of T'*T, where T is a
C              banded block Toeplitz matrix.
C
C   task =  7: Compute R and the least-squares solution of T*X = B for a
C              banded block Toeplitz matrix T and given right hand side
C              matrix B.
C
C   task =  8: Compute R and the minimum norm solution of T'*Y = C for a
C              banded block Toeplitz matrix and given right hand side
C              matrix C.
C
C   task =  9: Compute R, X and Y for a banded block Toeplitz matrix.
C
C   task = 10: Compute the least-squares solution of T*X = B for a block
C              Toeplitz matrix T and a given right hand side matrix B.
C
C   task = 11: Compute the minimum norm solution of T'*Y = C for a given
C              right hand side matrix C.
C
C   task = 12: Compute X and Y.
C
C   task = 13: Compute the matrix-vector products X = T*B.
C
C Input parameters:
C   task  - integer option to determine the computation to perform as
C           described above.
C   TC    - real M*K-by-L matrix, containing the first block column of
C           the block Toeplitz matrix T with K-by-L blocks.
C   TR    - real K-by-N*L matrix, containing the first block row of the
C           block Toeplitz matrix T with K-by-L blocks. If the first
C           block of TR differs from the first block of TC a warning
C           will be displayed that the first block of TC will be used
C           as block diagonal element of the matrix T.
C   B     - real right hand side M*K-by-NCB (task < 13) / N*L-by-NCB
C           (task = 13) matrix.
C   C     - real right hand side N*L-by-NCC matrix.
C   tol   - real tolerance used to estimate the numerical rank of T.
C           (See Section METHOD of the SLICOT Library routine MB02JX,
C           argument TOL1.)
C   m     - integer containing the number of block rows of T.
C   n     - integer containing the number of block columns of T.
C   TCB   - real (ML+1)*K-by-L matrix, containing the leading nonzero
C           blocks of the first block column of the banded block
C           Toeplitz matrix T with K-by-L blocks.
C   TRB   - real K-by-(NU+1)*L matrix, containing the leading nonzero
C           blocks of the first block row of the banded block Toeplitz
C           matrix T with K-by-L blocks. If the first block of TRB
C           differs from the first block of TCB a warning will be
C           displayed that the first block of TCB will be used as block
C           diagonal element of the matrix T.
C
C Output parameters:
C
C   Q     - the M*K-by-RNK orthogonal factor satisfying T = Q*R' /
C           T*E = Q*R', where RNK is the numerical rank of T. Note that
C           if task<>5 then RNK = MIN(M*K,N*L).
C   R     - the N*L-by-RNK Cholesky factor of T'*T / (T*E)'*(T*E)
C           (task < 6),
C           the MIN(ML+NU+1,N)*L-by-MIN(M*K,N*L) matrix containing the
C           Cholesky factor of T'*T in banded storage format
C           (6 <= task <= 9).
C   X     - the N*L-by-NCB least-squares solution of T*X = B (task < 13),
C           the M*K-by-NCB matrix containing the matrix-vector products
C           T*B (task = 13).
C   Y     - the M*K-by-NCC minimum norm solution of T'*Y = C.
C   E     - MIN(M*K,N*L) vector recording the column pivoting performed.
C           If E(j) = k, then the j-th column of T*P was the k-th column
C           of T, where P is the permutation matrix.
C
C Comments:
C   1.    - If only least-squares/minimum norm solutions are desired
C           then task = 6..9 should be used for efficiency.
C   2.    - Note that for the computation of the least-squares/minimum
C           norm solution, a semi-normal equation approach is used which
C           does not yield numerically backward stable solutions.
C           Iterative refinement should be used if improved accuracy
C           is required.
C
C Contributor:
C   D. Kressner, TU Berlin, Aug. 2002.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Aug. 2002,
C   Apr. 2009, Dec. 2012.
C
C **********************************************************************
C
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C
C     .. Mex-file interface parameters ..
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C     .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mexEvalString, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsComplex, mxIsNumeric
C
C     .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         JOB
      INTEGER           INFO, K, L, LDB, LDC, LDQ, LDR, LDTC, LDTR,
     $                  LDWORK, LDX, LDY, M, ML, N, NU, RNK, S
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER, ALLOCATABLE ::          JPVT(:)
      DOUBLE PRECISION, ALLOCATABLE :: B(:,:), C(:,:), DWORK(:), Q(:,:),
     $                                 R(:,:), TC(:,:), TR(:,:), X(:,:),
     $                                 Y(:,:)
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           BADK, BADL, GETX, GETY
      INTEGER           I, IP, IPB, IPC, J, LEN, MKNL, NB, NCB, NCC,
     $                  NCQ, NCR, NCTC, NCTR, NCX, NCY, NLM, NRB, NRC,
     $                  NRQ, NRR, NRTC, NRTR, NRX, NRY, POSX, PT, TASK
      DOUBLE PRECISION  TEMP, TOL1, TOL2
C
C     .. External Functions ..
      INTEGER           ILAENV
      EXTERNAL          ILAENV
C
C     .. External subroutines ..
      EXTERNAL          DGEMM, DLACPY, DLASET, DTBTRS, DTRSM, MB02HD,
     $                  MB02ID, MB02JD, MB02JX, MB02KD
C
C     ..Intrinsic functions..
      INTRINSIC         DBLE, INT, MAX, MIN, MOD
C
C     Check for proper number of arguments.
C
      IF ( NRHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'FSTOEQ requires at least 3 input arguments' )
      ELSE IF ( NLHS.GT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'FSTOEQ requires at most 4 output arguments' )
      END IF
C
C     Check dimensions of input parameters and read/set scalar
C     parameters.
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
C
      IF ( TASK.LT.1 .OR. TASK.GT.13 ) THEN
         CALL mexErrMsgTxt
     $        ( 'The only admissible values of TASK are 1, 2, ..., 13.')
      END IF
C
      JOB = 'R'
      IF ( TASK.EQ.1 ) THEN
         IF ( NLHS.GT.1 )
     $      JOB = 'Q'
         NLM = 2
      ELSE IF ( TASK.EQ.2 ) THEN
         IF ( NLHS.GT.2 )
     $      JOB = 'Q'
         NLM = 3
      ELSE IF ( TASK.EQ.3 ) THEN
         IF ( NLHS.GT.2 )
     $      JOB = 'Q'
         NLM = 3
      ELSE IF ( TASK.EQ.4 ) THEN
         IF ( NLHS.GT.3 )
     $      JOB = 'Q'
         NLM = 4
      ELSE IF ( TASK.EQ.5 ) THEN
         IF ( NLHS.GT.2 )
     $      JOB = 'Q'
         NLM = 3
      ELSE IF ( TASK.EQ.6 ) THEN
         NLM = 1
      ELSE IF ( TASK.EQ.7 ) THEN
         NLM = 2
      ELSE IF ( TASK.EQ.8 ) THEN
         NLM = 2
      ELSE IF ( TASK.EQ.9 ) THEN
         NLM = 3
      ELSE IF ( TASK.EQ.10 ) THEN
         JOB = 'O'
         NLM = 1
      ELSE IF ( TASK.EQ.11 ) THEN
         JOB = 'U'
         NLM = 1
      ELSE IF ( TASK.EQ.12 ) THEN
         JOB = 'A'
         NLM = 2
      ELSE IF ( TASK.EQ.13 ) THEN
         NLM = 1
      END IF
C
      IF ( NLHS.GT.NLM ) THEN
         WRITE( TEXT, '('' FSTOEQ requires at most '',I4,
     $            '' output arguments'')' ) NLM
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( TASK.EQ.1 .AND. NRHS.GT.3 ) THEN
         CALL mexErrMsgTxt
     $           ( 'FSTOEQ requires at most 3 input arguments' )
      ELSE IF ( ( TASK.EQ.2 .OR. TASK.EQ.3 .OR. TASK.EQ.5 .OR.
     $            TASK.EQ.10 .OR. TASK.EQ.11 .OR. TASK.EQ.13 ) .AND.
     $          NRHS.GT.4 ) THEN
         CALL mexErrMsgTxt
     $           ( 'FSTOEQ requires at most 4 input arguments' )
      ELSE IF ( ( TASK.EQ.4 .OR. TASK.EQ.6 .OR. TASK.EQ.12 ) .AND.
     $          NRHS.GT.5 ) THEN
         CALL mexErrMsgTxt
     $           ( 'FSTOEQ requires at most 5 input arguments' )
      ELSE IF ( ( TASK.EQ.7 .OR. TASK.EQ.8 ) .AND. NRHS.GT.6 ) THEN
         CALL mexErrMsgTxt
     $           ( 'FSTOEQ requires at most 6 input arguments' )
      END IF
C
      IF ( ( TASK.EQ.2 .OR. TASK.EQ.3 .OR. TASK.EQ.5 .OR. TASK.EQ.10
     $       .OR. TASK.EQ.11 .OR. TASK.EQ.13 ) .AND.
     $          NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $           ( 'FSTOEQ requires at least 4 input arguments' )
      ELSE IF ( ( TASK.EQ.4 .OR. TASK.EQ.6 .OR. TASK.EQ.12 ) .AND.
     $          NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $           ( 'FSTOEQ requires at least 5 input arguments' )
      ELSE IF ( ( TASK.EQ.7 .OR. TASK.EQ.8 ) .AND. NRHS.LT.6 ) THEN
         CALL mexErrMsgTxt
     $           ( 'FSTOEQ requires at least 6 input arguments' )
      ELSE IF ( TASK.EQ.9 .AND. NRHS.LT.7 ) THEN
         CALL mexErrMsgTxt
     $           ( 'FSTOEQ requires at least 7 input arguments' )
      END IF
C
C     Check integer leading parameters if applicable.
C
      IF ( TASK.GE.6 .AND. TASK.LE.9 ) THEN
         IF ( mxGetM( PRHS(2) ).NE.1 .OR. mxGetN( PRHS(2) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'M must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(2) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'M must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), TEMP, 1 )
         M = TEMP
         IF ( mxGetM( PRHS(3) ).NE.1 .OR. mxGetN( PRHS(3) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'N must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(3) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'N must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TEMP, 1 )
         N = TEMP
         IP = 4
      ELSE
         IP = 2
      END IF
C
C     Check TC and TR.
C
      NRTC = mxGetM( PRHS(IP) )
      NCTC = mxGetN( PRHS(IP) )
      NRTR = mxGetM( PRHS(IP+1) )
      NCTR = mxGetN( PRHS(IP+1) )
      K = NRTR
      L = NCTC
      IF ( NRTC.LT.0 .OR. L.LT.0 ) THEN
         CALL mexErrMsgTxt( 'TC must be a matrix' )
      ELSE IF ( TASK.GE.6 .AND. TASK.LE.9 .AND. NRTC.EQ.0 ) THEN
         CALL mexErrMsgTxt( 'TC must have at least one row' )
      END IF
      IF ( NCTR.LT.0 .OR. K.LT.0 ) THEN
         CALL mexErrMsgTxt( 'TR must be a matrix' )
      ELSE IF ( TASK.GE.6 .AND. TASK.LE.9 .AND. NRTR.EQ.0 ) THEN
         CALL mexErrMsgTxt( 'TR must have at least one column' )
      END IF
      IF ( K.NE.0 ) THEN
         BADK = MOD( NRTC, K ).NE.0
      ELSE
         BADK = NRTC.NE.0
      END IF
      IF ( L.NE.0 ) THEN
         BADL = MOD( NCTR, L ).NE.0
      ELSE
         BADL = NCTR.NE.0
      END IF
      IF ( BADK .OR. BADL ) THEN
         CALL mexErrMsgTxt( 'Dimensions of TC and TR do not match' )
      END IF
C
      IF ( TASK.GE.6 .AND. TASK.LE.9 ) THEN
         ML = NRTC / K - 1
         IF ( ML.GE.M ) THEN
            CALL mexErrMsgTxt(
     $'The number of block rows of TC must not be greater than M' )
         ELSE IF ( ( ML + 1 )*K.LT.L ) THEN
            CALL mexErrMsgTxt(
     $'The number of columns in TC must not be greater than the number o
     $f nonzero rows' )
         END IF
         NU = NCTR / L - 1
         IF ( NU.GE.N ) THEN
            CALL mexErrMsgTxt(
     $'The number of block columns of TR must not be greater than N' )
         END IF
         IF ( ( M*K.LE.N*L.AND.( ( ML.LT.M - INT( ( M*K - 1 )/L ) - 1 )
     $        .OR.( ML.LT.M - INT( M*K/L ).AND.MOD( M*K, L ).LT.K ) ) )
     $        .OR.( M*K.GE.N*L .AND. ML*K.LT.N*( L - K ) )
     $        .OR.( M + NU )*L.LT.MIN( M*K, N*L ) ) THEN
            CALL mexErrMsgTxt(
     $'The first MIN(M*K,N*L) columns of T must have full rank.' )
         END IF
      ELSE
         IF ( K.NE.0 ) THEN
            M = NRTC / K
         ELSE
            M = NRTC
         END IF
         IF ( L.NE.0 ) THEN
            N = NCTR / L
         ELSE
            N = NCTR
         END IF
      END IF
C
C     Solving least-squares problems using seminormal equations is only
C     possible if T has full column rank, which induces M*K >= N*L.
C
      IF ( ( ( TASK.GE.2 .AND. TASK.LE.4 ) .OR.
     $       ( TASK.GE.7 .AND. TASK.LE.9 ) .OR.
     $       ( TASK.GE.10 .AND. TASK.LE.12 ) ) .AND. M*K.LT.N*L )
     $   CALL mexErrMsgTxt( 'T must have full column rank.' )
      IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'TC must be a real matrix' )
      ELSE IF ( mxIsNumeric( PRHS(IP+1) ).EQ.0 .OR.
     $          mxIsComplex( PRHS(IP+1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'TR must be a real matrix' )
      END IF
C
C     Check B.
C
      GETX = ( TASK.EQ.2 .OR. TASK.EQ.4 .OR. TASK.EQ.7 .OR. TASK.EQ.9
     $         .OR. TASK.EQ.10 .OR. TASK.EQ.12 .OR. TASK.EQ.13 )
      IF ( GETX ) THEN
         IF ( TASK.EQ.7 .OR. TASK.EQ.9 ) THEN
            IPB = 6
         ELSE
            IPB = 4
         END IF
         NRB = mxGetM( PRHS(IPB) )
         NCB = mxGetN( PRHS(IPB) )
         IF ( NRB.LT.0 .OR. NCB.LT.0 )
     $      CALL mexErrMsgTxt( 'B must be a matrix' )
         IF ( ( TASK.EQ.13.AND.NRB.NE.N*L ).OR.
     $        ( TASK.NE.13.AND.NRB.NE.M*K ) )
     $      CALL mexErrMsgTxt( 'Dimensions of T and B do not match.' )
         IF ( mxIsNumeric( PRHS(IPB) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IPB) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'B must be a real matrix' )
      ELSE
         NRB = 0
         NCB = 0
      END IF
C
C     Check C.
C
      GETY = ( TASK.EQ.3 .OR. TASK.EQ.4 .OR. TASK.EQ.8 .OR. TASK.EQ.9
     $         .OR. TASK.EQ.11 .OR. TASK.EQ.12 )
      IF ( GETY ) THEN
         IPC = NRHS
         NRC = mxGetM( PRHS(IPC) )
         NCC = mxGetN( PRHS(IPC) )
         IF ( NRC.LT.0 .OR. NCC.LT.0 )
     $      CALL mexErrMsgTxt( 'C must be a matrix' )
         IF ( NRC.NE.N*L )
     $      CALL mexErrMsgTxt( 'Dimensions of T and C do not match.' )
         IF ( mxIsNumeric( PRHS(IPC) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IPC) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'C must be a real matrix' )
      ELSE
         NRC = 0
         NCC = 0
      END IF
C
C     Check TOL.
C
      IF ( TASK.EQ.5 ) THEN
         IF ( mxGetM( PRHS(4) ).NE.1 .OR. mxGetN( PRHS(4) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(4) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real scalar' )
         END IF
      END IF
C
C     Determine the lenghts of working arrays.
C
      LDTC = MAX( 1, NRTC )
      LDTR = MAX( 1, NRTR )
      LDB = MAX( 1, NRB )
      LDC = MAX( 1, NRC )
      NCR = 0
      NRR = 0
      LDR = 1
      NCQ = 0
      NRQ = 0
      LDQ = 1
      NCY = 0
      NRY = 0
      LDY = 1
      MKNL = MIN( M*K, N*L )
      IF ( TASK.LE.4 ) THEN
         NRR = N*L
         LDR = MAX( 1, NRR )
         NCR = MKNL
      END IF
      IF ( TASK.EQ.5 ) THEN
         NRR = N*L
         LDR = MAX( 1, NRR )
         NCR = N*L
      END IF
      IF ( TASK.GE.6 .AND. TASK.LE.9 ) THEN
         NRR = MIN( ML + NU + 1, N )*L
         LDR = MAX( 1, NRR )
         NCR = MKNL
      END IF
      IF ( TASK.LE.4 .AND. NLHS.GE.NRHS-1 ) THEN
         NRQ = M*K
         LDQ = MAX( 1, NRQ )
         NCQ = MKNL
      END IF
      IF ( TASK.EQ.5 .AND. NLHS.GE.NRHS-1 ) THEN
         NRQ = M*K
         LDQ = MAX( 1, NRQ )
         NCQ = N*L
      END IF
      IF ( GETX ) THEN
C
C        Create temporary array X to hold intermediate results.
C
         IF ( TASK.LE.12 ) THEN
            NRX = N*L
            LDX = MAX( 1, NRX )
            NCX = NCB
         ELSE IF ( TASK.EQ.13 ) THEN
            NRX = M*K
            LDX = MAX( 1, NRX )
            NCX = NCB
         END IF
      END IF
      IF ( GETY ) THEN
C
C        Create temporary array Y to hold intermediate results.
C
         IF ( TASK.LE.12 ) THEN
            NRY = M*K
            LDY = MAX( 1, NRY )
            NCY = NCC
         END IF
      END IF
C
C     Workspace requirements to enable blocked algorithms.
C
      IF ( TASK.LE.4 ) THEN
         NB = MIN( ILAENV( 1, 'DGELQF', ' ', MKNL + L, L, -1, -1 ), L )
         LDWORK = MAX( 1, ( MAX( M*K, N*L ) + L )*NB
     $                    + ( M*K + N*L )*( L + 2*K )
     $                    + 6*L + M*K + N*L )
      ELSE IF ( TASK.EQ.5 ) THEN
         LDWORK = MAX( 3, ( M*K + ( N - 1 )*L )*( L + 2*K ) + 9*L
     $                    + MAX( M*K, ( N - 1 )*L ) )
      ELSE IF ( TASK.GE.6 .AND. TASK.LE.9 ) THEN
         NB = MIN( ILAENV( 1, 'DGELQF', ' ', ( ML + NU + 1 )*L,
     $                                       L, -1, -1 ), L )
         LDWORK = 1 + MAX( ( ML + NU + 1 )*L*L + ( 2*NU + L )*L*K,
     $                     ( 2*( L + K ) + MAX( NB , 1 ) )*
     $                     ( ML + NU + 1 )*L + 6*L )
      ELSE IF ( TASK.GE.10 .AND. TASK.LE.12 ) THEN
         NB = MIN( ILAENV( 1, 'DGELQF', ' ', MKNL + L, L, -1, -1 ), L )
         LDWORK = ( MAX(M*K,N*L) + L )*NB + 2*N*L*( L + K )
     $            + ( 6 + N )*L + M*K*( L + 1 ) + N*L*MAX( NCB, NCC )
     $            + 1
      ELSE IF ( TASK.EQ.13 ) THEN
         LDWORK = MAX(1, 2*( K*L + K*NCB + L*NCB + 1 )*( M + N ) )
      END IF
C
C     Allocate variable dimension local arrays.
C     !Fortran 90/95
C
      ALLOCATE ( DWORK(LDWORK), TC(LDTC,NCTC), TR(LDTR,NCTR) )
      IF ( TASK.LE.9 )
     $   ALLOCATE ( R(LDR,NCR) )
      IF ( TASK.LE.5 )
     $   ALLOCATE ( Q(LDQ,NCQ) )
      IF ( GETX.OR.GETY )
     $   ALLOCATE ( B(LDB,NCB), C(LDC,NCC) )
      IF ( GETX )
     $   ALLOCATE ( X(LDX,NCX) )
      IF ( ( GETY .AND. TASK.LE.9 ) .OR.
     $     ( TASK.GE.10 .AND. TASK.LE.12 ) )
     $   ALLOCATE ( Y(LDY,NCY) )
      IF ( TASK.EQ.5 )
     $   ALLOCATE ( JPVT(MKNL) )
      CALL DLASET( 'All', NRR, NCR, ZERO, ZERO, R, LDR )
C
C     Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TC, NRTC*NCTC )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP+1) ), TR, NRTR*NCTR )
      IF ( GETX )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(IPB) ), B, NRB*NCB )
      IF ( GETY )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(IPC) ), C, NRC*NCC )
      IF ( TASK.EQ.5 )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), TOL1, 1 )
C
C     Check whether the first blocks of TC and TR are equal.
C
      IF ( NRTC.GE.1 .AND. NCTR.GE.1 ) THEN
         DO 20  J = 1, L
            DO 10  I = 1, K
               IF ( TC(I,J).NE.TR(I,J) )
     $            GO TO 30
   10       CONTINUE
   20    CONTINUE
         GO TO 40
   30    CONTINUE
         I = mexEvalString('warning(
     $   ''Block column wins block diagonal conflict in TC/TR.'');');
   40    CONTINUE
      END IF
C
C     Do the actual computations.
C
      IF ( TASK.LE.4 ) THEN
         INFO = 0
         IF ( M*K.GT.0 .AND. N*L.GT.0 ) THEN
            S = ( MKNL + L - 1 ) / L
            CALL MB02JD( JOB, K, L, M, N, 0, S, TC, LDTC,
     $                   TR(1,MIN( NCTR, L+1 )), LDTR, Q, LDQ, R, LDR,
     $                   DWORK, LDWORK, INFO )
         END IF
         IF ( INFO.EQ.0 ) THEN
            IF ( GETX ) THEN
               CALL MB02KD( 'Column', 'Transpose', K, L, M, N, NCB, ONE,
     $                      ZERO, TC, LDTC, TR(1,MIN( NCTR, L+1 )),
     $                      LDTR, B, LDB, X, LDX, DWORK, LDWORK, INFO )
               CALL DTRSM( 'Left', 'Lower', 'No Transpose',
     $                     'NoUnitDiag', N*L, NCB, ONE, R, LDR, X, LDX )
               CALL DTRSM( 'Left', 'Lower', 'Transpose', 'NoUnitDiag',
     $                     N*L, NCB, ONE, R, LDR, X, LDX )
            END IF
            IF ( GETY ) THEN
               CALL DTRSM( 'Left', 'Lower', 'No Transpose',
     $                     'NoUnitDiag', N*L, NCC, ONE, R, LDR, C, LDC )
               CALL DTRSM( 'Left', 'Lower', 'Transpose', 'NoUnitDiag',
     $                      N*L, NCC, ONE, R, LDR, C, LDC )
               CALL MB02KD( 'Column', 'No Transpose', K, L, M, N, NCC,
     $                      ONE, ZERO, TC, LDTC, TR(1,MIN( NCTR, L+1 )),
     $                      LDTR, C, LDC, Y, LDY, DWORK, LDWORK, INFO )
            END IF
         END IF
      ELSE IF ( TASK.EQ.5 ) THEN
         TOL2 = -1.0D0
         CALL MB02JX( JOB, K, L, M, N, TC, LDTC, TR(1,MIN( NCTR, L+1 )),
     $                LDTR, RNK, Q, LDQ, R, LDR, JPVT, TOL1, TOL2,
     $                DWORK, LDWORK, INFO )
      ELSE IF ( TASK.GE.6 .AND. TASK.LE.9 ) THEN
         S = ( MKNL + L - 1 ) / L
         CALL MB02HD( 'No Structure', K, L, M, ML, N, NU, 0, S, TC,
     $                LDTC, TR(1,MIN( NCTR, L+1 )), LDTR, R, LDR, DWORK,
     $                LDWORK, INFO )
         IF ( INFO.EQ.0 ) THEN
            IF ( GETX.OR.GETY ) THEN
C
C              Put a flipped copy of TC in DWORK.
C
               PT = ML*K + 1
               DO 50  I = 1, ( ML + 1 )*K*L, K*L
                  CALL DLACPY( 'All', K, L, TC(PT,1), LDTC, DWORK(I),
     $                         K )
                  PT = PT - K
   50          CONTINUE
            END IF
            IF ( GETX .AND. NCB.GT.0 ) THEN
C
C              Overwrite X with T'*B.
C
               PT = ML*K*L + 1
               LEN = NU*L
               CALL DLASET( 'All', N*L, NCB, ZERO, ZERO, X, LDX )
               DO 60  I = 1, M
                  POSX = MAX(1,(I-ML-1)*L+1)
                  IF ( POSX.LE.N*L ) THEN
                     CALL DGEMM( 'Transpose', 'No Transpose',
     $                           MIN( MIN(I,ML+1)*L, N*L-POSX+1 ), NCB,
     $                           K, ONE, DWORK(PT), K, B((I-1)*K + 1,1),
     $                           LDB, ONE, X(POSX,1), LDX )
                  END IF
                  IF ( LEN.GT.0 ) THEN
                     CALL DGEMM( 'Transpose', 'No Transpose', LEN, NCB,
     $                           K, ONE, TR(1,L+1), LDTR,
     $                           B((I-1)*K + 1,1), LDB, ONE, X(I*L+1,1),
     $                           LDX )
                  END IF
                  IF ( PT.GT.1 )
     $               PT = PT - K*L
                  IF ( I.GE.N-NU )
     $               LEN = LEN - L
   60          CONTINUE
               CALL DTBTRS( 'Lower', 'No Transpose', 'NoUnitDiag', N*L,
     $                      MIN( ML + NU + 1, N )*L - 1, NCB, R, LDR, X,
     $                      LDX, INFO )
               CALL DTBTRS( 'Lower', 'Transpose', 'NoUnitDiag', N*L,
     $                      MIN( ML + NU + 1, N )*L - 1, NCB, R, LDR, X,
     $                      LDX, INFO )
            END IF
C
            IF ( GETY .AND. NCC.GT.0 ) THEN
               CALL DTBTRS( 'Lower', 'No Transpose', 'NoUnitDiag', N*L,
     $                      MIN( ML + NU + 1, N )*L - 1, NCC, R, LDR, C,
     $                      LDC, INFO )
               CALL DTBTRS( 'Lower', 'Transpose', 'NoUnitDiag', N*L,
     $                      MIN( ML + NU + 1, N )*L - 1, NCC, R, LDR, C,
     $                      LDC, INFO )
C
C              Overwrite Y with T*C.
C
               PT = ML*K*L + 1
               LEN = NU*L
               CALL DLASET( 'All', M*K, NCC, ZERO, ZERO, Y, LDY )
               DO 70  I = 1, M
                  POSX = MAX(1,(I-ML-1)*L+1)
                  IF ( POSX.LE.N*L ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', K, NCC,
     $                           MIN( MIN(I,ML+1)*L, N*L-POSX+1 ), ONE,
     $                           DWORK(PT), K, C(POSX,1), LDC,
     $                           ONE, Y((I-1)*K + 1,1), LDY )
                  END IF
                  IF ( LEN.GT.0 ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', K, NCC,
     $                           LEN, ONE, TR(1,L+1), LDTR, C(I*L+1,1),
     $                           LDC, ONE, Y((I-1)*K + 1,1), LDY )
                  END IF
                  IF ( PT.GT.1 )
     $               PT = PT - K*L
                  IF ( I.GE.N-NU )
     $               LEN = LEN - L
   70          CONTINUE
            END IF
         END IF
      ELSE IF ( TASK.GE.10 .AND. TASK.LE.12 ) THEN
         IF ( GETY )
     $      CALL DLACPY( 'All', NRC, NCC, C, LDC, Y, LDY )
         CALL MB02ID( JOB, K, L, M, N, NCB, NCC, TC, LDTC,
     $                TR(1,MIN( NCTR, L+1 )), LDTR, B, LDB, Y, LDY,
     $                DWORK, LDWORK, INFO )
         IF ( INFO.EQ.0 ) THEN
            IF ( GETX )
     $         CALL DLACPY( 'All', NRX, NCX, B, LDB, X, LDX )
         END IF
      ELSE IF ( TASK.EQ.13 ) THEN
         CALL MB02KD( 'Column', 'No Transpose', K, L, M, N, NCB, ONE,
     $                ZERO, TC, LDTC, TR(1,MIN( NCTR, L+1 )), LDTR, B,
     $                LDB, X, LDX, DWORK, LDWORK, INFO )
      END IF
C
C     Copy output to MATLAB workspace.
C
      IP = 1
      IF ( INFO.EQ.0 ) THEN
         IF ( TASK.LE.4 .AND. NLHS.GE.NRHS-1 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( NRQ, NCQ, 0 )
            CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(IP) ), NRQ*NCQ )
            IP = IP + 1
         END IF
         IF ( TASK.EQ.5 .AND. NLHS.GE.NRHS-1 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( NRQ, RNK, 0 )
            CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(IP) ), NRQ*RNK )
            IP = IP + 1
         END IF
         IF ( TASK.LE.9 .AND. TASK.NE.5 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( NRR, NCR, 0 )
            IF ( MIN( NRR, NCR-1 ).GE.1 .AND. TASK.LE.4 ) THEN
               CALL DLASET( 'Upper part', NRR-1, NCR-1, ZERO, ZERO,
     $                      R(1,2), LDR )
            END IF
            CALL mxCopyReal8ToPtr( R, mxGetPr( PLHS(IP) ), NRR*NCR )
            IP = IP + 1
         END IF
         IF ( TASK.EQ.5 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( NRR, RNK, 0 )
            CALL mxCopyReal8ToPtr( R, mxGetPr( PLHS(IP) ), NRR*RNK )
            IP = IP + 1
         END IF
         IF ( GETX .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( NRX, NCX, 0 )
            CALL mxCopyReal8ToPtr( X, mxGetPr( PLHS(IP) ), NRX*NCX )
            IP = IP + 1
         END IF
         IF ( GETY .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( NRY, NCY, 0 )
            CALL mxCopyReal8ToPtr( Y, mxGetPr( PLHS(IP) ), NRY*NCY )
            IP = IP + 1
         END IF
         IF ( TASK.EQ.5 .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( 1, MKNL, 0 )
            DO 80 I = 1, MKNL
               DWORK(I) = DBLE( JPVT(I) )
   80       CONTINUE
            CALL mxCopyReal8ToPtr( DWORK, mxGetPr( PLHS(IP) ), MKNL )
            IP = IP + 1
         END IF
      END IF
C
C     Deallocate local arrays.
C     !Fortran 90/95
C
      DEALLOCATE( DWORK, TC, TR )
      IF ( TASK.LE.9 )
     $   DEALLOCATE( R )
      IF ( TASK.LE.5 )
     $   DEALLOCATE( Q )
      IF ( GETX.OR.GETY )
     $   DEALLOCATE( B, C )
      IF ( GETX )
     $   DEALLOCATE( X )
      IF ( ( GETY .AND. TASK.LE.9 ) .OR.
     $     ( TASK.GE.10 .AND. TASK.LE.12 ) )
     $   DEALLOCATE( Y )
      IF ( TASK.EQ.5 )
     $   DEALLOCATE( JPVT )
C
C     Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( TASK.LE.4 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB02JD'')' )
     $             INFO
         ELSE IF ( TASK.LE.5 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB02JX'')' )
     $             INFO
         ELSE IF ( TASK.LE.9 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB02HD'')' )
     $             INFO
         ELSE IF ( TASK.LE.12 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB02ID'')' )
     $             INFO
         ELSE IF ( TASK.EQ.13 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB02KD'')' )
     $             INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of FSTOEQ ***
      END
