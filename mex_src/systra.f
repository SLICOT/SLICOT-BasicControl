C SYSTRA.F - Gateway function for computing system transformations
C            with scaling, block diagonal decomposition, or Schur form
C            reduction of the state matrix, using the SLICOT routines
C            TB01ID, TB01KD, TB01LD, and TB01WD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [A,B,C(,num,ev,U)] = systra(task,A,B,C(,job,par))
C   [A,B,C(,num,ev,U)] = systra(task,A,B,C(,flag,par))
C
C   task = 1 :        [A,B,C,U] = systra(1,A,B,C,job,par)
C   task = 2 :     [A,B,C,ev,U] = systra(2,A,B,C)
C   task = 3 : [A,B,C,num,ev,U] = systra(3,A,B,C,flag,par)
C   task = 4 : [A,B,C,num,ev,U] = systra(4,A,B,C,flag,par)
C
C Purpose:
C   To apply a specified similarity transformation to a given system
C              [ A  B ]
C   matrix S = [      ]. The possible transformations are:
C              [ C  0 ]
C
C task = 1                                            [ D  0 ]
C   Do balance on the matrix S with a diagonal matrix [      ];
C                                                     [ 0  I ]
C task = 2
C   Reduce A to Schur form using an orthogonal matrix U, and transform
C   B and C accordingly,
C                                [* . . . * |* . *]--
C                                [  .     . |.   .] |
C  [ U'| 0 ][ A | B ][ U | 0 ]   [    .   . |.   .] N
C  [---|---][---|---][---|---] = [      . . |.   .] |  ;      (1)
C  [ 0 | I ][ C | 0 ][ 0 | I ]   [        * |* . *]--
C                                [----------|-----]
C                                [* . . . * |     ]-|-
C                                [.       . |  0  ] P
C                                [* . . . * |     ]-|-
C                                 |-- N --|  |-M-|
C task = 3
C   Reduce A to real Schur form with eigenvalues in required order, by
C   using an orthogonal matrix U, and transform B and C accordingly,
C
C                               |--[* . * * . * |* . *]--
C                              num [  . . .   . |* . *] |
C                               |--[    * * . * |.   .] |
C  [ U'| 0 ][ A | B ][ U | 0 ]     [      * . * |.   .] N
C  [---|---][---|---][---|---] =   [        . . |.   .] |  ;   (2)
C  [ 0 | I ][ C | 0 ][ 0 | I ]     [          * |* . *]--
C                                  [------------|-----]
C                                  [* . . . . * |     ]-|-
C                                  [.         . |  0  ] P
C                                  [* . . . . * |     ]-|-
C                                   |--- N ---|  |-M-|
C task = 4
C  Reduce A to a block diagonal form with eigenvalues in required order,
C  using a nonsingular real matrix U, and transform B and C accordingly,
C
C                                 |--[* . *       |* . *]--
C                                num [  . .       |* . *] |
C     -1                          |--[    *       |.   .] |
C  [ U  | 0 ][ A | B ][ U | 0 ]      [      * . * |.   .] N
C  [----|---][---|---][---|---] =    [        . . |.   .] | .   (3)
C  [ 0  | I ][ C | 0 ][ 0 | I ]      [          * |* . *]--
C                                    [------------|-----]
C                                    [* . . . . * |     ]-|-
C                                    [.         . |  0  ] P
C                                    [* . . . . * |     ]-|-
C                                     |--- N ---|  |-M-|
C
C Input parameters:
C   task  - integer option to indicate which transformation is needed:
C           = 1 : do balancing;
C           = 2 : compute the transformation (1);
C           = 3 : compute the transformation (2);
C           = 4 : compute the transformation (3).
C   A     - real N-by-N state matrix.
C   B     - real N-by-M input matrix.
C   C     - real P-by-N output matrix.
C   job   - (for task = 1 only: optional) integer option parameter:
C           job = 1  : balancing involves the matrix A only.
C           job = 2  : balancing involves the matrices A and B.
C           job = 3  : balancing involves the matrices A and C.
C           otherwise, balancing involves all matrices A, B, and C.
C           default:  job = 0;
C   flag  - (optional) integer vector containing options.
C           task = 1, 2 : flag is not used.
C           task = 3, 4 : flag is a vector of length 2
C              flag(1) = 0 : the system is continuous-time (default),
C                            otherwise, the system is discrete-time.
C              flag(2) = 0 : reorder the stable eigenvalues of A on the
C                            top left diagonal block (default),
C                            otherwise, reorder the unstable eigenvalues
C                            on the top left diagonal block.
C   par   - (optional) real parameter specifying the following values:
C           task = 1    : the maximum allowed reduction in the 1-norm of
C                         S if zero rows or columns are encountered;
C                         default:  par = 10;
C           task = 2    : not used.
C           task = 3, 4 : the stability or instability boundary for the
C                         eigenvalues of A of interest. For the
C                         discrete-time case, par >= 0 represents the
C                         boundary value for the moduli of eigenvalues.
C             default:    -sqrt(epsilon_machine) for continuous-time;
C                      1.0-sqrt(epsilon_machine) for discrete-time.
C
C Output parameters:
C   A     - transformed state matrix.
C   B     - transformed input matrix.
C   C     - transformed output matrix.
C   num   - (optional) integer counting for the number of the required
C           eigenvalues, or, equivalently, the size of the diagonal
C           block containing the required eigenvalues. If task = 1, 2,
C           num is not available.
C   ev    - (optional) complex vector of length N, containing the
C           eigenvalues of A, possibly in the required order.
C           When task = 1, ev is not available.
C   U     - (optional) nonsingular matrix which performs the
C           transformations.
C           task = 1    : U is N-by-1 storing the diagonal elements
C                         of D.
C           task = 2, 3 : U is orthogonal.
C           task = 4    : U is nonsingular.
C
C Contributor:
C   H. Xu, TU Chemnitz, FR Germany, Dec. 1998
C
C Revisions:
C   V. Sima, Katholieke Univ. Leuven, Belgium, May 1999, Apr. 2009,
C   Dec. 2012.
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
      CHARACTER         DICO, JOB, STDOM
      INTEGER           INFO, LDA, LDB, LDC, LDU, LDWORK, M, N, P
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), DWORK(:),
     $                                 U(:,:), V(:), WR(:), WI(:)
      COMPLEX*16,       ALLOCATABLE :: EV(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           FLAG(2), I, IJOB, IP, ISIZE, NUM, TASK
      DOUBLE PRECISION  FLAGR(2), PAR, RNUM, TEMP
C
C .. External functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C
C .. External subroutines ..
      EXTERNAL          TB01ID, TB01KD, TB01LD, TB01WD
C
C ..Intrinsic functions..
      INTRINSIC         DCMPLX, MAX, SQRT
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSTRA requires at least 4 input arguments' )
      ELSE IF ( NLHS.GT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSTRA requires at most 6 output arguments' )
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
     $        ( 'TASK has 1 ... 4 the only admissible values' )
      END IF
C
C   A(NxN), B(NxM), C(PxN)
C
      N = mxGetM( PRHS(2) )
      M = mxGetN( PRHS(3) )
      P = mxGetM( PRHS(4) )
C
      IF ( mxGetN( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'B must have the same row dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
C
      IF ( mxGetN( PRHS(4) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'C must have the same column dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
C
      FLAG(1) = 0
      IF ( TASK.EQ.1 ) THEN
         IJOB = 0
         PAR  = ZERO
      ELSE IF ( TASK.GE.3 ) THEN
         FLAG(2) = 0
         PAR = -SQRT( DLAMCH( 'Epsilon' ) )
      END IF
      IF ( NRHS.GE.5 ) THEN
         ISIZE = mxGetM( PRHS(5) ) * mxGetN( PRHS(5) )
         IF ( TASK.EQ.1 ) THEN
C
C   job
C
            IF ( ISIZE.GT.1 )
     $         CALL mexErrMsgTxt( 'JOB must be a scalar' )
            IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(5) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), TEMP, ISIZE )
            IF ( ISIZE.GT.0 ) IJOB = TEMP
         ELSE IF ( TASK.GE.3 ) THEN
C
C   flag
C
            IF ( ISIZE.GT.2 )
     $         CALL mexErrMsgTxt
     $              ( 'FLAG must be a vector with at most 2 elements' )
            IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(5) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'FLAG must be an integer vector' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), FLAGR, ISIZE )
            IF ( ISIZE.GT.0 ) FLAG(1) = FLAGR(1)
            IF ( ISIZE.GT.1 ) FLAG(2) = FLAGR(2)
            IF ( FLAG(1).NE.0 )
     $         PAR = ONE + PAR
         END IF
C
         IF ( NRHS.GE.6 ) THEN
            IF ( TASK.NE.2 ) THEN
               ISIZE = mxGetM( PRHS(6) ) * mxGetN( PRHS(6) )
C
C   par
C
               IF ( ISIZE.GT.1 )
     $            CALL mexErrMsgTxt( 'PAR must be a scalar' )
               IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $              mxIsComplex( PRHS(6) ).EQ.1 ) THEN
                  CALL mexErrMsgTxt( 'PAR must be a real scalar' )
               END IF
               CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), TEMP, ISIZE )
               IF ( ISIZE.GT.0 ) PAR = TEMP
            END IF
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms
C in DGEES. LDWORK >= N*MAX( 3, M, P ) will allow to use BLAS 3.
C
      LDWORK = MAX( 1, 3*N )
C
      LDA = MAX( 1, N )
      LDB = LDA
      LDC = MAX( 1, P )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE( A( LDA, N ), B( LDB, M ), C( LDC, N ) )
C
      IF ( TASK.EQ.1 ) THEN
         ALLOCATE( V( N ) )
      ELSE
         LDU = LDA
         ALLOCATE( U( LDU, N ), WR( N ), WI( N ), DWORK( LDWORK ) )
         ALLOCATE( EV( N ) )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), C, P*N )
C
C Do the actual computations.
C
      IF ( TASK.EQ.1 ) THEN
         IF ( IJOB.EQ.1 ) THEN
            JOB = 'N'
         ELSE IF ( IJOB.EQ.2 ) THEN
            JOB = 'B'
         ELSE IF ( IJOB.EQ.3 ) THEN
            JOB = 'C'
         ELSE
            JOB = 'A'
         END IF
         CALL TB01ID( JOB, N, M, P, PAR, A, LDA, B, LDB, C, LDC, V,
     $                INFO )
      ELSE IF ( TASK.EQ.2 ) THEN
         CALL TB01WD( N, M, P, A, LDA, B, LDB, C, LDC, U, LDU, WR, WI,
     $                DWORK, LDWORK, INFO )
      ELSE
         IF ( FLAG(1).EQ.0 ) THEN
            DICO = 'C'
         ELSE
            DICO = 'D'
         END IF
         IF ( FLAG(2).EQ.0 ) THEN
            STDOM = 'S'
         ELSE
            STDOM = 'U'
         END IF
C
         IF ( TASK.EQ.3 ) THEN
            CALL TB01LD( DICO, STDOM, 'G', N, M, P, PAR, A, LDA, B, LDB,
     $                   C, LDC, NUM, U, LDU, WR, WI, DWORK, LDWORK,
     $                   INFO )
         ELSE
            CALL TB01KD( DICO, STDOM, 'G', N, M, P, PAR, A, LDA, B, LDB,
     $                   C, LDC, NUM, U, LDU, WR, WI, DWORK, LDWORK,
     $                   INFO )
         END IF
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         IF ( NLHS.GE.1 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(1) ), N*N )
         END IF
         IF ( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( N, M, 0 )
            CALL mxCopyReal8ToPtr( B, mxGetPr( PLHS(2) ), N*M )
         END IF
         IF ( NLHS.GE.3 ) THEN
            PLHS(3) = mxCreateDoubleMatrix( P, N, 0 )
            CALL mxCopyReal8ToPtr( C, mxGetPr( PLHS(3) ), P*N )
         END IF
         IF ( NLHS.GE.4 ) THEN
            IP = 4
            IF ( TASK.EQ.1 ) THEN
               PLHS(4) = mxCreateDoubleMatrix( N, 1, 0 )
               CALL mxCopyReal8ToPtr( V, mxGetPr( PLHS(4) ), N )
            ELSE IF ( TASK.GE.3 ) THEN
               PLHS(4) = mxCreateDoubleMatrix( 1, 1, 0 )
               RNUM = NUM
               CALL mxCopyReal8ToPtr( RNUM, mxGetPr( PLHS(4) ), 1 )
               IP = IP + 1
            END IF
            IF ( TASK.GE.2 ) THEN
               IF ( NLHS.GE.IP ) THEN
                  DO 10 I = 1, N
                     EV(I) = DCMPLX( WR(I), WI(I) )
   10             CONTINUE
                  PLHS(IP) = mxCreateDoubleMatrix( N, 1, 1 )
                  CALL mxCopyComplex16ToPtr( EV, mxGetPr( PLHS(IP) ),
     $                                       mxGetPi( PLHS(IP) ), N )
                  IP = IP + 1
                  IF ( NLHS.GE.IP ) THEN
                     PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
                     CALL mxCopyReal8ToPtr( U, mxGetPr( PLHS(IP) ),
     $                                      N*N )
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      IF ( TASK.EQ.1 ) THEN
         DEALLOCATE( A, B, C, V )
      ELSE
         DEALLOCATE( A, B, C, U, WR, WI, EV, DWORK )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( TASK.EQ.1 ) THEN
            WRITE( TEXT,'('' INFO = '',I4,'' ON EXIT FROM TB01ID'')' )
     $         INFO
         ELSE IF ( TASK.EQ.2 ) THEN
            WRITE( TEXT,'('' INFO = '',I4,'' ON EXIT FROM TB01WD'')' )
     $         INFO
         ELSE IF ( TASK.EQ.3 ) THEN
            WRITE( TEXT,'('' INFO = '',I4,'' ON EXIT FROM TB01LD'')' )
     $         INFO
         ELSE
            WRITE( TEXT,'('' INFO = '',I4,'' ON EXIT FROM TB01KD'')' )
     $         INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C
C *** Last line of SYSTRA ***
      END
