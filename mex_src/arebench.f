C AREBENCH.F - Gateway function to generate the benchmark examples
C             for algebraic Riccati equations using SLICOT routines
C             BB01AD, for continuous-time case, and BB02AD, for
C             discrete-time case.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [X,A,G,Q(,B,C)(,S)(,parval)] = arebench(dico,nr1,nr2 ...
C                                           (,flag,param,opt,filnam))
C
C Purpose:
C   To generate benchmark examples for continuous-time algebraic
C   Riccati equations (CARE)
C
C     0  =  Q + A'X + XA - XGX,
C
C   where A,G,Q,X are real N-by-N matrices, Q and G are symmetric and
C   may be given in factored form
C                   -1 T                         T
C      (I)   G = B R  B  ,           (II)   Q = C W C ,
C   where C is P-by-N, W P-by-P, B N-by-M, and R M-by-M, where W
C   and R are symmetric,
C
C   or
C
C   to generate benchmark examples for the discrete-time algebraic
C   Riccati equations (DARE)
C            T            T               T    -1  T       T
C     0  =  A X A - X - (A X B + S) (R + B X B)  (B X A + S ) + Q,
C
C   where the matrices Q and R are symmetric and Q may be given in
C   factored form (II).
C   If R is nonsingular and S = 0, the DARE can be rewritten as
C                  T             -1
C     0  =  X  -  A X (I_n + G X)  A  -  Q,
C
C   where I_n is the N-by-N identity matrix and G is factored as in (I).
C
C Input parameters:
C   dico   - integer option to indicate continuous- or discrete-time:
C            = 1 : continuous-time example.
C            = 2 : discrete-time example.
C   nr1    - the group of examples:
C            = 1 : parameter-free problems of fixed size.
C            = 2 : parameter-dependent problems of fixed size.
C            = 3 : parameter-free problems of scalable size.
C            = 4 : parameter-dependent problems of scalable size.
C   nr2    - the number of the example in group nr1.
C            Let NEXi be the number of examples in group i. Currently,
C            NEX1 =  6, NEX2 = 9, NEX3 = 2, NEX4 = 4, for
C                                                     continuous-time,
C            NEX1 = 13, NEX2 = 5, NEX3 = 0, NEX4 = 1, for discrete-time.
C            1 <= nr1 <= 4;
C            1 <= nr2 <= NEXi , where i = nr1.
C   flag   - (optional) vector containing options:
C               flag(1) = 1  : G is returned.
C               flag(1) = 0  : G is returned in factored form, i.e.,
C                              B and R from (I) are returned.
C               flag(2) = 1  : Q is returned.
C               flag(2) = 0  : Q is returned in factored form, i.e.,
C                              C and W from (II) are returned.
C            If dico = 1, it has size 2; if dico = 2, it has size 3:
C               flag(3) = 1  : The coefficient matrix S of the DARE
C                              is returned in array S.
C               flag(3) = 0  : The coefficient matrix S of the DARE
C                              is not returned.
C   param  - (optional) real vector of parameter values:
C            = epsilon, if dico = 1, nr1 = 2;
C            = [q,r], if dico = 1, for 4.1 (i.e., nr1 = 4, nr2 = 1);
C            = [a,b,c,beta1,beta2,gamma1,gamma2], if dico = 1, for 4.2;
C            = [mu,delta,kappa], if dico = 1, nr1 = 4, nr2 = 3;
C            = epsilon, if dico = 2, nr1 = 2, nr2 = 2:4;
C            = [tau,D,K,r], if dico = 2, nr1 = 2, nr2 = 5;
C            = R, if dico = 2, nr1 = 2, nr2 = 1 or nr1 = 4, nr2 = 1;
C   opt    - (optional) integer vector for additional
C            options/parameters:
C            if dico = 1,
C            = 1, (when nr1 = 2, nr2 = 9): generates the CARE for
C                  optimal state feedback (default);
C            = 2, (when nr1 = 2, nr2 = 9): generates the
C                  Kalman filter CARE.
C            = the number of vehicles, if nr1 = 3, nr2 = 1.
C            = the order of the matrix A, for 3.2, 4.1 or 4.2.
C            = the dimension of the second-order system, i.e., the
C              order of the stiffness matrix for 4.3 or 4.4.
C            if dico = 2,
C            = the order of the output matrix A, for 4.1.
C   filnam - (optional) if nr1 = nr2 = 4 and opt differs from the
C            default dimension of the second-order system (211), a
C            character string containing the name of the data file
C            to be used.
C Output parameters:
C   X      - real symmetric solution of Riccati equation (an exact
C            solution is available for the continuous-time examples
C            1.1, 1.2, 2.1, 2.3-2.6, 3.2, and for the discrete-time
C            examples 1.1, 1.3, 1.4, 2.1, 2.3-2.5, 4.1;
C            otherwise, X is an empty array).
C   A      - real n-by-n coefficient matrix.
C   G      - real symmetric n-by-n matrix.
C            If flag(1) = 1, then a non-factored coefficient matrix G
C            of the ARE is returned.
C            If flag(1) = 0, then array G contains the
C            'control weighting matrix' R from (I).
C   Q      - real symmetric n-by-n matrix.
C            If flag(2) = 1, then a non-factored coefficient matrix Q
C            of the ARE is returned.
C            If flag(2) = 0, then array Q contains the
C            'output weighting matrix' W from (II).
C   B      - real n-by-m input matrix.
C            If flag(1) = 0, this array contains the coefficient
C            matrix B of the ARE.
C   C      - real p-by-n output matrix.
C            If flag(2) = 0, this array contains the coefficient
C            matrix C of the ARE.
C   S      - if (dico = 2, flag(3) = 1), then this array contains the
C            coefficient matrix S of the DARE.
C   parval - the values used for the problem parameters:
C            either the values given in param, or default values,
C            if param is not used.
C
C Contributor:
C   D. Sima, University of Bucharest, August 2001.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Nov. 2001, Oct. 2004,
C   Apr. 2009, Dec. 2012, May 2016, Jan. 2017.
C
C **********************************************************************
C
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Mex-file interface parameters ..
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxGetString, mxIsChar,
     $                  mxIsNumeric, mxIsComplex
C
C     .. Parameters ..
C     . # of examples available , # of examples with fixed size. .
      INTEGER           NEXC1, NEXC2, NEXC3, NEXC4, NEXD1, NEXD2, NEXD3,
     $                  NEXD4, NMX
      PARAMETER         ( NEXC1 =  6, NEXC2 = 9, NEXC3 = 2, NEXC4 = 4,
     $                    NEXD1 = 13, NEXD2 = 5, NEXD3 = 0, NEXD4 = 1,
     $                    NMX   = 13 )
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         DEF
      INTEGER           INFO, LDA, LDB, LDC, LDG, LDQ, LDS, LDWORK, LDX,
     $                  M, N, P
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:),   B(:,:), C(:,:),
     $                                 DWORK(:), G(:,:), Q(:,:), S(:,:),
     $                                 X(:,:)
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           BPAR(7), VEC(10)
      CHARACTER*255     CHPAR
      CHARACTER*120     TEXT
      INTEGER           DICO, FLAG(3), I, IE, IP, IPAR(4), ISIZE, L,
     $                  MDEFC(2,NMX), MDEFD(2,NMX), NDEFC(4,NMX),
     $                  NDEFD(4,NMX), NEX(4), NEXC(4), NEXD(4), NR(2),
     $                  PDEFC(2,NMX), PDEFD(2,NMX)
      DOUBLE PRECISION  DPAR(7), FLAGR(3), ITMP(7), TEMP
C
C .. External functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C
C .. External subroutines ..
      EXTERNAL          BB01AD, BB02AD, DLACPY
C
C ..Intrinsic functions..
      INTRINSIC         MAX
C
C     .. Data Statements ..
C     . default values for dimensions .
      DATA (NEXC(I), I = 1, 4) /NEXC1, NEXC2, NEXC3, NEXC4/
      DATA (NEXD(I), I = 1, 4) /NEXD1, NEXD2, NEXD3, NEXD4/
      DATA (NDEFC(1,I), I = 1, NEXC1) /2, 2, 4, 8, 9, 30/
      DATA (NDEFC(2,I), I = 1, NEXC2) /2, 2, 2, 2, 2, 3, 4, 4, 55/
      DATA (NDEFC(3,I), I = 1, NEXC3) /20, 64/
      DATA (NDEFC(4,I), I = 1, NEXC4) /21, 100, 30, 211/
      DATA (MDEFC(1,I), I = 1, NEXC1) /1, 1, 2, 2, 3, 3/
      DATA (MDEFC(2,I), I = 1, NEXC2) /1, 2, 1, 2, 1, 3, 1, 1, 2/
      DATA (PDEFC(1,I), I = 1, NEXC1) /2, 2, 4, 8, 9, 5/
      DATA (PDEFC(2,I), I = 1, NEXC2) /1, 1, 2, 2, 2, 3, 2, 1, 10/
      DATA (NDEFD(1,I), I = 1, NEXD1) /2, 2, 2, 3, 4, 4, 4, 5, 6, 9,
     $                                 11, 13, 26/
      DATA (NDEFD(2,I), I = 1, NEXD2) /2, 2, 2, 3, 4/
      DATA (NDEFD(4,I), I = 1, NEXD4) /100/
      DATA (MDEFD(1,I), I = 1, NEXD1) /1, 2, 1, 2, 2, 2, 4, 2, 2, 3,
     $                                 2, 2, 6/
      DATA (MDEFD(2,I), I = 1, NEXD2) /1, 2, 1, 3, 1/
      DATA (PDEFD(1,I), I = 1, NEXD1) /1, 2, 2, 3, 4, 4, 4, 5, 2, 2,
     $                                 4, 4, 12/
      DATA (PDEFD(2,I), I = 1, NEXD2) /2, 2, 2, 3, 1/
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'AREBENCH requires at least 3 input arguments' )
      ELSE IF ( NRHS.GT.6  ) THEN
         CALL mexErrMsgTxt
     $        ( 'AREBENCH requires at most 6 input arguments' )
      ELSE IF ( NLHS.GT.8 ) THEN
         CALL mexErrMsgTxt
     $        ( 'AREBENCH requires at most 8 output arguments' )
      END IF
C
      DEF = 'D'
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   dico
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'DICO must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'DICO must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      DICO = TEMP
      IF ( DICO.LT.1 .OR. DICO.GT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'DICO has 1 or 2 the only admissible values' )
      END IF
C
C   nr1, nr2
C
      IF ( mxGetM( PRHS(2) ).NE.1 .OR. mxGetN( PRHS(2) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'NR1 must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'NR1 must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), TEMP, 1 )
      NR(1) = TEMP
      IF ( NR(1).LT.1 .OR. NR(1).GT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'NR1 has 1,2,3, or 4 the only admissible values' )
      END IF
      IF ( DICO.EQ.2 .AND. NR(1).EQ.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'No available examples for this section' )
      END IF
C
      IF ( DICO.EQ.1 ) THEN
          NEX(1) = NEXC1
          NEX(2) = NEXC2
          NEX(3) = NEXC3
          NEX(4) = NEXC4
      ELSE
          NEX(1) = NEXD1
          NEX(2) = NEXD2
          NEX(3) = NEXD3
          NEX(4) = NEXD4
      END IF
C
      IF ( mxGetM( PRHS(3) ).NE.1 .OR. mxGetN( PRHS(3) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'NR2 must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'NR2 must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TEMP, 1 )
      NR(2) = TEMP
      IF ( NR(2).LT.1 .OR. NR(2).GT.NEX(NR(1)) ) THEN
         WRITE( TEXT, '( " NR2 has 1..", I4,
     $          " the only admisible values" )' ) NEX(NR(1))
         CALL mexErrMsgTxt ( TEXT )
      END IF
C
C   flag
C
      DO 10 I = 1, 7
         BPAR(I) = .TRUE.
   10 CONTINUE
      IP = 3
      IF ( NRHS.GT.IP ) THEN
        IP = IP + 1
        ISIZE = mxGetM( PRHS(IP) ) * mxGetN( PRHS(IP) )
        IF ( DICO.EQ.1 ) THEN
           IF ( ISIZE.NE.2 ) THEN
              CALL mexErrMsgTxt
     $              ( 'FLAG must be a vector with 2 elements' )
           END IF
        ELSE
           IF ( ISIZE.NE.3 ) THEN
              CALL mexErrMsgTxt
     $              ( 'FLAG must be a vector with 3 elements' )
           END IF
        END IF
        IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
           CALL mexErrMsgTxt( 'FLAG must be a real vector' )
        END IF
        CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), FLAGR, ISIZE )
        FLAG = FLAGR
C
        IF ( DICO.EQ.1 ) THEN
           IF ( FLAG(1).EQ.0 )
     $        BPAR(1) = .FALSE.
           IF ( FLAG(2).EQ.0 )
     $        BPAR(4) = .FALSE.
        ELSE
           IF ( FLAG(1).EQ.0 )
     $        BPAR(4) = .FALSE.
           IF ( FLAG(2).EQ.0 )
     $        BPAR(1) = .FALSE.
           IF ( FLAG(3).EQ.0 )
     $        BPAR(7) = .FALSE.
        END IF
      END IF
C
C  param
C
      IF ( NRHS.GT.IP .AND. ( NR(1).NE.4 .OR. NR(2).NE.4 ) ) THEN
        IP = IP + 1
        IF ( mxGetM( PRHS(IP) ).NE.0 .AND. mxGetN( PRHS(IP) ).NE.0 )
     $  THEN
           L = 0
           IF ( DICO.EQ.1 ) THEN
              IF ( NR(1).EQ.2 ) THEN
                 L = 1	
              ELSE IF ( NR(1).EQ.4 ) THEN
                 IF ( NR(2).EQ.1 ) THEN
                    L = 2
                 ELSE IF ( NR(2).EQ.2 ) THEN
                    L = 7
                 ELSE IF ( NR(2).EQ.3 ) THEN
                    L = 3
                 END IF
              END IF
           ELSE IF ( ( NR(1).EQ.2 .AND. NR(2).LE.4 ) .OR.
     $               ( NR(1).EQ.4 .AND. NR(2).EQ.1 ) ) THEN
              L = 1	
           ELSE IF ( NR(1).EQ.2 .AND. NR(2).EQ.5 ) THEN
              L = 4
           END IF
           IF ( L.NE.0 ) THEN
              IF ( mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) ).NE.L ) THEN
                 WRITE( TEXT,
     $                 '(''PARAM must be a vector of length '', I2)') L
                 CALL mexErrMsgTxt( TEXT )
              END IF
              IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $             mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
                 CALL mexErrMsgTxt( 'PARAM must be a real vector' )
              END IF
              CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), DPAR, L )
           ELSE
              CALL mexErrMsgTxt
     $             ( 'PARAM must not be used for this example' )
           END IF
        ELSE
           IF ( NRHS.EQ.IP )
     $        CALL mexErrMsgTxt( 'PARAM must be a real vector' )
        END IF
      END IF
C
C opt
C
      IF ( NRHS.GT.IP ) THEN
        IP = IP + 1
        IF ( mxGetM( PRHS(IP) ).NE.0 .AND. mxGetN( PRHS(IP) ).NE.0 )
     $  THEN
           L = 0
           IF ( DICO.EQ.1 ) THEN
              IF ( NR(1).EQ.2 .AND. NR(2).EQ.9 ) THEN
                 L = 2
              ELSE IF ( NR(1).EQ.3 .AND. NR(2).LE.2 .OR.
     $                  NR(1).EQ.4 .AND. NR(2).LE.4 ) THEN
                 L = 1
              END IF
           ELSE IF ( NR(1).EQ.4 .AND. NR(2).EQ.1 ) THEN
              L = 1
           END IF
           IF ( L.NE.0 ) THEN
              IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $             mxGetN( PRHS(IP) ).NE.1 ) THEN
                 CALL mexErrMsgTxt( 'OPT must be a scalar' )
              END IF
              IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $             mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
                 CALL mexErrMsgTxt( 'OPT must be an integer scalar' )
              END IF
              CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
              IPAR(1) = TEMP
              IF ( L.EQ.2 .AND. IPAR(1).NE.1 .AND. IPAR(1).NE.2 ) THEN
                 CALL mexErrMsgTxt
     $               ( 'OPT has 1 and 2 the only admissible values' )
              END IF
              IF ( NR(1).EQ.4 .AND. NR(2).EQ.4 .AND.
     $             IPAR(1).NE.NDEFC(4,4) ) THEN
C
C                Specify the file name; if it is '', the default name
C                is used.
C
                 IF ( NRHS.GT.IP ) THEN
                    IP = IP + 1
                    IF ( mxIsChar( PRHS(IP) ).NE.1 )
     $                 CALL mexErrMsgTxt( 'File name must be a string' )
                    L = mxGetN( PRHS(IP) )
                    IF ( L.EQ.0 ) THEN
                       CALL mexPrintf( 'Default file name is used' )
                       WRITE (CHPAR(1:11), '(A,I1,A,I1,A)') 'BB01',
     $                                      NR(1), '0', NR(2) , '.dat'
                       IPAR(4) = 11
                    ELSE
                       IE = mxGetString( PRHS(IP), CHPAR, L )
                       IF ( IE.NE.0 )
     $                    CALL mexErrMsgTxt
     $                               ( 'File name must be a string' )
                       IPAR(4) = L
                    END IF
                 ELSE
                    WRITE (CHPAR(1:11), '(A,I1,A,I1,A)') 'BB01',
     $                                   NR(1), '0', NR(2) , '.dat'
                    IPAR(4) = 11
                 END IF
              END IF
           ELSE
              CALL mexErrMsgTxt
     $             ( 'OPT must not be used for this example' )
           END IF
        ELSE
           CALL mexErrMsgTxt( 'OPT must be an integer scalar' )
        END IF
      END IF
C
      IF ( IP.GT.4 )
     $   DEF = 'N'
C
C Determine array dimensions.
C
      IF ( DICO.EQ.1 ) THEN
         N = NDEFC( NR(1), NR(2) )
         IF ( NR(1).LE.2  ) THEN
            M = MDEFC( NR(1), NR(2) )
            P = PDEFC( NR(1), NR(2) )
         ELSE
            IF ( NR(1).EQ.3 ) THEN
               IF ( NR(2).EQ.1 ) THEN
                  IF ( LSAME( DEF, 'D' ) ) THEN
                     M = N
                     P = N - 1
                     N = 2*N - 1
                  ELSE
                     M = IPAR(1)
                     P = IPAR(1) - 1
                     N = 2*IPAR(1) - 1
                  END IF
               ELSE IF ( NR(2).EQ.2 ) THEN
                  IF ( LSAME( DEF, 'N' ) )
     $               N = IPAR(1)
                  M = N
                  P = N
               END IF
            ELSE IF ( NR(1).EQ.4 ) THEN
               IF ( NR(2).LE.2 ) THEN
                  IF ( LSAME( DEF, 'N' ) )
     $               N = IPAR(1)
                  M = 1
                  P = 1
               ELSE IF ( NR(2).EQ.3 ) THEN
                  IF ( LSAME( DEF, 'D' ) ) THEN
                     N = 2*N
                     M = 2
                     P = N
                  ELSE
                     N = 2*IPAR(1)
                     M = 2
                     P = N
                  END IF
               ELSE IF ( NR(2).EQ.4 ) THEN
                  IF ( LSAME( DEF, 'D' ) ) THEN
                     M = N
                     P = N
                     N = 2*N - 1
                  ELSE
                     N = 2*IPAR(1) - 1
                     M = IPAR(1)
                     P = IPAR(1)
                  END IF
               END IF
            END IF
         END IF
      ELSE
         N = NDEFD( NR(1), NR(2) )
         IF ( NR(1).LE.2 ) THEN
            M = MDEFD( NR(1), NR(2) )
            P = PDEFD( NR(1), NR(2) )
         ELSE IF ( NR(1).EQ.4 ) THEN
            IF ( NR(2).EQ.1 ) THEN
               IF ( LSAME( DEF, 'D' ) ) THEN
                  N = 100
               ELSE
                  N = IPAR(1)
               END IF
               M = 1
               P = N
            END IF
         END IF
      END IF
C
C Determine the lengths of working arrays.
C
      LDA = N
      LDG = LDA
      LDQ = LDA
      LDB = LDA
      LDC = P
      LDX = LDA
C
      IF ( DICO.EQ.1 ) THEN
         LDWORK = N*MAX( 4, N )
      ELSE
         LDS = LDA
         IF ( .NOT.BPAR(4) ) LDG = M
         IF ( .NOT.BPAR(1) ) LDQ = P
         LDWORK = N*N
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE( A(LDA,N), B(LDB,M), C(LDC,N), G(LDG,N), Q(LDQ,N),
     $          X(LDX,N), DWORK(LDWORK) )
      IF ( DICO.EQ.2 )
     $   ALLOCATE( S(LDS,M) )
C
C Do the actual computations.
C
      IF ( DICO.EQ.1 ) THEN
         CALL BB01AD( DEF, NR, DPAR, IPAR, BPAR, CHPAR, VEC, N, M, P,
     $                A, LDA, B, LDB, C, LDC, G, LDG, Q, LDQ, X, LDX,
     $                DWORK, LDWORK, INFO )
      ELSE
         CALL BB02AD( DEF, NR, DPAR, IPAR, BPAR, CHPAR, VEC, N, M, P,
     $                A, LDA, B, LDB, C, LDC, Q, LDQ, G, LDG, S, LDS,
     $                X, LDX, DWORK, LDWORK, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         CALL mexPrintf( CHPAR )
         IF ( NLHS.GE.1 ) THEN
            IF ( ( DICO.EQ.1 .AND. VEC(9) ) .OR.
     $           ( DICO.EQ.2 .AND. VEC(10) ) ) THEN
               PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
               CALL mxCopyReal8ToPtr( X, mxGetPr( PLHS(1) ), N*N )
            ELSE
               PLHS(1) = mxCreateDoubleMatrix( 0, 0, 0 )
            END IF
         END IF
         IF ( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(2) ), N*N )
         END IF
         IF ( NLHS.GE.3 ) THEN
            IF ( FLAG(1).EQ.1 ) THEN
               PLHS(3) = mxCreateDoubleMatrix( N, N, 0 )
               CALL mxCopyReal8ToPtr( G, mxGetPr( PLHS(3) ), N*N )
            ELSE
               IF ( M.LT.LDG )
     $            CALL DLACPY( 'Full', M, M, G, LDG, G, M )
               PLHS(3) = mxCreateDoubleMatrix( M, M, 0 )
               CALL mxCopyReal8ToPtr( G, mxGetPr( PLHS(3) ), M*M )
            END IF
         END IF
         IF ( NLHS.GE.4 ) THEN
            IF ( FLAG(2).EQ.1 ) THEN
               PLHS(4) = mxCreateDoubleMatrix( N, N, 0 )
               CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(4) ), N*N )
            ELSE
               IF ( P.LT.LDQ )
     $            CALL DLACPY( 'Full', P, P, Q, LDQ, Q, P )
               PLHS(4) = mxCreateDoubleMatrix( P, P, 0 )
               CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(4) ), P*P )
            END IF
         END IF
         IP = 5
         IF ( ( NLHS.GE.IP ) .AND. ( FLAG(1).EQ.0 ) ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, M, 0 )
            CALL mxCopyReal8ToPtr( B, mxGetPr( PLHS(IP) ), N*M )
            IP = IP + 1
         END IF
         IF ( ( NLHS.GE.IP ) .AND. ( FLAG(2).EQ.0 ) ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( P, N, 0 )
            CALL mxCopyReal8ToPtr( C, mxGetPr( PLHS(IP) ), N*P )
            IP = IP + 1
         END IF
         IF ( NLHS.GE.IP .AND. DICO.EQ.2 .AND. BPAR(7) ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, M, 0 )
            CALL mxCopyReal8ToPtr( S, mxGetPr( PLHS(IP) ), N*M )
            IP = IP + 1
         END IF
         IF ( NLHS.GE.IP ) THEN
            IF ( DICO.EQ.1 ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( 7, 1, 0 )
               CALL mxCopyReal8ToPtr( DPAR, mxGetPr( PLHS(IP) ), 7 )
            ELSE
               PLHS(IP) = mxCreateDoubleMatrix( 4, 1, 0 )
               CALL mxCopyReal8ToPtr( DPAR, mxGetPr( PLHS(IP) ), 4 )
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( A, B, C, G, Q, X, DWORK )
      IF ( DICO.EQ.2 )
     $   DEALLOCATE( S )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( DICO.EQ.1 ) THEN
            WRITE( TEXT,'( " INFO = ", I4, " ON EXIT FROM BB01AD" )' )
     $             INFO
         ELSE
            WRITE( TEXT,'( " INFO = ", I4, " ON EXIT FROM BB02AD" )' )
     $             INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C
C *** Last line of AREBENCH ***
      END
