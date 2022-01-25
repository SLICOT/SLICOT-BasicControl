#include "fintrf.h"
C ARECOND.F - Gateway function for estimating the conditioning of
C             Lyapunov and algebraic Riccati equations using SLICOT
C             routines SB03QD, SB03SD, SB02QD, and SB02SD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   If fact = 0,
C   [(rcond(,sep))(,ferr)(,T(,U))] =
C            arecond(eq,job,reduced,fact,A,C,X(,flag))       |eq| = 1;
C   [(rcond(,sep))(,ferr)(,T(,U))] =
C            arecond(eq,job,reduced,fact,A,Q,G,X(,flag))     |eq| = 2;
C   If fact = 1 and reduced = 1,
C   [(rcond(,sep))(,ferr)] =
C            arecond(eq,job,reduced,fact,T,C,X(,flag))       |eq| = 1;
C   [(rcond(,sep))(,ferr)] =
C            arecond(eq,job,reduced,fact,T,Q,G,X(,flag))     |eq| = 2;
C   If fact = 1 and reduced = 0,
C   [(rcond(,sep))(,ferr)] =
C            arecond(eq,job,reduced,fact,A,T,U,C,X(,flag))   |eq| = 1;
C   [(rcond(,sep))(,ferr)] =
C            arecond(eq,job,reduced,fact,A,T,U,Q,G,X(,flag)) |eq| = 2.
C
C Purpose:
C   1. To estimate the conditioning and compute an error bound on the
C   solution of the real continuous-time Lyapunov matrix equation
C
C         op(A)'*X + X*op(A) = C,                                    (1)
C
C   or of the real discrete-time Lyapunov matrix equation
C
C         op(A)'*X*op(A) - X = C,                                    (2)
C
C   where op(A) = A or A' (A**T) and C is symmetric (C = C').
C
C   2. To estimate the conditioning and compute an error bound on the
C   solution of the real continuous-time matrix algebraic Riccati
C   equation
C
C         op(A)'*X + X*op(A) + Q - X*G*X = 0,                        (3)
C
C   or of the real discrete-time matrix algebraic Riccati equation
C                                 -1
C         X = op(A)'*X*(I_n + G*X)  *op(A) + Q,                      (4)
C
C   where Q and G are symmetric. Let Ac denote the closed-loop matrix,
C   hence, in the continuous-time case,
C
C         Ac = A - G*X,          if op(A) = A,  or
C         Ac = A - X*G,          if op(A) = A',
C
C   and, in the discrete-time case,
C
C        Ac = inv(I_n + G*X)*A,  if op(A) = A,  or
C        Ac = A*inv(I_n + X*G),  if op(A) = A'.
C
C   Riccati equations (3) and (4) are equivalent to the following
C   reduced equations, respectively,
C               _   _         _   _ _ _
C        op(T)'*X + X*op(T) + Q + X*G*X = 0,
C        _          _                _ _ _         _
C        X = op(T)'*X*op(T) + op(T)'*X*G*X*op(T) + Q,
C         _           _               _
C   where X = U'*X*U, Q = U'*Q*U, and G = U'*G*U, with U the
C   orthogonal matrix reducing Ac to a real Schur form, T = U'*Ac*U.
C   Similar, simpler formulas stand for Lyapunov equations (1) and (2).
C
C Input parameters:
C   eq     - integer option to indicate the equation type:
C            = -1 : continuous-time Lyapunov equation (1);
C            =  1 : discrete-time Lyapunov equation (2);
C            = -2 : continuous-time Riccati equation (3);
C            =  2 : discrete-time Riccati equation (4).
C   job    - integer option to indicate the calculation to be performed:
C            =  1 : compute the reciprocal condition number and
C                   separation;
C            =  2 : compute the error bound only;
C            =  3 : compute the reciprocal condition number, the
C                   separation, and the error bound.
C   reduced- integer option specifying whether or not the original
C            Lyapunov equations should be solved in the iterative
C            estimation process (also for (3) and (4)), as follows:
C            =  0 :  solve the original Lyapunov equations, updating
C                    the right-hand sides and solutions with the
C                    matrix U, e.g., RHS <-- U'*RHS*U;
C            =  1 :  solve reduced Lyapunov equations only, without
C                    updating the right-hand sides and solutions.
C                    This scheme is faster, but sometimes slightly
C                    less accurate.
C   fact   - integer option specifying whether or not the real Schur
C            factorization of the matrix A or Ac is supplied on entry,
C            as follows:
C            =  0 :  the Schur factorization of A (if |eq| = 1) or Ac
C                    (if |eq| = 2) will be computed and the factors can
C                    be stored in the output parameters T and U (if
C                    reduced = 0);
C            =  1 :  the input parameters T and U (if reduced = 0)
C                    contain the factors from the real Schur
C                    factorization of the matrix A (if |eq| = 1) or Ac
C                    (if |eq| = 2).
C   A      - if fact = 0 or reduced = 0, the real n-by-n system state
C            matrix A.
C   T      - if fact = 1, a real n-by-n Schur form of A, if |eq| = 1,
C            or Ac, if |eq| = 2.
C   U        If fact = 1 and reduced = 0, the n-by-n orthogonal matrix
C            U which reduced the matrix A or Ac to a real Schur form T.
C   C      - real symmetric n-by-n right-hand side matrix of the
C            original equation, if reduced = 0, or of the reduced
C            equation (with matrix T), if reduced = 1.
C   X      - real symmetric n-by-n solution matrix of original equation,
C            if reduced = 0, or of the reduced equation (with matrix T),
C            if reduced = 1.
C   Q      - real symmetric n-by-n state weighting matrix of the Riccati
C            equation. Matrix Q must correspond to that of the reduced
C            equation, if reduced = 1.
C   G      - real symmetric n-by-n matrix of the Riccati equation.
C            Matrix G must correspond to that of the reduced equation,
C            if reduced = 1.
C   flag   - (optional) vector containing options.
C            flag(1) specifies the form of op(A) to be used, as follows:
C                    = 0 : op(A) = A;
C                    = 1 : op(A) = A'.
C            flag(2) specifies which part of the symmetric matrix C
C                    (if |eq| = 1), or the symmetric matrices Q and G
C                    (if |eq| = 2), is to be used, as follows:
C                    = 0 : upper triangular part;
C                    = 1 : lower triangular part.
C            Default:      flag(1:2) = [0,0].
C
C Output parameters:
C   rcond  - (optional) if job = 1 or job = 3, an estimate of the
C            reciprocal condition number of the Lyapunov or Riccati
C            equation.
C   sep    - (optional) if job = 1 or job = 3, an estimate of the
C            separation, i.e., of the quantity
C               sep(op(A),-op(A)'),   if eq = -1;
C               sepd(op(A),op(A)'),   if eq =  1;
C               sep(op(Ac),-op(Ac)'), if eq = -2;
C               sepd(op(Ac),op(Ac)'), if eq =  2.
C   ferr   - (optional) if job = 2 or job = 3, an estimated forward
C            error bound for the solution X. If Xtrue is the true
C            solution, ferr bounds the magnitude of the largest entry in
C            (X - Xtrue) divided by the magnitude of the largest entry
C            in X.
C   T      - (optional) if fact = 0, the n-by-n real Schur form matrix
C            corresponding to A or Ac.
C   U      - (optional) if fact = 0 and reduced = 0, the n-by-n real
C            orthogonal matrix which reduced the matrix A or Ac to the
C            real Schur form T.
C
C Contributor:
C   V. Sima, Katholieke Univ. Leuven, Belgium, March 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Feb. 2005,
C   Apr. 2009, July 2012, Jul. 2016.
C
C **********************************************************************
C
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C
C .. Mex-file interface parameters ..
      mwPointer         PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      mwPointer         mxCalloc, mxCreateDoubleMatrix,
     $                  mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         FACT, JOB, LYAPUN, TRANA, UPLO
      INTEGER           INFO, LDA, LDC, LDG, LDQ, LDT, LDU, LDWORK, LDX,
     $                  N
      DOUBLE PRECISION  FERR, RCOND, SCALE, SEP
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      mwPointer         A, C, DWORK, G, Q, T, U, X
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           IWARN, JOBC, JOBE, NOFACT, UPDATE
      CHARACTER*1       MATN(6)
      CHARACTER*6       SUBNAM
      CHARACTER*120     TEXT
      INTEGER           FLAG(2), I, IA, IAEQ, IC, IEQ, IG, IP, IPAR,
     $                  IT, ITMP, IU, IX, LIWORK, LWA, NN
      DOUBLE PRECISION  FLAGR(2), TEMP
C
C .. External functions ..
      DOUBLE PRECISION  DLANSY
      EXTERNAL          DLANSY
C
C .. External subroutines ..
      EXTERNAL          SB02QD, SB02SD, SB03QD, SB03SD
C
C ..Intrinsic functions..
      INTRINSIC         ABS, MAX
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.7 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ARECOND requires at least 7 input arguments' )
      ELSE IF ( NLHS.GT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ARECOND requires at most 5 output arguments' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   eq
C
      IPAR = 7
      IF ( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 )
     $      THEN
         CALL mexErrMsgTxt( 'EQ must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'EQ must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      IEQ  = TEMP
      IAEQ = ABS( IEQ )
      IF ( IAEQ.GT.2 .OR. IAEQ.EQ.0 ) THEN
         CALL mexErrMsgTxt
     $        ( 'EQ has -2, -1, 1, or 2 the only admissible values' )
      END IF
      IF ( IAEQ.EQ.2 )
     $   IPAR = IPAR + 1
C
C   job
C
      IF ( mxGetM( PRHS(2) ).NE.1 .OR. mxGetN( PRHS(2) ).NE.1 )
     $      THEN
         CALL mexErrMsgTxt( 'JOB must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), TEMP, 1 )
      ITMP = TEMP
      IF ( ITMP.LT.1 .OR. ITMP.GT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'JOB has 1, 2, or 3 the only admissible values' )
      END IF
      JOBC = ITMP.EQ.1
      JOBE = ITMP.EQ.2
      IF ( JOBC ) THEN
         JOB = 'C'
      ELSE IF ( JOBE ) THEN
         JOB = 'E'
      ELSE
         JOB = 'B'
      END IF
C
C   reduced
C
      IF ( mxGetM( PRHS(3) ).NE.1 .OR. mxGetN( PRHS(3) ).NE.1 )
     $      THEN
         CALL mexErrMsgTxt( 'REDUCED must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'REDUCED must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TEMP, 1 )
      ITMP = TEMP
      IF ( ITMP.LT.0 .OR. ITMP.GT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'REDUCED has 0 or 1 the only admissible values' )
      END IF
      UPDATE = ITMP.EQ.0
      IF ( UPDATE ) THEN
         LYAPUN = 'O'
      ELSE
         LYAPUN = 'R'
      END IF
C
C   fact
C
      IF ( mxGetM( PRHS(4) ).NE.1 .OR. mxGetN( PRHS(4) ).NE.1 )
     $      THEN
         CALL mexErrMsgTxt( 'FACT must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'FACT must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), TEMP, 1 )
      ITMP = TEMP
      IF ( ITMP.LT.0 .OR. ITMP.GT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'FACT has 0 or 1 the only admissible values' )
      END IF
      NOFACT = ITMP.EQ.0
      IF ( NOFACT ) THEN
         FACT = 'N'
      ELSE
         FACT = 'F'
         IF ( UPDATE )
     $      IPAR = IPAR + 2
      END IF
C
      IF ( NRHS.LT.IPAR ) THEN
         WRITE( TEXT, '(''ARECOND requires at least '', I2,
     $                  '' input arguments'')' ) IPAR
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      IA = 5
      IF ( .NOT.NOFACT ) THEN
         IF ( UPDATE ) THEN
            IT = 6
            IU = 7
            MATN(1) = 'A'
            MATN(2) = 'T'
            MATN(3) = 'U'
         ELSE
            IT = 5
            IU = 5
            MATN(1) = 'T'
         END IF
      ELSE
         IT = 5
         IU = 5
         MATN(1) = 'A'
      END IF
      IC = IU + 1
      IF ( IAEQ.EQ.1 ) THEN
         MATN(IC-4) = 'C'
      ELSE
         IG = IC + 1
         MATN(IC-4) = 'Q'
         MATN(IC-3) = 'G'
      END IF
      IX = IC + IAEQ
      IP = IX + 1
      MATN(IX-4) = 'X'
C
C   A(NxN), (T(NxN), (U(NxN),)) C/Q(NxN), (G(NxN),) X(NxN)
C
      N = mxGetM( PRHS(IA) )
C
      IF ( mxGetN( PRHS(IA) ).NE.N ) THEN
         WRITE( TEXT, '(A1,'' must be a square matrix'')' ) MATN(IA-4)
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( mxIsNumeric( PRHS(IA) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(IA) ).EQ.1 ) THEN
         WRITE( TEXT, '(A1,'' must be a real matrix'')' ) MATN(IA-4)
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      DO 10 I = IA + 1, IX
         IF ( mxGetM( PRHS(I) ).NE.N .OR.
     $        mxGetN( PRHS(I) ).NE.N ) THEN
            WRITE( TEXT, '(A1,'' must have the same size as '', A1 )' )
     $          MATN(I-4), MATN(IA-4)
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(I) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(I) ).EQ.1 ) THEN
            WRITE( TEXT, '(A1,'' must be a real matrix'')' ) MATN(I-4)
            CALL mexErrMsgTxt( TEXT )
         END IF
   10 CONTINUE
C
      FLAG(1) = 0
      FLAG(2) = 0
C
      IF ( NRHS.GE.IP ) THEN
C
C   flag
C
         ITMP = mxGetM( PRHS(IP) ) * mxGetN( PRHS(IP) )
         IF ( ITMP.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'FLAG must be a vector with at most 2 elements' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'FLAG must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), FLAGR, ITMP )
         IF ( ITMP.GT.0 ) FLAG(1) = FLAGR(1)
         IF ( ITMP.GT.1 ) FLAG(2) = FLAGR(2)
      END IF
C
      IF ( FLAG(1).EQ.0 ) THEN
         TRANA = 'N'
      ELSE
         TRANA = 'T'
      END IF
C
      IF ( FLAG(2).EQ.0 ) THEN
         UPLO = 'U'
      ELSE
         UPLO = 'L'
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms
C in DGEES.
C
      LDA = MAX( 1, N )
      LDC = LDA
      LDQ = LDA
      LDG = LDA
      LDT = LDA
      LDU = LDA
      LDX = LDA
C
      NN = N*N
      IF( IEQ.EQ.-1 ) THEN
         IF( JOBC ) THEN
            LDWORK = 2*NN
         ELSE
            LDWORK = 3*NN
         END IF
         IF( .NOT.( JOBC .OR. UPDATE ) )
     $      LDWORK = LDWORK + N - 1
         IF( NOFACT )
     $      LDWORK = MAX( LDWORK, 5*N )
C
      ELSE IF( IEQ.EQ.1 ) THEN
         LDWORK = MAX( 3, 2*NN ) + NN
         IF( JOBC .AND. NOFACT ) THEN
            LDWORK = MAX( LDWORK, 5*N )
         ELSE IF( .NOT.JOBC ) THEN
            LDWORK = LDWORK + 2*N
         END IF
C
      ELSE IF( IEQ.EQ.-2 ) THEN
         IF( UPDATE .AND. .NOT.JOBC ) THEN
            LWA = NN
         ELSE
            LWA = 0
         END IF
C
         IF( NOFACT ) THEN
            IF( JOBC ) THEN
               LDWORK = MAX( 5*N, 2*NN )
            ELSE
               LDWORK = MAX( LWA + 5*N, 4*NN )
            END IF
         ELSE
            IF( JOBC ) THEN
               LDWORK = 2*NN
            ELSE
               LDWORK = 4*NN
            END IF
         END IF
C
      ELSE IF( IEQ.EQ.2 ) THEN
         IF( UPDATE ) THEN
            LWA = NN
         ELSE
            LWA = 0
         END IF
C
         LDWORK = MAX( 3, 2*NN ) + NN
         IF( .NOT.JOBC ) THEN
            LDWORK = LDWORK + NN
            IF( .NOT.UPDATE )
     $         LDWORK = LDWORK + N
         END IF
         IF( NOFACT )
     $      LDWORK = MAX( LWA + 5*N, LDWORK )
      END IF
C
      LDWORK = MAX( 1, LDWORK )
      LIWORK = NN
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE( IWORK(LIWORK) )
      DWORK = mxCalloc( LDWORK, 8 )
      T     = mxCalloc(  LDT*N, 8 )
      X     = mxCalloc(  LDX*N, 8 )
C
      IF ( UPDATE .OR. NOFACT ) THEN
         A = mxCalloc( LDA*N, 8 )
         IF ( UPDATE )
     $      U = mxCalloc( LDU*N, 8 )
      END IF
C
      IF ( IAEQ.EQ.1 ) THEN
         C = mxCalloc( LDC*N, 8 )
      ELSE
         G = mxCalloc( LDG*N, 8 )
         Q = mxCalloc( LDQ*N, 8 )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      IF ( UPDATE .OR. NOFACT ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IA) ), %VAL( A ),
     $                          N*N )
         IF ( .NOT.NOFACT )
     $      CALL mxCopyPtrToReal8( mxGetPr( PRHS(IU) ), %VAL( U ),
     $                             N*N )
      END IF
C
      IF ( .NOT.NOFACT )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(IT) ), %VAL( T ),
     $                          N*N )
C
      IF ( IAEQ.EQ.1 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IC) ), %VAL( C ),
     $                          N*N )
         SCALE = ONE
      ELSE
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IC) ), %VAL( Q ),
     $                          N*N )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IG) ), %VAL( G ),
     $                          N*N )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(IX) ), %VAL( X ), N*N )
C
C Do the actual computations.
C
      IF ( IEQ.EQ.-1 ) THEN
         CALL SB03QD( JOB, FACT, TRANA, UPLO, LYAPUN, N, SCALE,
     $                %VAL( A ), LDA, %VAL( T ), LDT, %VAL( U ), LDU,
     $                %VAL( C ), LDC, %VAL( X ), LDX, SEP, RCOND, FERR,
     $                IWORK, %VAL( DWORK ), LDWORK, INFO )
         SUBNAM = 'SB03QD'
C
      ELSE IF ( IEQ.EQ.1 ) THEN
         CALL SB03SD( JOB, FACT, TRANA, UPLO, LYAPUN, N, SCALE,
     $                %VAL( A ), LDA, %VAL( T ), LDT, %VAL( U ), LDU,
     $                %VAL( C ), LDC, %VAL( X ), LDX, SEP, RCOND, FERR,
     $                IWORK, %VAL( DWORK ), LDWORK, INFO )
         SUBNAM = 'SB03SD'
C
      ELSE IF ( IEQ.EQ.-2 ) THEN
         CALL SB02QD( JOB, FACT, TRANA, UPLO, LYAPUN, N, %VAL( A ), LDA,
     $                %VAL( T ), LDT, %VAL( U ), LDU, %VAL( G ), LDG,
     $                %VAL( Q ), LDQ, %VAL( X ), LDX, SEP, RCOND, FERR,
     $                IWORK, %VAL( DWORK ), LDWORK, INFO )
         SUBNAM = 'SB02QD'
C
      ELSE
         CALL SB02SD( JOB, FACT, TRANA, UPLO, LYAPUN, N, %VAL( A ), LDA,
     $                %VAL( T ), LDT,  %VAL( U ), LDU, %VAL( G ), LDG,
     $                %VAL( Q ), LDQ, %VAL( X ), LDX, SEP, RCOND, FERR,
     $                IWORK, %VAL( DWORK ), LDWORK, INFO )
         SUBNAM = 'SB02SD'
C
      END IF
C
      IWARN = INFO.EQ.N+1
      IF ( IWARN ) THEN
         WRITE( TEXT,
     $      '('' Warning: INFO = '', I4, '' ON EXIT FROM '', A6)' )
     $      INFO, SUBNAM
         INFO = 0
      ELSE IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '', I4, '' ON EXIT FROM  '', A6)' )
     $      INFO, SUBNAM
         GO TO 20
      END IF
C
      IF ( N.EQ.0 .OR.
     $     DLANSY( '1-norm', UPLO, N, %VAL( X ), LDX, FLAGR ).EQ.ZERO )
     $   SEP = ZERO
C
C Copy output to MATLAB workspace.
C
      IP = 1
      IF ( .NOT.JOBE ) THEN
         PLHS(1) = mxCreateDoubleMatrix( 1, 1, 0 )
         CALL mxCopyReal8ToPtr( RCOND, mxGetPr( PLHS(1) ), 1 )
         IF ( NLHS.GT.1 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( SEP, mxGetPr( PLHS(2) ), 1 )
            IP = 3
         END IF
      END IF
C
      IF ( .NOT.JOBC .AND. NLHS.GE.IP ) THEN
         PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
         CALL mxCopyReal8ToPtr( FERR, mxGetPr( PLHS(IP) ), 1 )
         IP = IP + 1
      END IF
      IF ( NOFACT .AND. NLHS.GE.IP ) THEN
         PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
         CALL mxCopyReal8ToPtr( %VAL( T ), mxGetPr( PLHS(IP) ),
     $                          N*N )
         IP = IP + 1
         IF ( UPDATE .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( U ), mxGetPr( PLHS(IP) ),
     $                             N*N )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
   20 CONTINUE
C
      DEALLOCATE(  IWORK )
      CALL mxFree( DWORK )
      CALL mxFree( T     )
      CALL mxFree( X     )
      IF ( UPDATE .OR. NOFACT ) THEN
         CALL mxFree( A )
         IF ( UPDATE )
     $      CALL mxFree( U )
      END IF
C
      IF ( IAEQ.EQ.1 ) THEN
         CALL mxFree( C )
      ELSE
         CALL mxFree( G )
         CALL mxFree( Q )
      END IF
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
C *** Last line of ARECOND ***
      END
