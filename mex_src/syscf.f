C SYSCF.F - Gateway function for SLICOT coprime factorization routines
C           SB08CD.F, SB08DD.F, SB08ED.F and SB08FD.F.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C MATLAB CALL:
C   [Af,Bf,Cf,Df,nr] = syscf(meth,A,B,C,D,tol,discr,alpha,beta)
C
C Purpose:
C   To compute for a state-space system (A,B,C,D) with the
C   transfer-function matrix G a left coprime factorization (LCF) or
C   a right coprime factorization (RCF)
C                   -1                    -1
C         LCF: G = M  N,  or  RCF:  G = NM  ,
C
C   as specified by the method parameter meth.
C   The resulting (Af,Bf,Cf,Df) is the state-space system having
C   the transfer-function matrix [ N M ] for a LCF or [ N ] for a RCF.
C                                                     [ M ]
C   A minimal realization of M of order nr can be recovered from the
C   resulting state-space matrices (see description of nr).
C
C Input parameters:
C   meth    - type of right/left coprime factorization (LCF/RCF):
C               = 1 : RCF with inner denominator;
C               = 2 : LCF with inner denominator;
C               = 3 : RCF with ALPHA stability degree;
C               = 4 : LCF with ALPHA stability degree.
C   A,B,
C   C,D     - state-space system matrices.
C   tol     - (optional) controllability/observability tolerance
C             for computing right/left coprime factorizations.
C                Default controllability tolerance
C                   tol = epsilon_machine*max(norm(A),norm(B));
C                default  observability tolerance
C                   tol = epsilon_machine*max(norm(A),norm(C)).
C   discr   - (optional) type of system:
C                = 0 : continuous-time (default);
C                = 1 : discrete-time.
C   alpha   - (optional) stability degree for LCF/RCF with prescribed
C                stability degree.
C                Default: -0.05  for continuous-time;
C                          0.95  for discrete-time.
C   beta    - (optional) stability margin for the LCF/RCF with prescribed
C                stability degree.
C                Default: beta = alpha.
C
C Output parameters:
C   Af, Bf,                                   ( N )                  -1
C   Cf, Df  - matrices of the compound system ( M ) such that G = N*M
C             is a RCF  if meth = 1 or meth = 3, or
C             matrices of the compound system ( N M ) such that
C                  -1
C             G = M  *N is a LCF  if meth = 2 or meth = 4.
C   nr      - order of a minimal realization of M;
C             for a p-by-m transfer-function matrix G, the minimal
C             realization of M can be recovered as the state-space model
C               ( Af(si,si), Bf(si,ii), Cf(oi,si), Df(oi,ii) )
C             where si, ii and oi are the index vectors for state, input
C             and output vector, respectively, defined as follows
C                si = 1:nr, ii = m+1:m+p, oi = 1:p for a LCF;
C                si = nr+1:end; ii = 1:m, oi = p+1:p+m for a RCF.
C
C Contributor:
C A. Varga, DLR-Oberpfaffenhofen; March 1999.
C
C Revisions:
C A. Varga, DLR-Oberpfaffenhofen; March 2003, Sept. 2005.
C V. Sima, Research Institute for Informatics, Nov. 2005, Apr. 2009,
C Dec. 2012.
C
C *********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ONE, TWO, ZERO, P05
      PARAMETER         ( ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0,
     $                    P05 = 0.05D0 )
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
      INTEGER           INFO, IWARN, LDA, LDB, LDBR, LDC, LDD, LDDR,
     $                  LWORK, M, N, NQ, NR, P
      DOUBLE PRECISION  ALPHA(2), TOL
C
C .. Allocatable local arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE:: A(:,:), B(:,:), C(:,:), D(:,:),
     $                                BR(:,:), DR(:,:), WORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           DISCR, LEFT
      INTEGER           I, LW1, LW2, LW3, LW4, METH, MP, M1,
     $                  N1, N2, N3, PM, P1
      DOUBLE PRECISION  DUM
C
C .. External functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C
C .. External Subroutines ..
      EXTERNAL          DLACPY, SB08CD, SB08DD, SB08ED, SB08FD
C
C .. Intrinsic functions ..
      INTRINSIC         DBLE, MAX
C
C Check for proper number of arguments
C
      IF( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSCF requires at least 5 input arguments' )
      ELSE IF( NLHS.GT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSCF has at most 5 output argument' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters
C
C  meth
C
      METH = 1
      IF( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be a scalar' )
      END IF
      IF( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be an integer scalar 0 or 1' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), DUM, 1 )
      METH = DUM
      IF( METH.LE.0 .OR. METH.GT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'METH has 1 ... 4 the only admissible values' )
      END IF
      LEFT = METH.EQ.2 .OR. METH.EQ.4
C
C   A(NxN), B(NxM), C(PxN), D(PxM)
C
      N  = mxGetM( PRHS(2) )
      M  = mxGetN( PRHS(3) )
      P  = mxGetM( PRHS(4) )
      N1 = mxGetN( PRHS(2) )
      N2 = mxGetM( PRHS(3) )
      N3 = mxGetN( PRHS(4) )
      P1 = mxGetM( PRHS(5) )
      M1 = mxGetN( PRHS(5) )
C
      IF( (N1.NE.N) ) THEN
         CALL mexErrMsgTxt
     $        ( 'A must be a square matrix' )
      END IF
      IF( N2.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'B must have the same row dimension as A' )
      END IF
      IF( N3.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'C must have the same column dimension as A ' )
      END IF
      IF( M1.NE.M ) THEN
         CALL mexErrMsgTxt
     $        ( 'D must have the same column dimension as B' )
      END IF
      IF( P1.NE.P ) THEN
         CALL mexErrMsgTxt
     $        ( 'D must have the same row dimension as C' )
      END IF
      IF( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a real matrix' )
      END IF
C
C   tol
C
      TOL = ZERO
      IF( NRHS.GT.5 ) THEN
         IF( mxGetM( PRHS(6) ).NE.1 .OR. mxGetN( PRHS(6) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), TOL, 1 )
      END IF
C
C   discr
C
      DICO = 'C'
      IF( NRHS.GT.6 ) THEN
         IF( mxGetM( PRHS(7) ).NE.1 .OR. mxGetN( PRHS(7) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'DISCR must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(7) ).EQ.1 ) THEN
           CALL mexErrMsgTxt( 'DISCR must be an integer scalar 0 or 1' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), DUM, 1 )
         IF( DUM.NE.ZERO ) DICO = 'D'
      END IF
      DISCR = LSAME( DICO, 'D' )
C
C   alpha
C
      IF( NRHS.GT.7 ) THEN
         IF( mxGetM( PRHS(8) ).NE.1 .OR. mxGetN( PRHS(8) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(8) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(8) ), ALPHA(1), 1 )
      ELSE
         ALPHA(1) = -P05
         IF( DISCR ) ALPHA(1) = ONE + ALPHA(1)
      END IF
C
C   beta
C
      IF( NRHS.GT.8 ) THEN
         IF( mxGetM( PRHS(9) ).NE.1 .OR. mxGetN( PRHS(9) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'BETA must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(9) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(9) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'BETA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(9) ), ALPHA(2), 1 )
      ELSE
         ALPHA(2) = ALPHA(1)
      END IF
C
C Determine the lenghts of working arrays
C
      LW1   = MAX( N*(N+5), 5*M, 4*P )
      LW2   = N*P + MAX( N*(N+5), 5*P, 4*M )
      LW3   = MAX( N*(N+5), M*(M+2), 4*M, 4*P )
      LW4   = P*N + MAX( N*(N+5),P*(P+2),4*P,4*M )
      LWORK = MAX( 1, LW1, LW2, LW3, LW4 )
C
      LDA  = MAX( 1, N )
      LDB  = MAX( 1, N )
      LDC  = MAX( 1, P+M )
      LDD  = MAX( 1, P+M )
      LDBR = LDB
      LDDR = MAX( 1, P )
C
C Allocate variable dimension local arrays
C !Fortran 90/95
C
      ALLOCATE ( A( LDA, N ), B( LDB, M+P ), C( LDC, N ), D( LDD, M+P ),
     $           BR( LDBR, P ), DR( LDDR, P ), WORK( LWORK ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), WORK, P*N )
      CALL DLACPY( 'F', P, N, WORK, MAX( 1, P ), C, LDC )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), WORK, P*M )
      CALL DLACPY( 'F', P, M, WORK, MAX( 1, P ), D, LDD )
C
C Do the actual computations
C
C RCFID
C
      IF( METH.EQ.1 )
     $   CALL SB08DD( DICO, N, M, P, A, LDA, B, LDB, C, LDC,
     $                D, LDD, NQ, NR, C(P+1,1), LDC, D(P+1,1), LDD,
     $                TOL, WORK, LWORK, IWARN, INFO )
C
C LCFID
C
      IF( METH.EQ.2 )
     $   CALL SB08CD( DICO, N, M, P, A, LDA, B, LDB, C, LDC,
     $                D, LDD, NQ, NR, BR, LDBR, DR, LDDR, TOL, WORK,
     $                LWORK, IWARN, INFO )
C
C RCFS
C
      IF( METH.EQ.3 )
     $   CALL SB08FD( DICO, N, M, P, ALPHA, A, LDA, B, LDB, C, LDC,
     $                D, LDD, NQ, NR, C(P+1,1), LDC, D(P+1,1), LDD,
     $                TOL, WORK, LWORK, IWARN, INFO )
C
C LCFS
C
      IF( METH.EQ.4 )
     $   CALL SB08ED( DICO, N, M, P, ALPHA, A, LDA, B, LDB, C, LDC,
     $                D, LDD, NQ, NR, BR, LDBR, DR, LDDR, TOL, WORK,
     $                LWORK, IWARN, INFO )
C
C Copy output to MATLAB workspace
C
      IF( INFO.EQ.0 ) THEN
         IF( LEFT ) THEN
            PM = P
            MP = M + P
         ELSE
            PM = P + M
            MP = M
         END IF
         IF( NLHS.GE.1 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( NQ, NQ, 0 )
            IF( NQ.LT.N .AND. NQ.GT.0 )
     $         CALL DLACPY( 'F', NQ, NQ, A, LDA, A, NQ )
            CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(1) ), NQ*NQ )
         END IF
         IF( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( NQ, MP, 0 )
            CALL DLACPY( 'F', NQ, M, B, LDB, WORK, NQ )
            IF( LEFT )
     $         CALL DLACPY( 'F', NQ, P, BR, LDBR, WORK(NQ*M+1), NQ )
            CALL mxCopyReal8ToPtr( WORK, mxGetPr( PLHS(2) ), NQ*MP )
         END IF
         IF( NLHS.GE.3 ) THEN
            PLHS(3) = mxCreateDoubleMatrix( PM, NQ, 0 )
            CALL DLACPY( 'F', PM, NQ, C, LDC, WORK, PM )
            CALL mxCopyReal8ToPtr( WORK, mxGetPr( PLHS(3) ), PM*NQ )
         END IF
         IF( NLHS.GE.4 ) THEN
            PLHS(4) = mxCreateDoubleMatrix( PM, MP, 0 )
            CALL DLACPY( 'F', PM, M, D, LDD, WORK, PM )
            IF( LEFT )
     $         CALL DLACPY( 'F', P, P, DR, LDDR, WORK(PM*M+1), PM )
            CALL mxCopyReal8ToPtr( WORK, mxGetPr( PLHS(4) ), PM*MP )
         END IF
C
         IF( NLHS.GE.5 ) THEN
            PLHS(5) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( DBLE( NR ), mxGetPr( PLHS(5) ), 1 )
         END IF
C
      END IF
C
C Deallocate local arrays
C !Fortran 90/95
C
      DEALLOCATE ( A, B, C, D, BR, DR, WORK )
C
C Error handling
C
      IF( INFO.NE.0 ) THEN
         IF( METH.EQ.1 ) THEN
            WRITE( TEXT,'( " INFO =", I4, " ON EXIT FROM SB08DD" )' )
     $             INFO
         ELSE IF( METH.EQ.2 ) THEN
            WRITE( TEXT,'( " INFO =", I4, " ON EXIT FROM SB08CD" )' )
     $             INFO
         ELSE IF( METH.EQ.3 ) THEN
            WRITE( TEXT,'( " INFO =", I4, " ON EXIT FROM SB08FD" )' )
     $             INFO
         ELSE IF( METH.EQ.4 ) THEN
            WRITE( TEXT,'( " INFO =", I4, " ON EXIT FROM SB08ED" )' )
     $             INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SYSCF ***
      END
