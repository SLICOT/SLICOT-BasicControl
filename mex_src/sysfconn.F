#include "fintrf.h"
C SYSFCONN.F - Gateway function to compute, for a given state-space
C              system (A,B,C,D), the closed-loop system (Ac,Bc,Cc,Dc)
C              corresponding to the output, or mixed output and state,
C              feedback control law, using SLICOT routines AB05RD and
C              AB05SD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [Ao,Bo,Co(,Do,rcondi)] = SYSFCONN(task,jobd,fbtype,alpha(,beta)
C                                     A,B,C(,D)(,G,H)(,F)(,K))
C
C   [Ao,Bo,Co(,Do,rcondi)] = SYSFCONN(1,jobd,fbtype,alpha,beta,A,B,C(,D)
C                                     G,H(,F)(,K))
C   [Ao,Bo,Co(,Do,rcondi)] = SYSFCONN(2,jobd,fbtype,alpha,A,B,C(,D)(,F))
C
C Purpose:
C   To compute a state-space model (A,B,C,D) for the various feedback
C   inter-connections of two systems, according to the value of task:
C
C   task = 1: To construct, for a given state space system (A,B,C,D),
C   the closed-loop system (Ac,Bc,Cc,Dc) corresponding to the mixed
C   output and state feedback control law
C
C          u = alpha*F*y + beta*K*x + G*v,
C          z = H*y.                                                  (1)
C
C   task = 2: To construct, for a given state space system (A,B,C,D),
C   the closed-loop system (Ac,Bc,Cc,Dc) corresponding to the output
C   feedback control law
C
C          u = alpha*F*y + v.                                        (2)
C
C Input parameters:
C   task   - integer specifying the computations to be performed.
C            task = 1 :  compute the closed-loop system corresponding to
C                        a mixed output and state feedback control law;
C            task = 2 :  compute the closed-loop system corresponding to
C                        an output feedback control law.
C   jobd   - integer specifying whether or not a non-zero matrix D
C            appears in the given state space model:
C            jobd = 0 :  D is assumed a zero matrix;
C            jobd = 1 :  D is present.
C   fbtype - integer specifying the type of the feedback law as follows:
C            fbtype = 1 :  Unitary output feedback (F = I);
C            fbtype = 2 :  General output feedback.
C   alpha  - the real coefficient alpha in the output feedback law.
C   beta   - the real coefficient beta in the state feedback law.
C   A      - the n-by-n state dynamics matrix A.
C   B      - the n-by-m input/state matrix B.
C   C      - the p-by-n state/output matrix C.
C   D      - (optional) if jobd = 1, the p-by-m input/output matrix D.
C            If jobd = 0, this parameter must not be given.
C   G      - the m-by-mv system input scaling matrix G.
C   H      - the pz-by-p system output scaling matrix H.
C   F      - (optional) if fbtype = 2, the m-by-p output feedback
C            matrix F. If fbtype = 1, then the feedback matrix is
C            assumed to be an m-by-m identity matrix. In this case,
C            or if alpha = 0, the parameter F must not be given.
C   K      - (optional) if beta <> 0, the m-by-n state feedback
C            matrix K. If beta = 0, the parameter K must not be given.
C
C Output parameters:
C   Ao     - the n-by-n state dynamics matrix Ac of the closed-loop
C            system.
C   Bo     - the n-by-mb input/state matrix Bc of the closed-loop
C            system, where mb = mv if task = 1, and mb = m, is task = 2.
C   Co     - the pc-by-n output/state matrix Cc of the closed-loop
C            system, where pc = pz if task = 1, and pc = p, is task = 2.
C   Do     - (optional) if jobd = 1, the pc-by-mb input/output
C            matrix Dc of the closed-loop system.
C   rcondi - (optional) the reciprocal condition number of the matrix
C            I - alpha*D*F.
C
C Method:
C
C   task = 1:  The matrices of the closed-loop system have the
C              expressions:
C
C     Ac = AC + beta*BC*K,      Bc = BC*G,
C     Cc = H*(CC + beta*DC*K),  Dc = H*DC*G,
C
C     where
C
C     AC = A + alpha*B*F*E*C,   BC = B + alpha*B*F*E*D,
C     CC = E*C,                 DC = E*D,
C
C     with E = (I - alpha*D*F)**-1.
C
C   task = 2:  The matrices of the closed-loop system have the
C              expressions:
C
C     Ac = A + alpha*B*F*E*C,  Bc = B + alpha*B*F*E*D,
C     Cc = E*C,                Dc = E*D,
C
C     where E = (I - alpha*D*F)**-1.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, July 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Feb. 2005,
C   Apr. 2009, Dec. 2010.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
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
      CHARACTER         FBTYPE, JOBD
      INTEGER           INFO, LDA, LDB, LDBC, LDC, LDCC, LDD, LDDC,
     $                  LDF, LDG, LDH, LDK, LDWORK, M, MV, N, P, PZ
      DOUBLE PRECISION  ALPHA, BETA, RCOND
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      mwPointer         A, B, BC, C, CC, D, DC, DWORK, F, G, H, K
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           TASK1, UNITF, WITHD
      INTEGER           IA, IB, IC, ID, IG, IH, IK, IP, ITMP, JF,
     $                  LIWORK, NC, ND, NDC, NF, NK
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          AB05RD, AB05SD
C .. Intrinsic functions ..
      INTRINSIC         MAX
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.7 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSFCONN requires at least 7 input arguments' )
      ELSE IF ( NLHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSFCONN requires at least 3 output arguments' )
      END IF
C
C   task, jobd, fbtype, alpha(, beta), A(nxn), B(n,m), C(p,n)(, D(p,m)),
C   G(m,mv), H(pz,p))(, F(m,p))(, K(m,n))
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR.
     $     mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'TASK must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'TASK must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      ITMP = TEMP
      IF ( ITMP.LT.1 .OR. ITMP.GT.2 ) THEN
         CALL mexErrMsgTxt
     $     ( 'TASK has 1 or 2 the only admissible values' )
      END IF
C
      TASK1 = ITMP.EQ.1
C
      IF ( mxGetM( PRHS(2) ).NE.1 .OR.
     $     mxGetN( PRHS(2) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'JOBD must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'JOBD must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), TEMP, 1 )
      ITMP = TEMP
      IF ( ITMP.LT.0 .OR. ITMP.GT.1 ) THEN
         CALL mexErrMsgTxt
     $     ( 'JOBD has 0 or 1 the only admissible values' )
      END IF
C
      WITHD = ITMP.EQ.1
C
      IF ( mxGetM( PRHS(3) ).NE.1 .OR.
     $     mxGetN( PRHS(3) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'FBTYPE must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'FBTYPE must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TEMP, 1 )
      ITMP = TEMP
      IF ( ITMP.LT.1 .OR. ITMP.GT.2 ) THEN
         CALL mexErrMsgTxt
     $     ( 'FBTYPE has 1 or 2 the only admissible values' )
      END IF
C
      UNITF = ITMP.EQ.1
C
      IF ( mxGetM( PRHS(4) ).NE.1 .OR.
     $     mxGetN( PRHS(4) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'ALPHA must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'ALPHA must be a real scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), ALPHA, 1 )
C
      IP = 5
      IF ( TASK1 ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'BETA must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'BETA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), BETA, 1 )
         IP = IP + 1
      END IF
C
      IA = IP
      IB = IA + 1
      IC = IB + 1
      IF ( WITHD ) THEN
         ID = IC + 1
      ELSE
         ID = IC
      END IF
      IF ( TASK1 ) THEN
         IG = ID + 1
         IH = IG + 1
         IF ( UNITF .OR. ALPHA.EQ.ZERO ) THEN
            JF = IH
         ELSE
            JF = IH + 1
         END IF
         IK = JF + 1
      ELSE
         JF = ID + 1
      END IF
C
      N = mxGetM( PRHS(IA) )
      M = mxGetN( PRHS(IB) )
      P = mxGetM( PRHS(IC) )
C
      IF ( mxGetN( PRHS(IA) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(IA) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(IA) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(IB) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same rowsize as A' )
      END IF
      IF ( mxIsNumeric( PRHS(IB) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(IB) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
C
      IF ( mxGetN( PRHS(IC) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'C must have the same columnsize as A')
      END IF
      IF ( mxIsNumeric( PRHS(IC) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(IC) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
C
      IF ( WITHD ) THEN
         IF ( mxGetM( PRHS(ID) ).NE.P ) THEN
            CALL mexErrMsgTxt( 'D must have the same rowsize as C' )
         END IF
         IF ( mxGetN( PRHS(ID) ).NE.M ) THEN
            CALL mexErrMsgTxt( 'D must have the same columnsize as B')
         END IF
         IF ( mxIsNumeric( PRHS(ID) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(ID) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'D must be a real matrix' )
         END IF
      END IF
C
      IF ( TASK1 ) THEN
         IF ( NRHS.LT.10 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SYSFCONN requires at least 10 input arguments' )
         END IF
C
         MV = mxGetN( PRHS(IG) )
         PZ = mxGetM( PRHS(IH) )
         IF ( mxGetM( PRHS(IG) ).NE.M ) THEN
            CALL mexErrMsgTxt(
     $           'The rowsize of G must equal the columnsize of B' )
         END IF
         IF ( mxIsNumeric( PRHS(IG) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IG) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'G must be a real matrix' )
         END IF
C
         IF ( mxGetN( PRHS(IH) ).NE.P ) THEN
            CALL mexErrMsgTxt(
     $           'The columnsize of H must equal the rowsize of C' )
         END IF
         IF ( mxIsNumeric( PRHS(IH) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IH) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'H must be a real matrix' )
         END IF
C
         IF ( BETA.NE.ZERO ) THEN
            IF ( mxGetM( PRHS(IK) ).NE.M ) THEN
               CALL mexErrMsgTxt(
     $              'The rowsize of K must equal the columnsize of B' )
            END IF
            IF ( mxGetN( PRHS(IK) ).NE.N ) THEN
               CALL mexErrMsgTxt(
     $           'K must have the same columnsize as A' )
            END IF
            IF ( mxIsNumeric( PRHS(IK) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IK) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'K must be a real matrix' )
            END IF
         END IF
      END IF
C
      IF ( .NOT.UNITF .AND. ALPHA.NE.ZERO ) THEN
         IF ( mxGetM( PRHS(JF) ).NE.M ) THEN
            CALL mexErrMsgTxt(
     $           'The rowsize of F must equal the columnsize of B' )
         END IF
         IF ( mxGetN( PRHS(JF) ).NE.P ) THEN
            CALL mexErrMsgTxt(
     $           'The columnsize of F must equal the rowsize of C' )
         END IF
         IF ( mxIsNumeric( PRHS(JF) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(JF) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'F must be a real matrix' )
         END IF
      END IF
C
      IF ( WITHD ) THEN
         JOBD = 'D'
      ELSE
         JOBD = 'Z'
      END IF
C
      IF ( UNITF ) THEN
         FBTYPE = 'I'
      ELSE
         FBTYPE = 'O'
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA = MAX( 1, N )
      LDB = LDA
      IF ( N.GT.0 ) THEN
         LDC = MAX( 1, P )
         NC  = N
      ELSE
         LDC = 1
         NC  = 1
      END IF
      IF ( WITHD ) THEN
         LDD = MAX( 1, P )
         ND  = M
      ELSE
         LDD = 1
         ND  = 1
      END IF
      IF ( UNITF .OR. ALPHA.EQ.ZERO ) THEN
         LDF = 1
         NF  = 1
      ELSE
         LDF = MAX( 1, M )
         NF  = P
      END IF
      IF ( TASK1 ) THEN
         LDG = MAX( 1, M )
         LDH = MAX( 1, PZ )
         IF ( BETA.NE.ZERO ) THEN
            LDK = MAX( 1, M )
            NK  = N
         ELSE
            LDK = 1
            NK  = 1
         END IF
         LDBC = LDA
         LDCC = LDH
         IF ( WITHD ) THEN
            LDDC = LDH
            NDC  = MV
         ELSE
            LDDC = 1
            NDC  = 1
         END IF
      END IF
C
C   ldwork
C   Use larger values, ldwork >= MAX( LDWORK, N*M, N*P ), for best
C   performance.
C
      IF ( WITHD ) THEN
         LIWORK = MAX( 1, 2*P )
         LDWORK = MAX( M, P*P + 4*P )
         IF ( TASK1 )
     $      LDWORK = MAX( LDWORK, P*MV )
      ELSE
         LIWORK = 1
         LDWORK = M
      END IF
      LDWORK = MAX( 1, LDWORK )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( IWORK( LIWORK ) )
      A     = mxCalloc(  LDA*N,  8 )
      B     = mxCalloc(  LDB*M,  8 )
      C     = mxCalloc(  LDC*NC, 8 )
      D     = mxCalloc(  LDD*ND, 8 )
      DWORK = mxCalloc( LDWORK,  8 )
      F     = mxCalloc( LDF*NF,  8 )
      IF ( TASK1 ) THEN
         BC = mxCalloc( LDBC*MV,  8 )
         CC = mxCalloc( LDCC*N,   8 )
         DC = mxCalloc( LDDC*NDC, 8 )
         G  = mxCalloc(  LDG*MV,  8 )
         H  = mxCalloc(  LDH*P,   8 )
         K  = mxCalloc(  LDK*NK,  8 )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(IA) ), %VAL( A ), N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(IB) ), %VAL( B ), N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(IC) ), %VAL( C ), P*N )
      IF ( WITHD )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(ID) ), %VAL( D ), P*M )
      IF ( TASK1 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IG) ), %VAL( G ),
     $                          M*MV )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IH) ), %VAL( H ),
     $                          PZ*P )
         IF ( BETA.NE.ZERO )
     $      CALL mxCopyPtrToReal8( mxGetPr( PRHS(IK) ), %VAL( K ),
     $                             M*N )
      END IF
      IF ( .NOT.UNITF .AND. ALPHA.NE.ZERO )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(JF) ), %VAL( F ),
     $                          M*P )
C
C Do the actual computations.
C
      IF ( TASK1 ) THEN
         CALL AB05RD( FBTYPE, JOBD, N, M, P, MV, PZ, ALPHA, BETA,
     $                %VAL( A ), LDA, %VAL( B ), LDB, %VAL( C ), LDC,
     $                %VAL( D ), LDD, %VAL( F ), LDF, %VAL( K ), LDK,
     $                %VAL( G ), LDG, %VAL( H ), LDH, RCOND, %VAL( BC ),
     $                LDBC, %VAL( CC ), LDCC, %VAL( DC ), LDDC, IWORK,
     $                %VAL( DWORK ), LDWORK, INFO )
      ELSE
         CALL AB05SD( FBTYPE, JOBD, N, M, P, ALPHA, %VAL( A ), LDA,
     $                %VAL( B ), LDB, %VAL( C ), LDC, %VAL( D ), LDD,
     $                %VAL( F ), LDF, RCOND, IWORK, %VAL( DWORK ),
     $                LDWORK, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      IP = 4
      IF ( INFO.EQ.0 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
         CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS(1) ), N*N )
         IF ( TASK1 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( N, MV, 0 )
            PLHS(3) = mxCreateDoubleMatrix( PZ, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( BC ), mxGetPr( PLHS(2) ),
     $                             N*MV )
            CALL mxCopyReal8ToPtr( %VAL( CC ), mxGetPr( PLHS(3) ),
     $                             PZ*N )
            IF ( NLHS.GE.IP .AND. WITHD ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( PZ, MV, 0 )
               CALL mxCopyReal8ToPtr( %VAL( DC ),
     $                                mxGetPr( PLHS(IP) ), PZ*MV )
               IP = IP + 1
            END IF
         ELSE
            PLHS(2) = mxCreateDoubleMatrix( N, M, 0 )
            PLHS(3) = mxCreateDoubleMatrix( P, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( B ), mxGetPr( PLHS(2) ),
     $                             N*M )
            CALL mxCopyReal8ToPtr( %VAL( C ), mxGetPr( PLHS(3) ),
     $                             P*N )
            IF ( NLHS.GE.IP .AND. WITHD ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( P, M, 0 )
               CALL mxCopyReal8ToPtr( %VAL( D ),
     $                                mxGetPr( PLHS(IP) ), P*M )
               IP = IP + 1
            END IF
         END IF
         IF ( NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( RCOND, mxGetPr( PLHS(IP) ), 1 )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( IWORK )
      CALL mxFree( A     )
      CALL mxFree( B     )
      CALL mxFree( C     )
      CALL mxFree( D     )
      CALL mxFree( DWORK )
      CALL mxFree( F     )
      IF ( TASK1 ) THEN
         CALL mxFree( BC )
         CALL mxFree( CC )
         CALL mxFree( DC )
         CALL mxFree( G )
         CALL mxFree( H )
         CALL mxFree( K )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( TASK1 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB05RD'')' )
     $         INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB05SD'')' )
     $         INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SYSFCONN ***
      END
