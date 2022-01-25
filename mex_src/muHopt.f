C MUHOPT.f   - Gateway function to compute the mu optimal or H_inf
C              controller using SLICOT routines SB10AD, SB10MD, AB04MD,
C              AB05MD, and AB07ND.
C
C RELEASE 2.0 of SLICOT Robust Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2004-2020 NICONET e.V.
C
C Matlab call:
C   If mu optimal controller is desired (job > 0):
C
C   [AK,BK,CK,DK(,mju,RCOND)] = muHopt(job,discr,A,B,C,D,ncon,nmeas,
C                gamma,omega,nblock,itype,ord(,qutol(,gtol(,actol))))
C
C   If H_inf optimal controller is only desired (job <= 0):
C
C   [AK,BK,CK,DK(,gammin,RCOND)] = muHopt(job,discr,A,B,C,D,ncon,
C                nmeas,gamma(,gtol(,actol)))
C
C Purpose:
C   To compute the matrices of the mu optimal or H_inf optimal
C   controller given the model in a state space. It also outputs the
C   mu norm of the closed loop system, if mu optimal controller is
C   desired, or the value of gamma reached in the H_inf synthesis
C   problem, if H_inf controller is only desired.
C   The discrete-time systems are handled via bilinear transformation
C   to continuous-time and then the controller obtained is discretised.
C   For the K step the SB10AD subroutine is employed, and the SB10MD
C   subroutine performs the D step.
C
C Input parameters:
C   job    - indicates the type of the controller as well as the strategy
C            for reducing the gamma value:
C            >  0 : mu optimal controller is desired;
C            <= 0 : H_inf optimal controller only is desired.
C            Specifically, job
C            = 0 : find suboptimal controller only;
C            and abs(job) specifies the strategy for reducing gamma:
C            = 1 : use bisection method for decreasing gamma from gamma
C                  to gammin until the closed-loop system leaves
C                  stability;
C            = 2 : scan from gamma to 0 trying to find the minimal gamma
C                  for which the closed-loop system retains stability;
C            = 3 : first bisection, then scanning.
C   discr  - indicates the type of the system, as follows:
C            = 0 : continuous-time system;
C            = 1 : discrete-time system.
C   A      - the n-by-n system state matrix A of the plant.
C   B      - the n-by-m system input matrix B of the plant.
C   C      - the p-by-n system output matrix C of the plant.
C   D      - the p-by-m system input/output matrix D of the plant.
C   ncon   - the number of control inputs.
C            p-nmeas >= ncon >= 0.
C   nmeas  - the number of measurements.
C            p-nmeas = m-ncon >= nmeas >= 0.
C   gamma  - the initial value of gamma on input. It is assumed that
C            gamma is sufficiently large so that the controller is
C            admissible.  gamma >= 0.
C   omega  - the vector of length lendat >= 2 with the frequencies.
C            They must be nonnegative, in increasing order, and
C            for discrete-time systems between 0 and pi.
C   nblock - the vector with the block structure of the uncertainty.
C            nblock(I) is the size of each block.
C   itype  - the vector of the same length as nblock indicating
C            the type of each block.
C            itype(I) = 1 indicates that the corresponding block is a
C            real block. THIS OPTION IS NOT SUPPORTED NOW.
C            itype(I) = 2 indicates that the corresponding block is a
C            complex block. THIS IS THE ONLY ALLOWED VALUE NOW!
C   ord    - the maximum order of each block in the D-fitting procedure.
C            1 <= ord <= lendat-1.
C   qutol  - (optional) the acceptable mean relative error between
C            the D(jw) and the frequency response of the estimated block
C            [ADi,BDi;CDi,DDi]. When it is reached, the result is
C            taken as good enough.
C            Default: qutol = 2.
C   gtol   - (optional) tolerance used for controlling the accuracy
C            of gamma and its distance to the estimated minimal possible
C            value of gamma.
C            If gtol <= 0, then sqrt(EPS) is used, where EPS is the
C            relative machine precision.
C            Default: gtol = 0.01.
C   actol  - (optional) upper bound for the poles of the closed-loop
C            system used for determining if it is stable.
C            actol <= 0 for stable systems.
C            Default: actol = 0.
C
C Output parameters:
C   AK     - the n-by-n controller state matrix AK.
C   BK     - the n-by-nmeas controller input matrix BK.
C   CK     - the ncon-by-n controller output matrix CK.
C   DK     - the ncon-by-nmeas controller input/output matrix DK.
C   mju    - (optional) the vector with the estimated upper bound of
C            the structured singular value for each frequency in omega
C            for the closed-loop system.
C   gammin - (optional) the estimated minimal admissible gamma.
C   RCOND  - (optional) for each successful J-th K step:
C            RCOND(J) contains the reciprocal condition number of the
C                     control transformation matrix;
C            RCOND(J+1) contains the reciprocal condition number of the
C                     measurement transformation matrix;
C            RCOND(J+2) contains an estimate of the reciprocal condition
C                     number of the X-Riccati equation;
C            RCOND(J+3) contains an estimate of the reciprocal condition
C                     number of the Y-Riccati equation.
C            If job = 0, only RCOND(1:4) are set by the first K step.
C
C Contributor:
C   A. Markovski, Technical University of Sofia, October 2003.
C
C Revisions:
C   V. Sima, April 2004, Sept. 2004, Mar. 2005, Apr. 2009, Dec. 2012.
C
C***********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      INTEGER           MAXIT, HNPTS
      PARAMETER         ( MAXIT = 15, HNPTS = 2048 )
      DOUBLE PRECISION  ZERO, ONE, TWO, P01
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                    P01  = 0.01D0 )
C
C .. Mex-file interface parameters ..
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Parameters used by SLICOT subroutines ..
      INTEGER           F, INFO, JOB, LBWORK, LDA, LDAC, LDAD, LDAE,
     $                  LDAK, LDB, LDBC, LDBD, LDBE, LDBK, LDC, LDCC,
     $                  LDCD, LDCE, LDCK, LDD, LDDC, LDDD, LDDE, LDDK,
     $                  LDWORK, LENDAT, LIWORK, LZWORK, M, M2, MNB, N,
     $                  N2E, NE, NEB, NP, NP1, NP2, NTEMP, ORD, TOTORD
      DOUBLE PRECISION  ACTOL, GAMMA, GTOL, QUTOL, TEMP
      DOUBLE PRECISION  RCOND(4*MAXIT)
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: BWORK(:)
      INTEGER,          ALLOCATABLE :: ITYPE(:), IWORK(:), NBLOCK(:)
      DOUBLE PRECISION, ALLOCATABLE :: A(:), AC(:), AD(:), AE(:), AK(:),
     $                                 AKB(:), B(:), BC(:), BD(:),
     $                                 BE(:), BK(:), BKB(:), C(:),
     $                                 CC(:), CD(:), CE(:), CK(:),
     $                                 CKB(:), D(:), DC(:), DD(:),
     $                                 DE(:), DK(:), DKB(:), DWORK(:),
     $                                 MJU(:), OMEGA(:), PMJU(:),
     $                                 RITYPE(:), RNBLCK(:)
      COMPLEX*16,       ALLOCATABLE :: ZWORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER         CONJOB
      CHARACTER*120     TEXT
      INTEGER           DISCR, I, IP, ITER, ITERB, LD1, LD2, LD3, LD4,
     $                  LI1, LI2, LI3, LI4, LW1, LW2, LW3, LW4, LW5,
     $                  LW6, LW7, LWA, LWB, LZD, LZM, M1, M11, MD, MIT,
     $                  MN, NA, NB, NC, ND, NIT, NNB, NP11
      DOUBLE PRECISION  PMPEAK, MUPEAK, TOL
C
C .. External functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C
C .. External subroutines ..
      EXTERNAL          AB04MD, AB05MD, AB07ND, DCOPY, SB10AD, SB10MD
C
C ..Intrinsic functions..
      INTRINSIC         ABS, COS, INT, MAX, MIN, SQRT
C
C Check for proper number of arguments.
C
      IF ( NRHS.EQ.0 ) THEN
         CALL mexErrMsgTxt( 'Matlab call: [AK,BK,CK,DK,MJU,RCOND]='//
     $ 'MUHOPT(JOB,DISCR,A,B,C,D,NCON,NMEAS,GAMMA,OMEGA,NBLOCK,ITYPE,'//
     $ 'ORD[,QUTOL[,GTOL[,ACTOL]]]) or [AK,BK,CK,DK,GAMMIN,RCOND]='//
     $ 'MUHOPT(JOB,DISCR,A,B,C,D,NCON,NMEAS,GAMMA[,GTOL[,ACTOL]])')
      ELSE IF ( NRHS.LT.9 ) THEN
         CALL mexErrMsgTxt
     $                 ( 'MUHOPT requires at least 9 input arguments' )
      ELSE IF ( NRHS.GT.16 ) THEN
         CALL mexErrMsgTxt
     $                 ( 'MUHOPT requires at most 16 input arguments' )
      ELSE IF ( NLHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $                 ( 'MUHOPT requires at least 4 output arguments' )
      END IF
C
C  job
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be a scalar' )
      ELSE IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $          mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
      END IF
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      JOB = INT( TEMP )
C
      IF ( ABS( JOB ).GT.3 ) THEN
         CALL mexErrMsgTxt( 'JOB has -3, -2, -1, 0, 1, 2, or 3 '//
     $                      'the only admissible values' )
      END IF
C
C     Determine the job.
C
      IF ( JOB.GT.0 ) THEN
C        mu controller desired.
         CONJOB = 'M'
      ELSE
C        H_inf controller desired.
         CONJOB = 'H'
         IF ( JOB.EQ.0 ) THEN
            JOB = 4
         ELSE
            JOB = ABS( JOB )
         END IF
      END IF
C
C     Recheck for proper number of arguments.
C
      IF ( CONJOB.EQ.'M' .AND. NRHS.LT.13 ) THEN
         CALL mexErrMsgTxt( 'MUHOPT requires at least 13 input '//
     $                  'arguments if mu optimal controller is desired')
      ELSE IF ( CONJOB.EQ.'H' .AND. NRHS.GT.11 ) THEN
         CALL mexErrMsgTxt( 'MUHOPT requires at most 11 input '//
     $               'arguments if H_inf optimal controller is desired')
      END IF
C
C  discr
C
      IF ( mxGetM( PRHS(2) ).NE.1 .OR. mxGetN( PRHS(2) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'DISCR must be a scalar' )
      ELSE IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $          mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'DISCR must be an integer scalar' )
      END IF
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), TEMP, 1 )
      DISCR = INT( TEMP )
C
      IF ( DISCR.LT.0 .OR. DISCR.GT.1 ) THEN
         CALL mexErrMsgTxt( 'DISCR must be 0 or 1' )
      END IF
C
C  A, B, C, D
C
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a numeric matrix' )
      END IF
C
      N  = mxGetM( PRHS(3) )
      NA = mxGetN( PRHS(3) )
C
      IF ( NA.NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
C
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a numeric matrix' )
      END IF
C
      NB = mxGetM( PRHS(4) )
      M  = mxGetN( PRHS(4) )
C
      IF ( NB.NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same row dimension as A' )
      END IF
C
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a numeric matrix' )
      END IF
C
      NP = mxGetM( PRHS(5) )
      NC = mxGetN( PRHS(5) )
C
      IF ( NC.NE.N ) THEN
         CALL mexErrMsgTxt
     $                 ( 'C must have the same column dimension as A' )
      END IF
C
      IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(6) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a numeric matrix' )
      END IF
C
      ND = mxGetM( PRHS(6) )
      MD = mxGetN( PRHS(6) )
C
      IF ( ND.NE.NP ) THEN
         CALL mexErrMsgTxt( 'D must have the same row dimension as C' )
      ELSE IF ( MD.NE.M ) THEN
         CALL mexErrMsgTxt
     $                 ( 'D must have the same column dimension as B' )
      END IF
C
C  ncon (M2), nmeas (NP2)
C
      IF ( mxGetM( PRHS(7) ).NE.1 .OR. mxGetN( PRHS(7) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'NCON must be a scalar' )
      ELSE IF ( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $          mxIsComplex( PRHS(7) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'NCON must be an integer scalar' )
      END IF
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), TEMP, 1 )
      M2 = INT( TEMP )
C
      IF ( mxGetM( PRHS(8) ).NE.1 .OR. mxGetN( PRHS(8) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'NMEAS must be a scalar' )
      ELSE IF ( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $          mxIsComplex( PRHS(8) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'NMEAS must be an integer scalar' )
      END IF
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(8) ), TEMP, 1 )
      NP2 = INT( TEMP )
      M1  = M  - M2
      NP1 = NP - NP2
C
      IF ( M1.NE.NP1 ) THEN
         CALL mexErrMsgTxt( 'M - NCON must be equal to P - NMEAS' )
      END IF
C
C  gamma
C
      IF ( mxGetM( PRHS(9) ).NE.1 .OR. mxGetN( PRHS(9) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'GAMMA must be a scalar' )
      ELSE IF ( mxIsNumeric( PRHS(9) ).EQ.0 .OR.
     $          mxIsComplex( PRHS(9) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'GAMMA must be a real scalar' )
      END IF
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(9) ), GAMMA, 1 )
C
      IF ( CONJOB.EQ.'M' ) THEN
C
C     The mu controller case.
C
C  omega
C
         IF ( mxIsNumeric( PRHS(10) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(10) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'OMEGA must be a numeric vector' )
         ELSE IF ( MIN( mxGetM( PRHS(10) ), mxGetN( PRHS(10) ) ).NE.1 )
     $        THEN
            CALL mexErrMsgTxt( 'OMEGA must be a vector' )
         END IF
C
         LENDAT = mxGetM( PRHS(10) )*mxGetN( PRHS(10) )
         IF ( LENDAT.LE.1 ) THEN
            CALL mexErrMsgTxt( 'OMEGA must have at least 2 elements' )
         END IF
C
C  nblock
C
         IF ( mxIsNumeric( PRHS(11) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(11) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'NBLOCK must be a numeric vector' )
         END IF
C
         MNB = mxGetM( PRHS(11) )
         NNB = mxGetN( PRHS(11) )
         MN  = MIN( MNB, NNB )
         MNB = MAX( MNB, NNB )
         NNB = MN
         IF ( NNB.NE.1 ) THEN
            CALL mexErrMsgTxt( 'NBLOCK must be a vector' )
         END IF
C
C  itype
C
         IF ( mxIsNumeric( PRHS(12) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(12) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ITYPE must be a numeric vector' )
         END IF
C
         MIT = mxGetM( PRHS(12) )
         NIT = mxGetN( PRHS(12) )
         MN  = MIN( MIT, NIT )
         MIT = MAX( MIT, NIT )
         NIT = MN
C
         IF ( NIT.NE.1 ) THEN
            CALL mexErrMsgTxt( 'ITYPE must be a vector' )
         ELSE IF ( MNB.NE.MIT ) THEN
            CALL mexErrMsgTxt
     $                     ( 'ITYPE must have the same size as NBLOCK' )
         END IF
C
C  ord
C
         IF ( mxGetM( PRHS(13) ).NE.1 .OR.
     $        mxGetN( PRHS(13) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ORD must be a scalar' )
         ELSE IF ( mxIsNumeric( PRHS(13) ).EQ.0 .OR.
     $             mxIsComplex( PRHS(13) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ORD must be an integer scalar' )
         END IF
C
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(13) ), TEMP, 1 )
         ORD = INT( TEMP )
C
         IF ( ORD.LT.1 ) THEN
            CALL mexErrMsgTxt( 'ORD must be at least 1' )
         ELSE IF ( ORD.GE.LENDAT ) THEN
            WRITE( TEXT, '(''ORD must be less than LENDAT - 1 = '',
     $             I7)' ) LENDAT - 1
            CALL mexErrMsgTxt( TEXT )
         END IF
C
C  qutol
C
         IF ( NRHS.GE.14 ) THEN
            IF ( mxGetM( PRHS(14) ).NE.1 .OR.
     $           mxGetN( PRHS(14) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'QUTOL must be a scalar' )
            ELSE IF ( mxIsNumeric( PRHS(14) ).EQ.0 .OR.
     $                mxIsComplex( PRHS(14) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'QUTOL must be a real scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(14) ), QUTOL, 1 )
         ELSE
            QUTOL = TWO
         END IF
         IP = 15
      ELSE
         IP = 10
      END IF
C
C  gtol
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'GTOL must be a scalar' )
         ELSE IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $             mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'GTOL must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), GTOL, 1 )
         IP = IP + 1
      ELSE
         GTOL = P01
      END IF
C
C  actol
C
      IF ( NRHS.EQ.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ACTOL must be a scalar' )
         ELSE IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $             mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ACTOL must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr(PRHS(IP) ), ACTOL, 1 )
      ELSE
         ACTOL = ZERO
      END IF
C
C     Set the default tolerance.
C
      TOL = SQRT( DLAMCH( 'Epsilon' ) )
C
C Determine the lenghts of working arrays.
C
C     The original system.
C
      LDA = MAX( 1, N )
      LDB = LDA
      LDC = MAX( 1, NP )
      LDD = LDC
C
      IF ( CONJOB.EQ.'M' ) THEN
C
C        The scaling system.
C
         TOTORD = NP1*ORD
C
         F    = MAX( M2, NP2 )
         LDAD = MAX( 1, TOTORD )
         LDBD = LDAD
         LDCD = MAX( 1, NP1 + F )
         LDDD = LDCD
C
C        The extended system.
C
         NE   = N + 2*TOTORD
         LDAE = MAX( 1, NE )
         LDBE = LDAE
         LDCE = LDCD
         LDDE = LDCE
      ELSE
         NE   = N
         LDAE = LDA
      END IF
C
C     The closed-loop system.
C
      N2E  = 2*NE
      LDAC = MAX( 1, N2E )
      LDBC = LDAC
      LDCC = MAX( 1, NP1 )
      LDDC = LDCC
C
C     The controller.
C
      LDAK = LDAE
      LDBK = LDAK
      LDCK = MAX( 1, M2 )
      LDDK = LDCK
C
C     LBWORK
C
C     ..SB10AD..
      LBWORK = N2E
C
C     LIWORK
C
C     ..SB10AD..
      LI1 = MAX( 2*MAX( NE, M - M2, NP - NP2, M2, NP2 ), NE*NE )
C
      IF ( CONJOB.EQ.'M' ) THEN
C        ..SB10MD..
         LI2 = MAX( N2E, 4*MNB - 2, NP1 )
         IF ( QUTOL.GE.ZERO )
     $      LI2 = MAX( LI2, 2*ORD + 1  )
C        ..AB07ND..
         LI3 = 2*M
      ELSE
         LI2 = 0
         LI3 = 0
      END IF
C
C     ..AB04MD..
      LI4 = NE
C
      LIWORK = MAX( LI1, LI2, LI3, LI4 )
C
C     LDWORK
C
C     ..AB04MD..
      LD1 = LDAK
C
C     ..SB10AD..
      NP11 = NP1 - M2
      M11  = M1 - NP2
      LW1 = NE*M + NP*NE + NP*M + M2*M2 + NP2*NP2
      LW2 = MAX( ( NE + NP1 + 1 )*( NE + M2 ) +
     $             MAX( 3*( NE + M2 ) + NE + NP1, 5*( NE + M2 ) ),
     $           ( NE + NP2 )*( NE + M1 + 1 ) +
     $             MAX( 3*( NE + NP2 ) + NE + M1, 5*( NE + NP2 ) ),
     $           M2 + NP1*NP1 + MAX( NP1*MAX( NE, M1 ), 3*M2 + NP1,
     $                               5*M2 ),
     $           NP2 + M1*M1 +  MAX( MAX( NE, NP1 )*M1, 3*NP2 + M1,
     $                               5*NP2 ) )
      LW3 = MAX( NP11*M1 + MAX( 4*MIN( NP11, M1 ) + MAX( NP11, M1 ),
     $                          6*MIN( NP11, M1 ) ),
     $           NP1*M11 + MAX( 4*MIN( NP1, M11 ) + MAX( NP1, M11 ),
     $                          6*MIN( NP1, M11 ) ) )
      LW4 = 2*M*M + NP*NP + 2*M*NE + M*NP + 2*NE*NP
      LW5 = 2*NE*NE + M*NE + NE*NP
      LW6 = MAX( M*M + MAX( 2*M1, 3*NE*NE +
     $                      MAX( NE*M, 10*NE*NE + 12*NE + 5 ) ),
     $           NP*NP + MAX( 2*NP1, 3*NE*NE +
     $                        MAX( NE*NP, 10*NE*NE + 12*NE + 5 ) ) )
      LW7 = M2*NP2 + NP2*NP2 + M2*M2 +
     $      MAX( NP11*NP11 + MAX( 2*NP11, ( NP11 + M11 )*NP2 ),
     $           M11*M11 + MAX( 2*M11, M11*M2 ), 3*NE,
     $           NE*( 2*NP2 + M2 ) +
     $           MAX( 2*NE*M2, M2*NP2 +
     $                         MAX( M2*M2 + 3*M2, NP2*( 2*NP2 +
     $                              M2 + MAX( NP2, NE ) ) ) ) )
      LD2 = LW1 + MAX( 1, LW2, LW3, LW4, LW5 + MAX( LW6, LW7 ) )
C
      IF ( CONJOB.EQ.'M' ) THEN
C
C        ..SB10MD..
         MN  = MIN( 2*LENDAT, 2*ORD + 1 )
         LWA = NP1*LENDAT + 2*MNB + NP1 - 1
         LWB = LENDAT*( NP1 + 2 ) + ORD*( ORD + 2 ) + 1
         LW1 = 2*LENDAT + 4*HNPTS
         LW2 =   LENDAT + 6*HNPTS
         LW3 = 2*LENDAT*( 2*ORD + 1 ) + MAX( 2*LENDAT, 2*ORD + 1 ) +
     $                                  MAX( MN + 6*ORD + 4, 2*MN + 1 )
         LW4 = MAX( ORD*ORD + 5*ORD, 6*ORD + 1 + MIN( 1, ORD ) )
C
         LD3 = LWA + MAX( N2E + MAX( N2E, NP1 - 1 ),
     $                    2*NP1*NP1*MNB - NP1*NP1 + 9*MNB*MNB +
     $                      NP1*MNB + 11*NP1 + 33*MNB - 11 )
         IF ( QUTOL.GE.ZERO ) THEN
            LD4 = LWB + MAX( LW1, LW2, LW3, LW4 )
         ELSE
            LD4 = 0
         END IF
C        .. AB05MD..
         LD4 = MAX( LD4, MAX( M, NP )*MAX( N + TOTORD, M, NP ) )
C        .. AB07ND..
         LD4 = MAX( LD4, 4*M )
      ELSE
         LD3 = 0
         LD4 = 0
      END IF
C
      LDWORK = MAX( LD1, LD2, LD3, LD4 )
C
C     LZWORK
C
      IF ( CONJOB.EQ.'M' ) THEN
C        ..SB10MD..
         LZM = MAX( NP1*NP1 + N2E*NP1 + N2E*N2E + 2*N2E,
     $              6*NP1*NP1*MNB + 13*NP1*NP1 + 6*MNB + 6*NP1 - 3 )
         IF ( QUTOL.GE.ZERO ) THEN
            LZD = MAX( LENDAT*( 2*ORD + 3 ), ORD*ORD + 3*ORD + 1 )
         ELSE
            LZD = 0
         END IF
         LZWORK = MAX( LZM, LZD )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( A(LDA*N), AC(LDAC*N2E), AK(LDAK*NE),
     $           B(LDB*M), BC(LDBC*M1),  BK(LDBK*NP2), BWORK(LBWORK),
     $           C(LDC*N), CC(LDCC*N2E), CK(LDCK*NE),
     $           D(LDD*M), DC(LDDC*M1),  DK(LDDK*NP2), DWORK(LDWORK),
     $           IWORK(LIWORK) )
C
      IF ( CONJOB.EQ.'M' ) THEN
         ALLOCATE ( AD(LDAD*TOTORD), AE(LDAE*NE),     AKB(LDAK*NE),
     $              BD(LDBD*(M1+F)), BE(LDBE*(M1+F)), BKB(LDBK*NP2),
     $              CD(LDCD*TOTORD), CE(LDCE*NE),     CKB(LDCK*NE),
     $              DD(LDDD*(M1+F)), DE(LDDE*(M1+F)), DKB(LDDK*NP2),
     $              ITYPE(MNB), MJU(LENDAT), NBLOCK(MNB), OMEGA(LENDAT),
     $              PMJU(LENDAT), RITYPE(MNB), RNBLCK(MNB),
     $              ZWORK(LZWORK) )
      END IF
C
C Copy right hand side arguments to local arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), C, NP*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), D, NP*M )
C
      IF ( CONJOB.EQ.'M' ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(10) ), OMEGA, LENDAT )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(11) ), RNBLCK, MNB )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(12) ), RITYPE, MNB )
C
         DO 10 I = 1, MNB
            NBLOCK(I) = INT( RNBLCK(I) )
            ITYPE(I)  = INT( RITYPE(I) )
   10    CONTINUE
      END IF
C
C Do the actual computations.
C
C     Transform to continuous-case, if needed.
C
      IF ( DISCR.EQ.1 ) THEN
         CALL AB04MD( 'D', N, M, NP, ONE, ONE, A, LDA, B, LDB, C, LDC,
     $                D, LDD, IWORK, DWORK, LDWORK, INFO )
C
         IF ( INFO.NE.0 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB04MD'')' )
     $             INFO
            GOTO 60
         END IF
C
         IF ( CONJOB.EQ.'M' ) THEN
C
            DO 20 I = 1, LENDAT
               OMEGA(I) = SQRT( ( ONE - COS( OMEGA(I) ) ) /
     $                          ( ONE + COS( OMEGA(I) ) + TOL ) )
   20       CONTINUE
C
         END IF
C
      END IF
C
C     Set parameters for the first K step.
C
      N2E  = 2*N
      NE   = N
      NEB  = N
      ITER = 1
C
C     First K step - makes use of the original system.
C
      CALL SB10AD( JOB, N, M, NP, M2, NP2, GAMMA, A, LDA, B, LDB, C,
     $             LDC, D, LDD, AK, LDA, BK, LDB, CK, LDCK, DK, LDDK,
     $             AC, LDAC, BC, LDBC, CC, LDCC, DC, LDDC, RCOND, GTOL,
     $             ACTOL, IWORK, LIWORK, DWORK, LDWORK, BWORK, LBWORK,
     $             INFO )
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB10AD'')' )
     $          INFO
         GOTO 60
      END IF
C
C     Skip the D step if H_inf controller only desired.
C
      IF ( CONJOB.EQ.'H' )
     $   GOTO 50
C
      PMPEAK = GAMMA + P01
C
C     Start the iteration process
C     --------- D step -------------------------------------------------
C
   30 CONTINUE
C
         CALL SB10MD( N2E, NP1, LENDAT, F, ORD, MNB, NBLOCK, ITYPE,
     $                QUTOL, AC, LDAC, BC, LDBC, CC, LDCC, DC, LDDC,
     $                OMEGA, TOTORD, AD, LDAD, BD, LDBD, CD, LDCD, DD,
     $                LDDD, MJU, IWORK, LIWORK, DWORK, LDWORK, ZWORK,
     $                LZWORK, INFO )
C
         IF ( INFO.NE.0 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB10MD'')' )
     $             INFO
            GOTO 60
         END IF
C
C        Check mu.
C
         MUPEAK = ZERO
C
         DO 40 I = 1, LENDAT
            MUPEAK = MAX( MUPEAK, MJU(I) )
   40    CONTINUE
C
         IF ( MUPEAK.GT.PMPEAK ) THEN
            IF ( ITER.NE.1 )
     $         CALL DCOPY( LENDAT, PMJU, 1, MJU, 1 )
            GOTO 50
         ELSE
C
C           Save the best controller.
C
            PMPEAK = MUPEAK
            ITERB  = ITER
            NEB    = NE
C
            CALL DCOPY( LENDAT, MJU, 1, PMJU, 1 )
            IF ( ITER.NE.1 ) THEN
               CALL DCOPY( NE*NE,  AKB, 1, AK, 1 )
               CALL DCOPY( NE*NP2, BKB, 1, BK, 1 )
               CALL DCOPY( M2*NE,  CKB, 1, CK, 1 )
               CALL DCOPY( M2*NP2, DKB, 1, DK, 1 )
            END IF
         END IF
C
         ITER = ITER + 1
         IF ( ITER.GT.MAXIT )
     $      GOTO 50
C
C        Dl*P.
C
         CALL AB05MD( 'U', 'N', N, M, NP, TOTORD, NP, A, LDA, B, LDB,
     $                C, LDC, D, LDD, AD, LDAD, BD, LDBD, CD, LDCD, DD,
     $                LDDD, NTEMP, AE, LDAE, BE, LDBE, CE, LDCE, DE,
     $                LDDE, DWORK, LDWORK, INFO )
C
C        inv(Dr).
C
         CALL AB07ND( TOTORD, M, AD, LDAD, BD, LDBD, CD, LDCD, DD, LDDD,
     $                TEMP, IWORK, DWORK, LDWORK, INFO )
C
         IF ( INFO.NE.0 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB07ND'')' )
     $             INFO
            GOTO 60
         END IF
C
C        Dl*P*inv(Dr).
C
         CALL AB05MD( 'U', 'O', TOTORD, M, M, NTEMP, NP, AD, LDAD, BD,
     $                LDBD, CD, LDCD, DD, LDDD, AE, LDAE, BE, LDBE,
     $                CE, LDCE, DE, LDDE, NE, AE, LDAE, BE, LDBE, CE,
     $                LDCE, DE, LDDE, DWORK, LDWORK, INFO )
         N2E  = 2*NE
         LDAK = MAX( 1, NE )
         LDBK = LDAK
C
C        --------- K step ----------------------------------------------
C
         CALL SB10AD( JOB, NE, M, NP, M2, NP2, GAMMA, AE, LDAE, BE,
     $                LDBE, CE, LDCE, DE, LDDE, AKB, LDAK, BKB, LDBK,
     $                CKB, LDCK, DKB, LDDK, AC, LDAC, BC, LDBC, CC,
     $                LDCC, DC, LDDC, RCOND( 4*ITER-3 ), GTOL, ACTOL,
     $                IWORK, LIWORK, DWORK, LDWORK, BWORK, LBWORK,
     $                INFO )
C
         IF ( INFO.NE.0 ) THEN
            INFO = 0
            GOTO 50
         END IF
C
      GOTO 30
C
C     ---- End of the mu synthesis -------------------------------------
C
   50 CONTINUE
C
C     Transform back to discrete time if needed.
C
      IF ( DISCR.EQ.1 ) THEN
         CALL AB04MD( 'C', NEB, NP2, M2, ONE, ONE, AK, LDAK, BK, LDBK,
     $                CK, LDCK, DK, LDDK, IWORK, DWORK, LDWORK, INFO )
C
         IF ( INFO.NE.0 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB04MD'')' )
     $             INFO
            GOTO 60
         END IF
      END IF
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( NEB, NEB, 0 )
      PLHS(2) = mxCreateDoubleMatrix( NEB, NP2, 0 )
      PLHS(3) = mxCreateDoubleMatrix( M2,  NEB, 0 )
      PLHS(4) = mxCreateDoubleMatrix( M2,  NP2, 0 )
      CALL mxCopyReal8ToPtr( AK, mxGetPr( PLHS(1) ), NEB*NEB )
      CALL mxCopyReal8ToPtr( BK, mxGetPr( PLHS(2) ), NEB*NP2 )
      CALL mxCopyReal8ToPtr( CK, mxGetPr( PLHS(3) ), M2*NEB )
      CALL mxCopyReal8ToPtr( DK, mxGetPr( PLHS(4) ), M2*NP2 )
C
      IF ( NLHS.GE.5 ) THEN
         IF ( CONJOB.EQ.'M' ) THEN
C           mu
            PLHS(5) = mxCreateDoubleMatrix( LENDAT, 1, 0 )
            CALL mxCopyReal8ToPtr( MJU, mxGetPr( PLHS(5) ), LENDAT )
            IF ( NLHS.GE.6 ) THEN
               PLHS(6) = mxCreateDoubleMatrix( ITERB*4, 1, 0 )
               CALL mxCopyReal8ToPtr( RCOND, mxGetPr( PLHS(6) ),
     $                                ITERB*4 )
            END IF
         ELSE
C           H_inf
            PLHS(5) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( GAMMA, mxGetPr( PLHS(5) ), 1 )
            IF ( NLHS.GE.6 ) THEN
               PLHS(6) = mxCreateDoubleMatrix( 4, 1, 0 )
               CALL mxCopyReal8ToPtr( RCOND, mxGetPr( PLHS(6) ), 4 )
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
   60 CONTINUE
C
      DEALLOCATE ( A, AC, AK, B, BC, BK, BWORK, C, CC, CK, D, DC, DK,
     $             DWORK, IWORK )
C
      IF ( CONJOB.EQ.'M' )
     $   DEALLOCATE ( AD, AE, AKB, BD, BE, BKB, CD, CE, CKB, DD, DE,
     $                DKB, ITYPE, MJU, NBLOCK, OMEGA, PMJU, RITYPE,
     $                RNBLCK, ZWORK )
C
C Error and warning handling ..
C
      IF ( INFO.NE.0 )
     $   CALL mexErrMsgTxt( TEXT )
C
      RETURN
C
C *** Last line of MUHOPT ***
      END
