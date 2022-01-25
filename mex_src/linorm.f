C LINORM.F - Gateway function for computation of the L_infinity norm of
C            a continuous-time or discrete-time system, in standard or
C            descriptor form, using SLICOT routine AB13DD.
C
C RELEASE 2.0 of SLICOT Robust Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [gpeak(,fpeak)] = linorm(A,E,B,C,D(,dico,systype,equil,fpeak0,tol))
C
C Purpose:
C   To compute the L-infinity norm of the continuous-time or
C   discrete-time system, either standard or in the descriptor form,
C
C                                     -1
C        G(lambda) = C*( lambda*E - A ) *B + D .
C
C Input parameters:
C   A       - the n-by-n system state matrix A.
C   E       - the n-by-n descriptor matrix E of the system, or an empty
C             matrix (i.e., with 0 rows and/or columns), in which case
C             E is taken as an identity matrix of order n (standard
C             system). The contents of E are ignored if it is non-empty
C             and systype = 1, but its order is then checked out.
C   B       - the n-by-m system input matrix B.
C   C       - the p-by-n system output matrix C.
C   D       - the p-by-m system matrix D.
C   dico    - (optional) specifies the type of the system:
C             = 1 : continuous-time system;
C             = 2 : discrete-time system.
C             Default: dico = 1.
C   systype - (optional) specifies whether or not the system is of
C             descriptor type:
C             = 0 : descriptor system;
C             = 1 : standard system (E = I).
C             Default: systype = 0.
C   equil   - (optional) specifies whether the user wishes to
C             preliminarily equilibrate the system (A,E,B,C) or (A,B,C):
C             = 1 : do not perform equilibration;
C             = 2 : perform equilibration (scaling).
C             Default: equil = 1.
C   fpeak0  - (optional) array of length 2 containing an estimate of the
C             frequency where the gain of the frequency response would
C             achieve its peak value. If fpeak0(2) = 0, the frequency is
C             infinite.
C             Default: fpeak0 = [0; 1].
C   tol     - (optional) tolerance used to set the accuracy in
C             determining the norm.
C             Default: sqrt(epsilon_machine) where epsilon_machine is
C             the relative machine precision.
C
C Output parameters:
C   gpeak   - array of length 2 containing the value of the L_infinity
C             norm of the system. If gpeak(2) = 0, the norm is infinite.
C   fpeak   - (optional) array of length 2 containing the frequency
C             for which the frequency response achieves its peak
C             value gpeak. If fpeak(2) = 0, the frequency is infinite.
C
C Contributor:
C   D. Sima, University of Bucharest, Romania, May 2001.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, June 2001,
C   May 2003, Apr. 2009, Dec. 2012, Apr. 2017.
C
C **********************************************************************
C
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO, HUNDRD
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0, HUNDRD = 1.0D2 )
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
      CHARACTER         DICO, EQUIL, JOBD, JOBE
      INTEGER           INFO, LCWORK, LDA, LDB, LDC, LDD, LDE, LDWORK,
     $                  M, N, P
      DOUBLE PRECISION  TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), D(:,:),
     $                                 DWORK(:), E(:,:)
      COMPLEX*16,       ALLOCATABLE :: CWORK(:)
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           NODYN, UNITE, USEPEN, WITHD
      CHARACTER*120     TEXT
      INTEGER           I, IA, IB, IC, ICOL, ID, IE, IP, IR, ISIZE,
     $                  ISYS, ISYSD, ITMP, J, LIWORK, LW, MINPM, MINWRK,
     $                  N2, N2PM, NN, PM
      DOUBLE PRECISION  BNORM, CNORM, TEMP
      DOUBLE PRECISION  FPEAK(2), FPK(2), GPEAK(2)
C
C .. External functions ..
      LOGICAL           LSAME, MA02HD
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, LSAME, MA02HD
C
C .. External subroutines ..
      EXTERNAL          AB13DD
C
C ..Intrinsic functions..
      INTRINSIC         MAX, MIN, SQRT
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'LINORM requires at least 5 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'LINORM requires at least 1 output argument' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   A(nxn), E(nxn), B(nxm), C(pxn), D(pxm)
C
      N = mxGetM( PRHS(1) )
      M = mxGetN( PRHS(3) )
      P = mxGetM( PRHS(4) )
C
      IF ( mxGetN( PRHS(1) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a numeric matrix' )
      END IF
      IF ( mxGetM( PRHS(2) )*mxGetN( PRHS(2) ).NE.0 ) THEN
         IF ( mxGetM( PRHS(2) ).NE.N .OR. mxGetN( PRHS(2) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $           ( 'E must be a square matrix of the same order as A' )
         END IF
         IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(2) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'E must be a numeric matrix' )
         END IF
         UNITE = .FALSE.
      ELSE
         UNITE = .TRUE.
      END IF
      IF ( mxGetM( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same row dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a numeric matrix' )
      END IF
      IF ( mxGetN( PRHS(4) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $       ('C must have the same column dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a numeric matrix' )
      END IF
      IF ( mxGetM( PRHS(5) ).NE.P ) THEN
         CALL mexErrMsgTxt( 'D must have the same row dimension as C' )
      END IF
      IF ( mxGetN( PRHS(5) ).NE.M ) THEN
         CALL mexErrMsgTxt
     $       ( 'D must have the same column dimension as B' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a numeric matrix' )
      END IF
C
C   dico, systype, equil, fpeak0, tol
C
      IP = 6
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'DICO must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'DICO must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ISYS = TEMP
         IF ( ISYS.NE.1 .AND. ISYS.NE.2 ) THEN
              CALL mexErrMsgTxt( 'DICO must be 1 or 2' )
         END IF
         IF ( ISYS.EQ.2 ) THEN
            DICO = 'D'
         ELSE
            DICO = 'C'
         END IF
         IP = IP + 1
      ELSE
         ISYS = 1
         DICO = 'C'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'SYSTYPE must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'SYSTYPE must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ISYSD = TEMP
         IF ( ISYSD.NE.0 .AND. ISYSD.NE.1 ) THEN
              CALL mexErrMsgTxt( 'SYSTYPE must be 0 or 1' )
         END IF
         IF ( ISYSD.EQ.1 .OR. UNITE ) THEN
            ISYSD = 1
            JOBE = 'I'
         ELSE
            JOBE = 'G'
         END IF
         IP = IP + 1
      ELSE
         IF ( UNITE ) THEN
            ISYSD = 1
            JOBE  = 'I'
         ELSE
            ISYSD = 0
            JOBE  = 'G'
         END IF
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'EQUIL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'EQUIL must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ITMP = TEMP
         IF ( ITMP.NE.1 .AND. ITMP.NE.2 ) THEN
              CALL mexErrMsgTxt( 'EQUIL must be 1 or 2' )
         END IF
         IF ( ITMP.EQ.2 ) THEN
            EQUIL = 'S'
         ELSE
            EQUIL = 'N'
         END IF
         IP = IP + 1
      ELSE
         ITMP  = 1
         EQUIL = 'N'
      END IF
C
      FPEAK(1) = ZERO
      FPEAK(2) = ONE
      IF ( NRHS.GE.IP ) THEN
         ISIZE = mxGetM( PRHS(IP) ) * mxGetN( PRHS(IP) )
         IF ( ISIZE.GT.2 ) THEN
            CALL mexErrMsgTxt
     $         ( 'FPEAK0 must be a vector with at most 2 elements' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'FPEAK0 must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), FPK, ISIZE )
         IF ( ISIZE.GT.0 ) FPEAK(1) = FPK(1)
         IF ( ISIZE.GT.1 ) FPEAK(2) = FPK(2)
         IP = IP + 1
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'TOL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'TOL must be a numeric scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TOL, 1 )
         IF ( TOL.LE.ZERO ) THEN
              CALL mexErrMsgTxt( 'TOL must be a positive scalar' )
         END IF
         TOL = MAX( HUNDRD*DLAMCH( 'Eps' ), TOL )
      ELSE
         TOL = SQRT( DLAMCH( 'Eps' ) )
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms
C in DGGES, etc.
C
      NN  = N*N
      LDA = MAX( 1, N )
      LDB = LDA
      LDC = MAX( 1, P )
      LDD = LDC
      IF ( ISYSD.EQ.0 ) THEN
         LDE  = LDA
         ICOL = N
      ELSE
         LDE  = 1
         ICOL = 1
      END IF
      LIWORK = N
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( A( LDA, N ), E( LDE, ICOL ), B( LDB, M ), C( LDC, N ),
     $           D( LDD, M ) )
      ALLOCATE ( IWORK( LIWORK ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), A, NN )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), C, P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), D, P*M )
      IF ( .NOT.UNITE .AND. LSAME( JOBE, 'G' ) ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), E, NN )
C
C Test whether E is the unit matrix.
C
         IF ( N.EQ.0 .OR. MA02HD( 'Full', N, N, ONE, E, LDE ) ) THEN
            JOBE  = 'I'
            ISYSD = 1
            DEALLOCATE ( E )
            LDE  = 1
            ICOL = 1
            ALLOCATE ( E( LDE, ICOL ) )
         END IF
      END IF
C
C Test whether D is the zero matrix.
C
      MINPM  = MIN( P, M )
      IF ( MINPM.EQ.0 .OR. MA02HD( 'Full', P, M, ZERO, D, LDD ) ) THEN
         JOBD = 'Z'
      ELSE
         JOBD = 'D'
      END IF
      WITHD = LSAME( JOBD, 'D' )
C
C Test whether B is the zero matrix.
C
      BNORM = DLANGE( '1-norm', N, M, B, LDB, A )
C
C Test whether C is the zero matrix.
C
      CNORM = DLANGE( '1-norm', P, N, C, LDC, A )
C
C Compute LDWORK and allocate DWORK .
C
      N2     = 2*N
      PM     = P  + M
      N2PM   = N2 + PM
      NODYN  = N.EQ.0 .OR. MIN( BNORM, CNORM ).EQ.ZERO
      USEPEN = ISYSD.EQ.0 .OR. ISYS.EQ.2
C
      ID = 6*MINPM
      IC = MAX( 4*MINPM + MAX( P, M ), ID )
      IF( MINPM.EQ.0 ) THEN
         MINWRK = 1
      ELSE IF( NODYN ) THEN
         IF( WITHD ) THEN
            MINWRK = P*M + IC
         ELSE
            MINWRK = 1
         END IF
      ELSE
         IF ( ISYS.EQ.2 ) THEN
            IB = 0
            IE = ID
         ELSE
            IB = N*( N + M )
            IF ( .NOT.WITHD )
     $         IB = IB + P*M
            IE = IC
         END IF
         IF ( WITHD ) THEN
            IR = P*M
            IF ( .NOT.USEPEN ) THEN
               MINWRK = P*P + M*M
               IR = IR + N*PM
            ELSE
               MINWRK = 0
            END IF
            MINWRK = MINWRK + IR + IC
            IR = IR + MINPM
         ELSE
            IR = 0
            MINWRK = 0
         END IF
         IR = IR + N*( N + PM )
         IF ( ISYSD.EQ.0 ) THEN
            IR = IR + NN
            IF ( ITMP.EQ.2 )
     $         MINWRK = MAX( MINWRK, IR + 9*N )
            MINWRK = MAX( MINWRK, IR + 4*N + MAX( M, 2*NN,
     $                                            N + IB + IE ) )
         ELSE
            MINWRK = MAX( MINWRK, IR + N + MAX( M, P, NN + N2,
     $                                          3*N + IB + IE ) )
         END IF
         LW = 0
         IF ( .NOT.USEPEN ) THEN
            LW = IR + 4*NN + 11*N
            IF ( WITHD )
     $         LW = LW + MAX( M, P ) + N*PM
         END IF
         IF ( USEPEN .OR. WITHD )
     $      LW = MAX( LW, IR + 6*N + N2PM*N2PM +
     $                    MAX( N2PM + PM, 8*( NN + N2 ) ) )
         MINWRK = MAX( 1, MINWRK, LW, IR + N2 + IE )
      END IF
C
      LDWORK = MINWRK
      ALLOCATE ( DWORK( LDWORK ) )
C
C Compute LCWORK and allocate CWORK .
C
      IF ( NODYN ) THEN
         LCWORK = 1
      ELSE
         LCWORK = MAX( 1, ( N + M )*( N + P ) + 2*MINPM + MAX( P, M ) )
      END IF
      ALLOCATE ( CWORK( LCWORK ) )
C
C Do the actual computations.
C
      CALL AB13DD( DICO, JOBE, EQUIL, JOBD, N, M, P, FPEAK, A, LDA,
     $             E, LDE, B, LDB, C, LDC, D, LDD, GPEAK, TOL, IWORK,
     $             DWORK, LDWORK, CWORK, LCWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( 2, 1, 0 )
      CALL mxCopyReal8ToPtr( GPEAK, mxGetPr( PLHS(1) ), 2 )
      IF( NLHS.GT.1 ) THEN
         PLHS(2) = mxCreateDoubleMatrix( 2, 1, 0 )
         CALL mxCopyReal8ToPtr( FPEAK, mxGetPr( PLHS(2) ), 2 )
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( A, E, B, C, D, IWORK, DWORK, CWORK )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB13DD'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of LINORM ***
      END
