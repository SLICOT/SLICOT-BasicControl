C HESSOL.F - Gateway function for analysing and solving a system of
C            linear equations with an upper Hessenberg coefficient
C            matrix using SLICOT routines MB02RD, MB02RZ, MB02SD,
C            MB02SZ, MB02TD, and MB02TZ.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [(X)(,rcnd,LU,ipiv,Hn)] = Hessol(job,H(,ipiv,Hn,mtype,normh,B,tran))
C
C   [rcnd]                = Hessol(-2,LU,ipiv,Hn(,mtype,normh))
C   [X]                   = Hessol(-1,LU,ipiv,mtype,B(,tran))
C   [LU,ipiv(,Hn)]        = Hessol( 0,H(,mtype,normh))
C   [X(,LU,ipiv,Hn)]      = Hessol( 1,H,mtype,normh,B(,tran))
C   [rcnd(,LU,ipiv,Hn)]   = Hessol( 2,H(,mtype,normh))
C   [X,rcnd(,LU,ipiv,Hn)] = Hessol( 3,H,mtype,normh,B(,tran))
C
C Purpose:
C  To perform analysis and solution of a system of linear equations
C     H * X = B,  H' * X = B  or  H**H * X = B,
C  with a real or complex upper Hessenberg coefficient matrix H.
C  An LU factorization with row pivoting, H = P*L*U, is used, where
C  P is a permutation matrix, L is lower triangular with unit diagonal
C  elements (and one nonzero subdiagonal), and U is upper triangular.
C
C Input parameters:
C   job    - option parameter indicating the task to be performed.
C            =-2 :  condition estimation using LU factorization of H;
C            =-1 :  solution of the system using LU factorization of H;
C            = 0 :  LU factorization of H;
C            = 1 :  LU factorization of H and solution of the system;
C            = 2 :  LU factorization of H and condition estimation;
C            = 3 :  LU factorization of H, condition estimation, and
C                   solution of the system.
C   H      - the n-by-n Hessenberg matrix H or its LU factorization
C            (if job < 0).
C   ipiv   - (optional) if job < 0, the pivot indices defining P in the
C            LU factorization; for 1 <= i <= n, row i of the matrix H
C            was interchanged with row ipiv(i) of H.
C   Hn     - if job = -2, the 1-norm or the infinity-norm of the initial
C            matrix H, as specified by normh.
C   mtype  - (optional) option parameter indicating the type of the
C            matrix H.
C            = 0 :  real matrix;
C            = 1 :  complex matrix.
C            Default:  mtype = 0.
C   normh  - (optional) if job = -2, or job >= 0, specifies whether the
C            1-norm or the infinity-norm reciprocal condition number is
C            required.
C            = 1 :  1-norm;
C            = 2 :  infinity-norm.
C            Default:  normh = 1.
C   B      - (optional) if |job| = 1, or job = 3, the n-by-m matrix B.
C   tran   - (optional) if |job| = 1, or job = 3,  option parameter
C            indicating whether the matrix H or its (conjugate) transpose
C            should be used.
C            = 0 :  use the matrix H;
C            = 1 :  use the matrix H';
C            = 2 :  use the matrix H**H (if mtype = 1).
C            Default:  tran = 0.
C
C Output parameters:
C   X      - if |job| = 1, or job = 3, the n-by-m solution matrix.
C   rcnd   - if |job| = 2, or job = 3, the reciprocal of the condition
C            number of the matrix H, approximating (in the norm defined
C            by normh) 1/(norm(H) * norm(inv(H))).
C   LU     - if job >= 0, the n-by-n matrix containing the factors
C            L and U from the factorization H = P*L*U.
C   ipiv   - (optional) if job >= 0, the pivot indices defining P in the
C            LU factorization; for 1 <= i <= n, row i of the matrix H
C            was interchanged with row ipiv(i) of H.
C   Hn     - if job >= 0, the 1-norm or the infinity-norm of the initial
C            matrix H, as specified by normh.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2002.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Jan. 2005,
C   Apr. 2009, Dec. 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C
C .. Mex-file interface parameters ..
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mxGetPi, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         NORM, TRANS
      INTEGER           INFO, LDB, LDH, N, NRH
      DOUBLE PRECISION  HNORM, RCND
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IPIV(:), IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: B(:,:),  H(:,:),  DWORK(:)
      COMPLEX*16,       ALLOCATABLE :: BZ(:,:), HZ(:,:), ZWORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           ICON, ISOL, ISREAL, SING
      INTEGER           I, IB, IJOB, INRM, IP, ITRAN, JOB, M, MB, MTYPE
      DOUBLE PRECISION  TEMP
C
C .. External functions ..
      INTEGER           IDAMAX, IZAMAX
      DOUBLE PRECISION  DLANHS, ZLANHS
      EXTERNAL          DLANHS, IDAMAX, IZAMAX, ZLANHS
C
C .. External subroutines ..
      EXTERNAL          MB02RD, MB02RZ, MB02SD, MB02SZ, MB02TD, MB02TZ
C
C .. Intrinsic functions ..
      INTRINSIC         ABS, MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HESSOL requires at least 2 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HESSOL requires at least 1 output argument' )
      END IF
C
C   job, H(nxn), (ipiv(n), Hn, mtype, (normh,) B(nxm), tran).
C
      M = mxGetM( PRHS(2) )
      N = mxGetN( PRHS(2) )
C
      IF ( M.NE.N ) THEN
         CALL mexErrMsgTxt( 'H must be a square matrix' )
      END IF
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR.
     $     mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      JOB  = TEMP
      IJOB = ABS( JOB )
      IF ( JOB.LT.-2 ) THEN
         CALL mexErrMsgTxt
     $           ( 'JOB must be larger than or equal to -2' )
      ELSE IF ( JOB.GT.3 ) THEN
         CALL mexErrMsgTxt
     $           ( 'JOB must be less than or equal to 3' )
      END IF
C
      IF ( ( JOB.EQ.3 .OR. JOB.EQ.0 ) .AND. NLHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HESSOL requires at least 2 output arguments' )
      END IF
C
      ISOL = IJOB.EQ.1 .OR. JOB.EQ.3
      ICON = IJOB.EQ.2 .OR. JOB.EQ.3
C
      IF ( JOB.LT.0 ) THEN
         IF ( JOB.EQ.-2 .AND. NRHS.LT.4 ) THEN
            CALL mexErrMsgTxt
     $           ( 'HESSOL requires at least 4 input arguments' )
         ELSE IF ( NRHS.LT.5 ) THEN
            CALL mexErrMsgTxt
     $           ( 'HESSOL requires at least 5 input arguments' )
         END IF
C
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(3) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'IPIV must be an integer vector' )
         END IF
         MB = mxGetM( PRHS(3) )
         IB = mxGetN( PRHS(3) )
         IF ( MAX( MB, IB ).NE.N .AND.
     $        MIN( MB, IB ).NE.MIN( N, 1 ) ) THEN
            CALL mexErrMsgTxt
     $         ( 'IPIV must have the length equal to the order of H' )
         END IF
         IP = 4
C
         IF ( JOB.EQ.-2 ) THEN
C
C   Hn
C
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'HN must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'HN must be a real scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), HNORM, 1 )
            IF ( HNORM.LT.ZERO ) THEN
               CALL mexErrMsgTxt
     $            ( 'HN must be greater than or equal to zero' )
            END IF
            IP = IP + 1
         END IF
C
      ELSE
         IP = 3
      END IF
C
      IF ( NRHS.GE.IP ) THEN
C
C   mtype
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'MTYPE must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'MTYPE must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         MTYPE = TEMP
         IF ( MTYPE.LT.0 .OR. MTYPE.GT.1 ) THEN
            CALL mexErrMsgTxt
     $         ( 'MTYPE has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
      ELSE
         MTYPE = 0
      END IF
      ISREAL = MTYPE.EQ.0
C
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 )
     $   CALL mexErrMsgTxt( 'H must be a numeric matrix' )
      IF ( mxIsComplex( PRHS(2) ).EQ.1 .AND. ISREAL )
     $   CALL mexErrMsgTxt( 'H must be a real matrix' )
C
      IF ( ( JOB.GE.0 .OR. ICON ) .AND. NRHS.GE.IP ) THEN
C
C   normh
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'NORMH must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'NORMH must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         INRM = TEMP
         IF ( INRM.LT.1 .OR. INRM.GT.2 ) THEN
            CALL mexErrMsgTxt
     $         ( 'NORMH has 1 or 2 the only admissible values' )
         END IF
         IP = IP + 1
      ELSE
         INRM = 1
      END IF
C
      IF ( ISOL ) THEN
         IF ( NRHS.GE.IP ) THEN
C
C   B
C
            MB  = mxGetM( PRHS(IP) )
            NRH = mxGetN( PRHS(IP) )
C
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 )
     $         CALL mexErrMsgTxt( 'B must be a numeric matrix' )
            IF ( mxIsComplex( PRHS(IP) ).EQ.1 .AND. ISREAL )
     $         CALL mexErrMsgTxt( 'B must be a real matrix' )
            IF ( MB.NE.N ) THEN
               CALL mexErrMsgTxt( 'B must have the same row size as H' )
            END IF
            IB = IP
            IP = IP + 1
         ELSE
            CALL mexErrMsgTxt
     $           ( 'HESSOL requires at least 5 input arguments' )
         END IF
C
         IF ( NRHS.GE.IP ) THEN
C
C   tran
C
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'TRAN must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'TRAN must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            ITRAN = TEMP
            IF ( ITRAN.LT.0 .OR. ITRAN.GT.2 ) THEN
               CALL mexErrMsgTxt
     $            ( 'TRAN has 0, 1, or 2 the only admissible values' )
            END IF
         ELSE
            ITRAN = 0
         END IF
      END IF
C
      IF ( INRM.EQ.1 ) THEN
         NORM = '1'
      ELSE
         NORM = 'I'
      END IF
C
      IF ( ISOL ) THEN
         IF ( ITRAN.EQ.0 ) THEN
            TRANS = 'N'
         ELSE IF ( ITRAN.EQ.1 ) THEN
            TRANS = 'T'
         ELSE
            TRANS = 'C'
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C
      LDH = MAX( 1, N )
      IF ( ISOL )
     $   LDB = LDH
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( IPIV( N ) )
      IF ( ISREAL ) THEN
         ALLOCATE ( H( LDH, N ) )
         IF ( ISOL )
     $      ALLOCATE ( B( LDB, NRH ) )
         IF ( ICON )
     $      ALLOCATE ( DWORK( 3*N ), IWORK( N ) )
      ELSE
         ALLOCATE ( HZ( LDH, N ) )
         IF ( ISOL )
     $      ALLOCATE ( BZ( LDB, NRH ) )
         IF ( ICON )
     $      ALLOCATE ( DWORK( N ), ZWORK( 2*N ) )
      END IF
      IF ( .NOT.ICON )
     $   ALLOCATE ( DWORK( N ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      IF ( ISREAL ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), H, N*N )
         IF ( ISOL )
     $      CALL mxCopyPtrToReal8( mxGetPr( PRHS(IB) ), B, N*NRH )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS(2) ),
     $                              mxGetPi( PRHS(2) ), HZ, N*N )
         IF ( ISOL )
     $      CALL mxCopyPtrToComplex16( mxGetPr( PRHS(IB) ),
     $                                 mxGetPi( PRHS(IB) ), BZ, N*NRH )
      END IF
      IF ( JOB.LT.0 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), DWORK, N*MIN( N, 1 )
     $                        )
         DO 10 I = 1, N
            IPIV(I) = DWORK(I)
   10    CONTINUE
      END IF
C
C Do the actual computations.
C
      SING = .FALSE.
      IF ( JOB.GE.0 ) THEN
         IF ( ISREAL ) THEN
            HNORM = DLANHS( NORM, N, H, LDH, DWORK )
         ELSE
            HNORM = ZLANHS( NORM, N, HZ, LDH, DWORK )
         END IF
C
         IF ( ISREAL ) THEN
            CALL MB02SD( N, H, LDH, IPIV, INFO )
         ELSE
            CALL MB02SZ( N, HZ, LDH, IPIV, INFO )
         END IF
         IF ( INFO.NE.0 ) THEN
            RCND = ZERO
            SING = .TRUE.
         END IF
         IF ( .NOT.SING .AND. JOB.GT.0 ) THEN
            IF ( ISOL ) THEN
               IF ( ISREAL ) THEN
                  CALL MB02RD( TRANS, N, NRH, H, LDH, IPIV, B, LDB,
     $                         INFO )
               ELSE
                  CALL MB02RZ( TRANS, N, NRH, HZ, LDH, IPIV, BZ, LDB,
     $                         INFO )
               END IF
            END IF
            IF ( ICON ) THEN
               IF ( ISREAL ) THEN
                  CALL MB02TD( NORM, N, HNORM, H, LDH, IPIV, RCND,
     $                         IWORK, DWORK, INFO )
               ELSE
                  CALL MB02TZ( NORM, N, HNORM, HZ, LDH, IPIV, RCND,
     $                         DWORK, ZWORK, INFO )
               END IF
            END IF
         END IF
      ELSE IF ( ISOL ) THEN
         IF ( ISREAL ) THEN
            INFO = IDAMAX( N, H, LDH+1 )
            IF ( INFO.GT.0 ) THEN
               IF ( H(INFO,INFO).EQ.ZERO ) THEN
                  SING = .TRUE.
                  GO TO 30
               END IF
            END IF
            CALL MB02RD( TRANS, N, NRH, H, LDH, IPIV, B, LDB, INFO )
         ELSE
            INFO = IZAMAX( N, HZ, LDH+1 )
            IF ( INFO.GT.0 ) THEN
               IF ( ABS( HZ(INFO,INFO) ).EQ.ZERO ) THEN
                  SING = .TRUE.
                  GO TO 30
               END IF
            END IF
            CALL MB02RZ( TRANS, N, NRH, HZ, LDH, IPIV, BZ, LDB, INFO )
         END IF
      ELSE
         IF ( ISREAL ) THEN
            CALL MB02TD( NORM, N, HNORM, H, LDH, IPIV, RCND, IWORK,
     $                   DWORK, INFO )
         ELSE
            CALL MB02TZ( NORM, N, HNORM, HZ, LDH, IPIV, RCND, DWORK,
     $                   ZWORK, INFO )
         END IF
      END IF
C
C Copy output to MATLAB workspace.
C
      IP = 1
      IF ( ISOL .AND. .NOT.SING ) THEN
         IF ( ISREAL ) THEN
            PLHS(1) = mxCreateDoubleMatrix( N, NRH, 0 )
            CALL mxCopyReal8ToPtr( B, mxGetPr( PLHS(1) ), N*NRH )
         ELSE
            PLHS(1) = mxCreateDoubleMatrix( N, NRH, 1 )
            CALL mxCopyComplex16ToPtr( BZ, mxGetPr( PLHS(1) ),
     $                                 mxGetPi( PLHS(1) ), N*NRH )
         END IF
         IP = IP + 1
      END IF
      IF ( ( JOB.EQ.3 .AND. NLHS.GE.IP ) .OR. IJOB.EQ.2 ) THEN
         PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
         CALL mxCopyReal8ToPtr( RCND, mxGetPr( PLHS(IP) ), 1 )
         IP = IP + 1
      END IF
      IF ( JOB.GE.0 .AND. NLHS.GE.IP ) THEN
         IF ( ISREAL ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( H, mxGetPr( PLHS(IP) ), N*N )
         ELSE
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 1 )
            CALL mxCopyComplex16ToPtr( HZ, mxGetPr( PLHS(IP) ),
     $                                 mxGetPi( PLHS(IP) ), N*N )
         END IF
         IP = IP + 1
         IF ( NLHS.GE.IP ) THEN
            DO 20 I = 1, N
               DWORK(I) = IPIV(I)
   20       CONTINUE
            PLHS(IP) = mxCreateDoubleMatrix( N, MIN( N, 1 ), 0 )
            CALL mxCopyReal8ToPtr( DWORK, mxGetPr( PLHS(IP) ),
     $                             N*MIN( N, 1 ) )
            IP = IP + 1
         END IF
         IF ( NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( HNORM, mxGetPr( PLHS(IP) ), 1 )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
   30 CONTINUE
C
      DEALLOCATE ( IPIV )
      IF ( ISREAL ) THEN
         DEALLOCATE ( H )
         IF ( ISOL )
     $      DEALLOCATE ( B )
         IF ( ICON )
     $      DEALLOCATE ( DWORK, IWORK )
      ELSE
         DEALLOCATE ( HZ )
         IF ( ISOL )
     $      DEALLOCATE ( BZ )
         IF ( ICON )
     $      DEALLOCATE ( DWORK, ZWORK )
      END IF
      IF ( .NOT.ICON )
     $   DEALLOCATE ( DWORK )
C
C Error and warning handling.
C
      IF ( SING ) THEN
         WRITE( TEXT, '('' Warning: Matrix H is singular.'')' )
         CALL mexPrintf( TEXT )
         IF ( ISOL ) THEN
            WRITE( TEXT, '('' The solution is infinite.'')' )
            CALL mexErrMsgTxt( TEXT )
         END IF
      END IF
C
      RETURN
C *** Last line of HESSOL ***
      END
