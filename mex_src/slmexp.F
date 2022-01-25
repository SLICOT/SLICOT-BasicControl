#include "fintrf.h"
C SLMEXP.F   - Gateway function to compute the exponential of a matrix
C              and, optionally, its integral, using SLICOT routines
C              MB05MD, MB05ND, and MB05OD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [F(,V,Y,Wr,Wi)(,mdig,idig)(,H)] = SLMEXP(meth,A(,delta)(,scale)
C                                                  (,ordPad)(,tol))
C
C   [F(,V,Y,Wr,Wi)] = SLMEXP(0,A(,delta,scale))
C   [F(,V,Y,sw,sw)] = SLMEXP(0,A(,delta,scale))
C   [F(,mdig,idig)] = SLMEXP(1,A(,delta,scale,ordPad))
C   [F(,H)]         = SLMEXP(2,A(,delta)(,tol))
C
C Purpose:
C   To compute F(delta) = exp(A*delta), where A is a real n-by-n matrix
C   and delta is a scalar value. Optionally (for meth = 2), the matrix
C   integral H(delta), defined by
C
C      H(delta) = Int[F(s) ds] from s = 0 to s = delta,
C
C   is returned. Other useful results can also be obtained if required.
C   When meth = 0, matrix A is assumed to be non-defective.
C
C Input parameters:
C   meth   - scalar indicating the method to be used for computing the
C            matrix exponential, as follows:
C            = 0 :  use an eigenvalue/eigenvector decomposition
C                   technique;
C            = 1 :  use a diagonal Pade approximant with scaling and
C                   squaring;
C            = 2 :  use a Pade approximation.
C   A      - the n-by-n matrix A. If delta = 0, the given A is not used.
C   delta  - (optional) the scalar value delta of the problem.
C            Default:  delta = 1.
C   scale  - (optional) if meth = 0 or 1, scalar indicating whether or
C            not the matrix should be diagonally scaled, as follows:
C            = 0 :  do not scale the matrix;
C            = 1 :  diagonally scale the matrix, i.e., replace A by
C                   D*A*D**(-1), where D is a diagonal matrix chosen to
C                   make the rows and columns of A more equal in norm.
C            Default:  scale = 1.
C   ordPad - (optional) if meth = 1, scalar specifying the order of the
C            diagonal Pade approximant. 1 <= ordPad <= 15. In the
C            absence of further information, ordPad should be set to 9.
C            Default:  ordPad = 9.
C   tol    - (optional) if meth = 2, real scalar indicating the
C            tolerance to be used in determining the order of the Pade
C            approximation to H(t), where t is a scale factor determined
C            internally. A reasonable value for tol may be sqrt(EPS),
C            where EPS is the machine precision (see LAPACK Library
C            routine DLAMCH).
C            Default:  tol = sqrt(EPS).
C
C Output parameters:
C   F      - the n-by-n solution matrix exp(A*delta).
C   V      - (optional) if meth = 0, the n-by-n eigenvector matrix
C            for A. If the k-th eigenvalue is real, the k-th column of
C            the eigenvector matrix holds the eigenvector corresponding
C            to the k-th eigenvalue.
C            Otherwise, the k-th and (k+1)-th eigenvalues form a
C            complex conjugate pair and the k-th and (k+1)-th columns
C            of the eigenvector matrix hold the real and imaginary
C            parts of the eigenvectors corresponding to these
C            eigenvalues as follows.
C            If p and q denote the k-th and (k+1)-th columns of the
C            eigenvector matrix, respectively, then the eigenvector
C            corresponding to the complex eigenvalue with positive
C            (negative) imaginary value is given by
C                                      2
C            p + q*j (p - q*j), where j  = -1.
C   Y      - (optional) if meth = 0, the n-by-n real matrix satisfying
C            the equation V*Y = F. If all eigenvalues are real, then Y
C            is the matrix product exp(Lambda*delta) times the inverse
C            of the (right) eigenvector matrix V of A, where Lambda is
C            the diagonal matrix of eigenvalues.
C   Wr,Wi  - (optional) if meth = 0, n-vectors of real and imaginary
C            parts, respectively, of the eigenvalues of the matrix A.
C            The eigenvalues are unordered except that complex conjugate
C            pairs of values appear consecutively, with the eigenvalue
C            having positive imaginary part first.
C   sw     - (optional) if meth = 0, but the matrix A is found to be
C            defective, sw is set to 1, and the calculations are redone
C            automatically for meth = 1 with ordPad = 9. Then, the
C            parameters V and Y (if specified) will contain mdig and
C            idig, respectively. Switching the methods also occurs when
C            sw is not present (the output list has less than four
C            parameters).
C   mdig   - (optional) if meth = 1, the minimal number of accurate
C            digits in the 1-norm of exp(A*delta).
C   idig   - (optional) if meth = 1, the number of accurate digits in
C            the 1-norm of exp(A*delta) at 95% confidence level.
C   H      - (optional) if meth = 2, the n-by-n matrix containing an
C            approximation to H(delta).
C
C Method
C   If meth = 0, an eigenvalue/eigenvector decomposition technique is
C   used, based on a modification of LAPACK Library routine DGEEV for
C   obtaining the right eigenvector matrix. A condition estimate is
C   then employed to determine if the matrix A is near defective and
C   hence the exponential solution is inaccurate. If A is defective, a
C   warning is returned, and meth = 1 is then used. If the output
C   parameter list in that case had 4 items, the parameters are then
C   F, mdig, idig, and sw (set to 1); for 5 items, sw is duplicated.
C
C   If meth = 1, the exponential of the matrix A is evaluated from a
C   diagonal Pade approximant. The algorithm exploits the identity
C
C       (exp[(2**-m)*A]) ** (2**m) = exp(A),
C
C   where m is an integer determined internally, to improve the accuracy
C   for matrices with large norms.
C
C   If meth = 2, a Pade approximation to H(t) for some small value of t
C   (where 0 < t <= delta) is used, and then F(t) is calculated from
C   H(t). Finally, the results are re-scaled to give F(delta) and
C   H(delta).
C
C Comments
C   The general recommended method for computing the matrix exponential
C   is meth = 1. Setting meth = 0 could be more efficient for
C   non-defective matrices.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009,
C   July 2012.
C
C **********************************************************************
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
      CHARACTER         BALANC
      INTEGER           IDIG, INFO, IWARN, LDA, LDEX, LDEXIN, LDV,
     $                  LDWORK, LDY, MDIG, N, NDIAG
      DOUBLE PRECISION  DELTA, TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      mwPointer         A, DWORK, EX, EXINT, V, VALI, VALR, Y
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           PADE
      CHARACTER*120     TEXT
      INTEGER           IMETH, IP, ISCAL
      DOUBLE PRECISION  TEMP
C
C .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C
C .. External subroutines ..
      EXTERNAL          DLASET, MB05MD, MB05ND, MB05OD
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN, SQRT
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SLMEXP requires at least 2 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SLMEXP requires at least 1 output argument' )
      END IF
C
C   meth, A(nxn) (, delta) (, scale) (, ordPad) (,tol)
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR.
     $     mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      IMETH = TEMP
      IF ( IMETH.LT.0 .OR. IMETH.GT.2 ) THEN
         CALL mexErrMsgTxt
     $     ( 'METH has 0, 1 or 2 the only admissible values' )
      END IF
C
      N =  mxGetM( PRHS(2) )
      IF ( mxGetN( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
C   delta
C
      IP = 3
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'DELTA must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'DELTA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), DELTA, 1 )
         IP = IP + 1
      ELSE
         DELTA = ONE
      END IF
C
      IF ( IMETH.LE.1 ) THEN
C
C   scale
C
         IF ( NRHS.GE.IP ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'SCALE must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'SCALE must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            ISCAL = TEMP
            IF ( ISCAL.LT.0 .OR. ISCAL.GT.1 ) THEN
               CALL mexErrMsgTxt
     $              ( 'SCALE has 0 or 1 the only admissible values')
            END IF
            IP = IP + 1
         ELSE
            ISCAL = 1
         END IF
C
         IF ( ISCAL.EQ.0 ) THEN
            BALANC = 'N'
         ELSE
            BALANC = 'S'
         END IF
C
C   ordPad
C
         IF ( NRHS.GE.IP .AND. IMETH.EQ.1 ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'ORDPAD must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'ORDPAD must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            NDIAG = TEMP
            IF ( NDIAG.LT.1 ) THEN
               CALL mexErrMsgTxt
     $              ( 'ORDPAD must be larger than or equal to 1' )
            END IF
            IF ( NDIAG.GT.15 ) THEN
               CALL mexErrMsgTxt
     $              ( 'ORDPAD must be less than or equal to 15' )
            END IF
         ELSE
            NDIAG = 9
         END IF
C
      ELSE
C
C   tol
C
         IF ( NRHS.GE.IP ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'TOL must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'TOL must be a real scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TOL, 1 )
         ELSE
            TOL = SQRT( DLAMCH( 'Epsilon' ) )
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger LDWORK for allowing block algorithms in DGEHRD, DORGHR
C (in MB05MD).
C
      LDA = MAX( 1, N )
      IF ( IMETH.EQ.0 ) THEN
         LDV = LDA
         LDY = LDA
      ELSE IF ( IMETH.EQ.2 ) THEN
         LDEX   = LDA
         LDEXIN = LDA
      END IF
C
C   ldwork
C
      IF ( IMETH.EQ.0 ) THEN
         LDWORK = MAX( 1, 4*N )
      ELSE IF ( IMETH.EQ.1 ) THEN
         IF ( N.GT.1 ) THEN
            LDWORK = N*( 2*N + NDIAG + 1 ) + NDIAG
         ELSE
            LDWORK = 1
         END IF
      ELSE
C        For optimum performance LDWORK should be larger (2*N*N).
         LDWORK = MAX( 1, N*( N + 1 ) )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( IWORK( N ) )
      A     = mxCalloc(  LDA*N, 8 )
      DWORK = mxCalloc( LDWORK, 8 )
      IF ( IMETH.EQ.0 ) THEN
         V    = mxCalloc( LDV*N, 8 )
         VALI = mxCalloc(     N, 8 )
         VALR = mxCalloc(     N, 8 )
         Y    = mxCalloc( LDY*N, 8 )
      ELSE IF ( IMETH.EQ.2 ) THEN
         EX    = mxCalloc(   LDEX*N, 8 )
         EXINT = mxCalloc( LDEXIN*N, 8 )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL(A), N*N )
C
C Do the actual computations.
C
      IWARN = 0
      PADE  = .FALSE.
      IF ( IMETH.EQ.0 ) THEN
         IF ( DELTA.EQ.ZERO )
     $      CALL DLASET( 'Full', N, N, ZERO, ONE, %VAL(A), LDA )
C
         CALL MB05MD( BALANC, N, DELTA, %VAL(A), LDA, %VAL(V), LDV,
     $                %VAL(Y), LDY, %VAL(VALR), %VAL(VALI), IWORK,
     $                %VAL(DWORK), LDWORK, INFO )
         IF ( INFO.GT.0 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB05MD.'',
     $                     '' MB05OD is used instead (ORDPAD = 9).'')'
     $           ) INFO
            CALL mexPrintf( TEXT )
            PADE  = .TRUE.
            NDIAG = 9
            IF ( N.GT.1 ) THEN
               LDWORK = N*( 2*N + NDIAG + 1 ) + NDIAG
            ELSE
               LDWORK = 1
            END IF
            CALL mxFree( DWORK )
            CALL mxFree( V     )
            CALL mxFree( VALI  )
            CALL mxFree( VALR  )
            CALL mxFree( Y     )
            DWORK = mxCalloc( LDWORK, 8 )
            CALL MB05OD( BALANC, N, NDIAG, DELTA, %VAL(A), LDA, MDIG,
     $                   IDIG, IWORK, %VAL(DWORK), LDWORK, IWARN, INFO )
         END IF
      ELSE IF ( IMETH.EQ.2 ) THEN
         CALL MB05ND( N, DELTA, %VAL(A), LDA, %VAL(EX), LDEX,
     $                %VAL(EXINT), LDEXIN, TOL, IWORK, %VAL(DWORK),
     $                LDWORK, INFO )
      ELSE
         CALL MB05OD( BALANC, N, NDIAG, DELTA, %VAL(A), LDA, MDIG, IDIG,
     $                IWORK, %VAL(DWORK), LDWORK, IWARN, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
         IF ( IMETH.LE.1 ) THEN
            CALL mxCopyReal8ToPtr( %VAL(A), mxGetPr( PLHS(1) ),
     $                             N*N )
         ELSE
            CALL mxCopyReal8ToPtr( %VAL(EX), mxGetPr( PLHS(1) ),
     $                             N*N )
         END IF
         IP = 2
C
         IF ( NLHS.GE.IP ) THEN
            IF ( IMETH.EQ.1 .OR. PADE ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
               TEMP = MDIG
               CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS(IP) ), 1 )
               IP = IP + 1
               IF ( NLHS.GE.IP ) THEN
                  PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
                  TEMP = IDIG
                  CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS(IP) ),
     $                                   1 )
                  IP = IP + 1
                  IF ( NLHS.GE.IP .AND. PADE ) THEN
                     PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
                     TEMP = ONE
                     CALL mxCopyReal8ToPtr( TEMP,
     $                                      mxGetPr( PLHS(IP) ),
     $                                      1 )
                     IP = IP + 1
                     IF ( NLHS.GE.IP ) THEN
                        PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
                        CALL mxCopyReal8ToPtr( TEMP,
     $                                         mxGetPr( PLHS(IP) ),
     $                                         1 )
                     END IF
                  END IF
               END IF
            ELSE IF ( IMETH.EQ.0 ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
               CALL mxCopyReal8ToPtr( %VAL(V), mxGetPr( PLHS(IP) ),
     $                                N*N )
               IP = IP + 1
               IF ( NLHS.GE.IP ) THEN
                  PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
                  CALL mxCopyReal8ToPtr( %VAL(Y),
     $                                   mxGetPr( PLHS(IP) ), N*N )
                  IP = IP + 1
                  IF ( NLHS.GE.IP ) THEN
                     PLHS(IP) = mxCreateDoubleMatrix( N,
     $                                                MIN( N, 1 ), 0 )
                     CALL mxCopyReal8ToPtr( %VAL(VALR),
     $                                      mxGetPr( PLHS(IP) ),
     $                                      N*MIN( N, 1 ) )
                     IP = IP + 1
                     IF ( NLHS.GE.IP ) THEN
                        PLHS(IP) = mxCreateDoubleMatrix( N,
     $                                                   MIN( N, 1 ), 0)
                        CALL mxCopyReal8ToPtr( %VAL(VALI),
     $                                         mxGetPr( PLHS(IP) ),
     $                                         N*MIN( N, 1 ) )
                     END IF
                  END IF
               END IF
            ELSE
               PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
               CALL mxCopyReal8ToPtr( %VAL(EXINT),
     $                                mxGetPr( PLHS(IP) ), N*N )
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE(  IWORK )
      CALL mxFree( A     )
      CALL mxFree( DWORK )
      IF ( IMETH.EQ.2 ) THEN
         CALL mxFree( EX    )
         CALL mxFree( EXINT )
      ELSE IF ( IMETH.EQ.0 .AND. .NOT.PADE ) THEN
         CALL mxFree( V    )
         CALL mxFree( VALI )
         CALL mxFree( VALR )
         CALL mxFree( Y    )
      END IF
C
C Error and warning handling.
C
      IF ( IWARN.NE.0 ) THEN
         WRITE( TEXT, '(''  IWARN = '',I4,'' ON EXIT FROM MB05OD'')'
     $        ) IWARN
         CALL mexPrintf( TEXT )
      END IF
C
      IF ( INFO.NE.0 ) THEN
         IF ( IMETH.EQ.1 .OR. PADE ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB05OD'')' )
     $             INFO
         ELSE IF ( IMETH.EQ.0 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB05MD'')' )
     $             INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB05ND'')' )
     $             INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SLMEXP ***
      END
