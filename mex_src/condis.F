#include "fintrf.h"
C CONDIS.F   - Gateway function to perform a transformation on the
C              parameters (A,B,C,D) of a system (equivalent to a
C              bilinear transformation of the corresponding transfer
C              function matrix), using SLICOT routine AB04MD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [Ao,Bo,Co,Do] = CONDIS(task,A,B,C,D(,alpha,beta))
C
C Purpose:
C   To perform a transformation on the parameters (A,B,C,D) of a
C   system, which is equivalent to a bilinear transformation of the
C   corresponding transfer function matrix.
C
C Input parameters:
C   task   - integer specifying the type of the transformation to be
C            performed.
C            task = 1 :  discrete-time   -> continuous-time;
C            task = 2 :  continuous-time -> discrete-time.
C   A      - the n-by-n state dynamics matrix A.
C   B      - the n-by-m input/state matrix B.
C   C      - the p-by-n state/output matrix C.
C   D      - the p-by-m input/output matrix D.
C   alpha, - nonzero parameters specifying the bilinear transformation.
C   beta     Recommended values for stable systems: alpha = 1, beta = 1.
C            Default:  alpha = 1, beta = 1.
C
C Output parameters:
C                                                         _
C   Ao     - the n-by-n transformed state dynamics matrix A.
C                                                       _
C   Bo     - the n-by-m transformed input/state  matrix B.
C                                                       _
C   Co     - the p-by-n transformed state/output matrix C.
C                                                       _
C   Do     - the p-by-m transformed input/output matrix D.
C
C Method:
C   The parameters of the discrete-time system are transformed into
C   the parameters of the continuous-time system (task = 1), or
C   vice-versa (task = 2) by the transformation:
C
C   1.  Discrete -> continuous
C       _                     -1
C       A = beta*(alpha*I + A)  * (A - alpha*I)
C       _                                     -1
C       B = sqrt(2*alpha*beta) * (alpha*I + A)  * B
C       _                                         -1
C       C = sqrt(2*alpha*beta) * C * (alpha*I + A)
C       _                        -1
C       D = D - C * (alpha*I + A)  * B
C
C   which is equivalent to the bilinear transformation
C
C                     z - alpha
C       z -> s = beta ---------
C                     z + alpha
C
C   of one transfer matrix onto the other.
C
C   2.  Continuous -> discrete
C       _                     -1
C       A = alpha*(beta*I - A)  * (beta*I + A)
C       _                                    -1
C       B = sqrt(2*alpha*beta) * (beta*I - A)  * B
C       _                                        -1
C       C = sqrt(2*alpha*beta) * C * (beta*I - A)
C       _                       -1
C       D = D + C * (beta*I - A)  * B
C
C   which is equivalent to the bilinear transformation
C
C                      beta + s
C       s -> z = alpha --------
C                      beta - s
C
C   of one transfer matrix onto the other.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, July 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009,
C   Dec. 2010.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C     .. Parameters ..
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
      CHARACTER         SYST
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDWORK, M, N, P
      DOUBLE PRECISION  ALPHA, BETA
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      mwPointer         A, B, C, D, DWORK
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           IP, ITASK
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          AB04MD
C .. Intrinsic functions ..
      INTRINSIC         MAX
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CONDIS requires at least 5 input arguments' )
      ELSE IF ( NLHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CONDIS requires at least 4 output arguments' )
      END IF
C
C   task, A(nxn), B(n,m), C(p,n), D(p,m)(, alpha, beta))
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
      ITASK = TEMP
      IF ( ITASK.LT.1 .OR. ITASK.GT.2 ) THEN
         CALL mexErrMsgTxt
     $     ( 'TASK has 1 or 2 the only admissible values' )
      END IF
C
      N = mxGetM( PRHS(2) )
      M = mxGetN( PRHS(3) )
      P = mxGetM( PRHS(4) )
      IF ( mxGetN( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same rowsize as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
C
      IF ( mxGetN( PRHS(4) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'C must have the same columnsize as A' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(5) ).NE.P ) THEN
         CALL mexErrMsgTxt( 'D must have the same rowsize as C' )
      END IF
      IF ( mxGetN( PRHS(5) ).NE.M ) THEN
         CALL mexErrMsgTxt( 'D must have the same columnsize as B' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a real matrix' )
      END IF
C
      IF ( ITASK.EQ.1 ) THEN
         SYST = 'D'
      ELSE
         SYST = 'C'
      END IF
C
C   alpha, beta
C
      IP = 6
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), ALPHA, 1 )
         IF ( ALPHA.EQ.ZERO ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be nonzero' )
         END IF
         IP = IP + 1
C
         IF ( NRHS.GE.IP ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'BETA must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'BETA must be a real scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), BETA, 1 )
            IF ( BETA.EQ.ZERO ) THEN
               CALL mexErrMsgTxt( 'BETA must be nonzero' )
            END IF
         ELSE
            BETA = ONE
         END IF
      ELSE
         ALPHA = ONE
         BETA  = ONE
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA = MAX( 1, N )
      LDB = LDA
      LDC = MAX( 1, P )
      LDD = LDC
C
C   ldwork
C   Use larger values, for optimum performance in AB04MD.
C
      LDWORK = LDA
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( IWORK( N ) )
      A     = mxCalloc(  LDA*N, 8 )
      B     = mxCalloc(  LDB*M, 8 )
      C     = mxCalloc(  LDC*N, 8 )
      D     = mxCalloc(  LDD*M, 8 )
      DWORK = mxCalloc( LDWORK, 8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( A ), N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL( B ), N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL( C ), P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), %VAL( D ), P*M )
C
C Do the actual computations.
C
      CALL AB04MD( SYST, N, M, P, ALPHA, BETA, %VAL( A ), LDA,
     $             %VAL( B ), LDB, %VAL( C ), LDC, %VAL( D ), LDD,
     $             IWORK, %VAL( DWORK ), LDWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
         CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS(1) ), N*N )
         PLHS(2) = mxCreateDoubleMatrix( N, M, 0 )
         CALL mxCopyReal8ToPtr( %VAL( B ), mxGetPr( PLHS(2) ), N*M )
         PLHS(3) = mxCreateDoubleMatrix( P, N, 0 )
         CALL mxCopyReal8ToPtr( %VAL( C ), mxGetPr( PLHS(3) ), P*N )
         PLHS(4) = mxCreateDoubleMatrix( P, M, 0 )
         CALL mxCopyReal8ToPtr( %VAL( D ), mxGetPr( PLHS(4) ), P*M )
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
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB04MD'')' )
     $         INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of CONDIS ***
      END
