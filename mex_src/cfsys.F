#include "fintrf.h"
C CFSYS.F    - Gateway function for constructing the state-space
C              representation for the system G = (A,B,C,D) from the
C              factors Q and R of its left or right coprime
C              factorization, where G, Q and R are the corresponding
C              transfer-function matrices, using SLICOT routines SB08GD
C              and SB08HD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [Ao,Bo,Co,Do(,rcnd)] = cfsys(task,A,B,C,D(,BR)(,CR),DR)
C   [Ao,Bo,Co,Do(,rcnd)] = cfsys(1,A,B,C,D,BR,DR)
C   [Ao,Bo,Co,Do(,rcnd)] = cfsys(2,A,B,C,D,CR,DR)
C
C Purpose:
C   To construct the state-space representation for the system
C   G = (A,B,C,D) from the factors Q = (AQR,BQ,CQR,DQ) and
C   R = (AQR,BR,CQR,DR) of its left coprime factorization,
C                   -1
C              G = R  * Q,
C
C   or from the factors Q = (AQR,BQR,CQ,DQ) and R = (AQR,BQR,CR,DR)
C   of its right coprime factorization,
C                       -1
C              G = Q * R  ,
C
C   where G, Q and R are the corresponding transfer-function matrices.
C
C Input parameters:
C   task   - integer specifying the computations to be performed.
C            = 1 :  compute the state-space representation from the
C                   factors Q and R of its left coprime factorization;
C            = 2 :  compute the state-space representation from the
C                   factors Q and R of its right coprime factorization.
C   A      - the n-by-n state dynamics matrix AQR of the systems
C            Q and R.
C   B      - if task = 1, the n-by-m input/state matrix BQ of the
C            system Q;
C          - if task = 2, the n-by-m input/state matrix BQR of the
C            systems Q and R.
C   C      - if task = 1, the p-by-n state/output matrix CQR of the
C            systems Q and R;
C          - if task = 2, the p-by-n state/output matrix CQ of the
C            system Q.
C   D      - the p-by-m input/output matrix DQ of the system Q.
C   BR     - if task = 1, the n-by-p input/state matrix BR of the
C            system R.
C   CR     - if task = 2, the m-by-n state/output matrix CR of the
C            system R.
C   DR     - if task = 1, the p-by-p input/output matrix DR of the
C            system R.
C          - if task = 2, the m-by-m input/output matrix DR of the
C            system R.
C
C Output parameters:
C   Ao     - the n-by-n state dynamics matrix of the system G.
C   Bo     - the n-by-m input/state matrix of the system G.
C   Co     - the p-by-n state/output matrix of the system G.
C   Do     - the p-by-m input/output matrix of the system G.
C   rcnd   - (optional) an estimate of the reciprocal condition number
C            of the matrix DR.
C
C Method:
C   If task = 1, the matrices of the state-space representation
C   G = (A,B,C,D) are computed by using the formulas:
C
C                      -1              -1
C     A = AQR - BR * DR  * CQR,  C = DR  * CQR,
C                      -1              -1
C     B = BQ  - BR * DR  * DQ,   D = DR  * DQ.
C
C   If task = 2, the matrices of the state-space representation
C   G = (A,B,C,D) are computed by using the formulas:
C
C                       -1                   -1
C     A = AQR - BQR * DR  * CR,  B = BQR * DR  ,
C                      -1                   -1
C     C = CQ  - DQ * DR  * CR,   D = DQ * DR  .
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Jan. 2005,
C   Apr. 2009, Dec. 2010.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
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
      INTEGER           INFO, LDA, LDB, LDBR, LDC, LDCR, LDD, LDDR, M,
     $                  N, P
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      mwPointer         A, B, BR, C, CR, D, DR, DWORK
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      CHARACTER*2       NAME
      INTEGER           ITASK, LDWORK, ND
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          SB08GD, SB08HD
C
C .. Intrinsic functions ..
      INTRINSIC         MAX
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.7 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CFSYS requires 7 input arguments' )
      ELSE IF ( NLHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CFSYS requires at least 4 output arguments' )
      END IF
C
C   task, A(nxn), B(nxm), C(pxn), D(pxm), BR(nxp)/CR(mxn),
C   DR(pxp)/DR(mxm).
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
     $    ( 'ITASK has 1 or 2 the only admissible values' )
      END IF
C
      N = mxGetM( PRHS(2) )
      M = mxGetN( PRHS(3) )
      P = mxGetM( PRHS(4) )
C
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF ( mxGetN( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF ( mxGetM( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same number of rows as A' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
      IF ( mxGetN( PRHS(4) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $      ( 'C must have the same number of columns as A' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a real matrix' )
      END IF
      IF ( mxGetM( PRHS(5) ).NE.P ) THEN
         CALL mexErrMsgTxt( 'D must have the same number of rows as C' )
      END IF
      IF ( mxGetN( PRHS(5) ).NE.M ) THEN
         CALL mexErrMsgTxt
     $      ( 'D must have the same number of columns as B' )
      END IF
C
      IF ( ITASK.EQ.1 ) THEN
         NAME = 'BR'
         ND   = P
      ELSE
         NAME = 'CR'
         ND   = M
      END IF
C
      IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(6) ).EQ.1 ) THEN
         WRITE( TEXT, '(A2, '' must be a real matrix'')' ) NAME
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(7) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'DR must be a real matrix' )
      END IF
      IF ( ITASK.EQ.1 ) THEN
         IF ( mxGetM( PRHS(6) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $         ( 'BR must have the same number of rows as A' )
         END IF
         IF ( mxGetN( PRHS(6) ).NE.P ) THEN
            CALL mexErrMsgTxt
     $     ( 'BR must have as many columns as the number of rows of C' )
         END IF
      ELSE
         IF ( mxGetM( PRHS(6) ).NE.M ) THEN
            CALL mexErrMsgTxt
     $     ( 'CR must have as many rows as the number of columns of B' )
         END IF
         IF ( mxGetN( PRHS(6) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $         ( 'CR must have the same number of columns as C' )
         END IF
      END IF
      IF ( mxGetM( PRHS(7) ).NE.ND .OR.
     $     mxGetN( PRHS(7) ).NE.ND ) THEN
         WRITE( TEXT, '('' DR must have the order '',I7)' ) ND
         CALL mexErrMsgTxt( TEXT )
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA = MAX( 1, N )
      LDB = LDA
      LDC = MAX( 1, P )
      LDD = LDC
      IF ( ITASK.EQ.1 ) THEN
         LDBR = LDA
      ELSE
         LDCR = MAX( 1, M )
      END IF
      LDDR   = MAX( 1, ND )
      LDWORK = MAX( 1, 4*ND )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( IWORK( ND ) )
      A     = mxCalloc(  LDA*N,  8 )
      B     = mxCalloc(  LDB*M,  8 )
      C     = mxCalloc(  LDC*N,  8 )
      D     = mxCalloc(  LDD*M,  8 )
      DR    = mxCalloc( LDDR*ND, 8 )
      DWORK = mxCalloc( LDWORK,  8 )
      IF ( ITASK.EQ.1 ) THEN
         BR = mxCalloc( LDBR*P, 8 )
      ELSE
         CR = mxCalloc( LDCR*N, 8 )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( A ), N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL( B ), N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL( C ), P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), %VAL( D ), P*M )
      IF ( ITASK.EQ.1 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), %VAL( BR ), N*P )
      ELSE
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), %VAL( CR ), M*N )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), %VAL( DR ), ND*ND )
C
C Do the actual computations.
C
      IF ( ITASK.EQ.1 ) THEN
         CALL SB08GD( N, M, P, %VAL( A ), LDA, %VAL( B ), LDB,
     $                %VAL( C ), LDC, %VAL( D ), LDD, %VAL( BR ), LDBR,
     $                %VAL( DR ), LDDR, IWORK, %VAL( DWORK ), INFO )
      ELSE
         CALL SB08HD( N, M, P, %VAL( A ), LDA, %VAL( B ), LDB,
     $                %VAL( C ), LDC, %VAL( D ), LDD, %VAL( CR ), LDCR,
     $                %VAL( DR ), LDDR, IWORK, %VAL( DWORK ), INFO )
      END IF
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
         IF ( NLHS.GE.5 ) THEN
            PLHS(5) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( %VAL( DWORK ), mxGetPr( PLHS(5) ),
     $                             1 )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE ( IWORK )
      CALL mxFree( A     )
      CALL mxFree( B     )
      CALL mxFree( C     )
      CALL mxFree( D     )
      CALL mxFree( DR    )
      CALL mxFree( DWORK )
      IF ( ITASK.EQ.1 ) THEN
         CALL mxFree( BR )
      ELSE
         CALL mxFree( CR )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( ITASK.EQ.1 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB08GD'')'
     $           ) INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB08HD'')'
     $           ) INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of CFSYS ***
      END
