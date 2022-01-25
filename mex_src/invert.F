#include "fintrf.h"
C INVERT.F   - Gateway function to compute the dual or inverse
C              of a linear (descriptor) system, using SLICOT routines
C              AB07MD, AB07ND, and AG07BD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [Ao(,Eo),Bo,Co(,Do)(,rcondD)] = INVERT(task,A,B,C(,D)(,E))
C
C   [Ao,Bo,Co(,Do)]        = INVERT(1,A,B,C(,D))
C   [Ao,Bo,Co,Do(,rcondD)] = INVERT(2,A,B,C,D)
C   [Ao,Eo,Bo,Co,Do]       = INVERT(3,A,B,C,D(,E))
C
C Purpose:
C   To compute the dual of a standard system or the inverse of a
C   standard system or a descriptor system, according to the value of
C   task:
C
C   task = 1:  To compute the dual of a standard system in state-space
C   representation, i.e., if the m-input/p-output system is (A,B,C,D),
C   its dual is simply the p-input/m-output system (A',C',B',D').
C
C   task = 2:  To compute the inverse of a standard system.
C   The matrices of the inverse system are computed with the formulas
C                    -1              -1         -1           -1
C        Ai = A - B*D  *C,  Bi = -B*D  ,  Ci = D  *C,  Di = D  .     (1)
C
C   task = 3:  To compute the inverse (Ai-lambda*Ei,Bi,Ci,Di) of a given
C   descriptor system (A-lambda*E,B,C,D). The matrices of the inverse
C   system are computed with the formulas
C
C             ( E  0 )        ( A  B )         (  0 )
C        Ei = (      ) , Ai = (      ) ,  Bi = (    ),
C             ( 0  0 )        ( C  D )         ( -I )
C
C        Ci = ( 0  I ),  Di = 0.                                     (2)
C
C Input parameters:
C   task   - integer specifying the computations to be performed.
C            task = 1 :  compute the dual of a standard system;
C            task = 2 :  compute the inverse of a standard system;
C            task = 3 :  compute the inverse of a descriptor system.
C   A      - the n-by-n state dynamics matrix A.
C   B      - the n-by-m input/state matrix B.
C   C      - the p-by-n state/output matrix C. If task >= 2, p = m.
C   D      - the p-by-m input/output matrix D. If task >= 2, p = m.
C            If task = 1 and D is zero, this argument may be omitted.
C   E      - the n-by-n descriptor matrix E. When E is an identity
C            matrix, this argument may be omitted.
C
C Output parameters:
C   Ao     - If task = 1, the n-by-n dual state dynamics matrix A'.
C            If task = 2, the n-by-n state matrix Ai of the inverse
C            system (1).
C            If task = 3, the (n+m)-by-(n+m) state matrix Ai of the
C            inverse system (2).
C   Eo     - If task = 3 the (n+m)-by-(n+m) descriptor matrix Ei of the
C            inverse system (2).
C   Bo     - If task = 1, the n-by-p dual input/state matrix C'.
C            If task = 2, the n-by-m input matrix Bi of the inverse
C            system (1).
C            If task = 3, the (n+m)-by-m input matrix Bi of the
C            inverse system (2).
C   Co     - If task = 1, the m-by-n dual state/output matrix B'.
C            If task = 2, the m-by-n output matrix Ci of the inverse
C            system (1).
C            If task = 3, the m-by-(n+m) output matrix Ci of the
C            inverse system (2).
C   Do     - If task = 1, the m-by-p dual input/output matrix D'.
C            If task = 2, the m-by-m input/output matrix Di of the
C            inverse system (1).
C            If task = 3, the m-by-m input/output (zero) matrix Di of
C            the inverse system (2).
C   rcondD - (optional) the estimated reciprocal condition number of the
C            input/output matrix D of the original system.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, July 2003.
C
C Revisions:
C   V. Sima, Feb. 2004, Apr. 2009, Dec. 2010, July 2012, Sep. 2012.
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
      CHARACTER         JOBD, JOBE
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDE, LDWORK, M, N, P
      DOUBLE PRECISION  RCOND
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      mwPointer         A, B, C, D, DW, DWORK, E
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           WITHD, WITHE
      INTEGER           I, IP, ITASK, J, LIWORK, MP, N1, ND, NM
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          AB07MD, AB07ND, AG07BD, DLACPY
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'INVERT requires at least 4 input arguments' )
      ELSE IF ( NLHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'INVERT requires at least 3 output arguments' )
      END IF
C
C   task, A(nxn), B(n,m), C(p,n)(, D(p,m))(, E(n,n))
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
      IF ( ITASK.LT.1 .OR. ITASK.GT.3 ) THEN
         CALL mexErrMsgTxt
     $     ( 'TASK has 1, 2, or 3 the only admissible values' )
      END IF
C
      NM = 3 + MIN( ITASK, 2 )
      IF ( NRHS.LT.NM ) THEN
         WRITE( TEXT,
     $      '(''INVERT requires at least '',I1,'' input arguments'')' )
     $      NM
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( ITASK.LT.3 )
     $   NM = NM - 1
      IF ( NLHS.LT.NM ) THEN
         WRITE( TEXT,
     $      '(''INVERT requires at least '',I1,'' output arguments'')' )
     $      NM
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      WITHD = NRHS.GE.5
      WITHE = NRHS.GE.6
C
      N = mxGetM( PRHS(2) )
      M = mxGetN( PRHS(3) )
      P = mxGetM( PRHS(4) )
      IF ( ITASK.GE.2 .AND. P.NE.M ) THEN
         CALL mexErrMsgTxt(
     $   'The number of rows of C must equal the number of columns of B'
     $                    )
      END IF
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
      IF ( WITHD ) THEN
         IF ( mxGetM( PRHS(5) ).NE.P ) THEN
            CALL mexErrMsgTxt( 'D must have the same rowsize as C' )
         END IF
         IF ( mxGetN( PRHS(5) ).NE.M ) THEN
            CALL mexErrMsgTxt( 'D must have the same columnsize as B' )
         END IF
         IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(5) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'D must be a real matrix' )
         END IF
      END IF
C
      IF ( WITHE ) THEN
         IF ( mxGetM( PRHS(6) ).NE.N .OR.
     $        mxGetN( PRHS(6) ).NE.N ) THEN
            CALL mexErrMsgTxt( 'E must have the same size as A' )
         END IF
         IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'E must be a real matrix' )
         END IF
      END IF
C
      IF ( ITASK.EQ.1 ) THEN
         IF ( WITHD ) THEN
            JOBD = 'D'
         ELSE
            JOBD = 'Z'
         END IF
      ELSE IF ( ITASK.EQ.3 ) THEN
         IF ( WITHE ) THEN
            JOBE = 'G'
         ELSE
            JOBE = 'I'
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C
      N1 = MAX( 1, N )
      IF ( ITASK.EQ.3 ) THEN
         NM = N + M
      ELSE
         NM = N
      END IF
      LDA = MAX( 1, NM )
      LDB = LDA
      IF ( ITASK.EQ.1 ) THEN
         MP = MAX( M, P )
         IF ( N.GT.0 ) THEN
            LDC = MAX( 1, MP )
         ELSE
            LDC = 1
         END IF
      ELSE
         MP  = M
         LDC = MAX( 1, MP )
         IF ( ITASK.EQ.3 )
     $      LDE = LDA
      END IF
      IF ( WITHD ) THEN
         LDD = MAX( 1, MP )
         ND  = MP
      ELSE
         LDD = 1
         ND  = 1
      END IF
C
C   ldwork
C   Use larger values, for optimum performance in AB07ND.
C
      IF ( ITASK.EQ.2 ) THEN
         LIWORK = 2*M
         LDWORK = MAX( 1, 4*M )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      A = mxCalloc(  LDA*NM, 8 )
      B = mxCalloc(  LDB*MP, 8 )
      C = mxCalloc(  LDC*NM, 8 )
      D = mxCalloc(  LDD*ND, 8 )
      IF ( ITASK.EQ.2 ) THEN
         ALLOCATE ( IWORK( LIWORK ) )
         DWORK = mxCalloc( LDWORK, 8 )
      ELSE IF ( ITASK.EQ.3 ) THEN
         E = mxCalloc( LDE*NM, 8 )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( A ), N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL( B ), N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL( C ), P*N )
      IF ( WITHD )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), %VAL( D ), P*M )
      IF ( WITHE )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), %VAL( E ), N*N )
C
C Do the actual computations.
C
      IF ( ITASK.EQ.1 ) THEN
         IF ( M.GT.P .AND. P.GT.0 ) THEN
            DW = mxCalloc( P*N, 8 )
            CALL mxCopyPtrToReal8( C,  %VAL( DW ), P*N )
            CALL DLACPY( 'Full', P, N, %VAL( DW ), P, %VAL( C ), LDC )
            CALL mxFree( DW )
C
            IF ( WITHD ) THEN
               DW = mxCalloc( P*M, 8 )
               CALL mxCopyPtrToReal8( D,  %VAL( DW ), P*M )
               CALL DLACPY( 'Full', P, M, %VAL( DW ), P, %VAL( D ),
     $                     LDD )
               CALL mxFree( DW )
            END IF
         END IF
         CALL AB07MD( JOBD, N, M, P, %VAL( A ), LDA, %VAL( B ), LDB,
     $                %VAL( C ), LDC, %VAL( D ), LDD, INFO )
         IF ( M.GT.0 .AND. M.LT.P ) THEN
            IF ( N.GT.0 )
     $         CALL DLACPY( 'Full', M, N, %VAL( C ), LDC, %VAL( C ), M )
            IF ( WITHD )
     $         CALL DLACPY( 'Full', M, P, %VAL( D ), LDD, %VAL( D ), M )
         END IF
      ELSE IF ( ITASK.EQ.2 ) THEN
         CALL AB07ND( N, M, %VAL( A ), LDA, %VAL( B ), LDB, %VAL( C ),
     $                LDC, %VAL( D ), LDD, RCOND, IWORK, %VAL( DWORK ),
     $                LDWORK, INFO )
      ELSE IF ( ITASK.EQ.3 ) THEN
         IF ( M.GT.0 ) THEN
            DW = mxCalloc( N*MAX( N, M ), 8 )
            CALL mxCopyPtrToReal8( A,  %VAL( DW ), N*N )
            CALL DLACPY( 'Full', N, N, %VAL( DW ), N1, %VAL( A ), LDA )
            CALL mxCopyPtrToReal8( E,  %VAL( DW ), N*N )
            CALL DLACPY( 'Full', N, N, %VAL( DW ), N1, %VAL( E ), LDE )
            CALL mxCopyPtrToReal8( B,  %VAL( DW ), N*M )
            CALL DLACPY( 'Full', N, M, %VAL( DW ), N1, %VAL( B ), LDB )
            CALL mxFree( DW )
         END IF
         CALL AG07BD( JOBE, N, M, %VAL( A ), LDA, %VAL( E ), LDE,
     $                %VAL( B ), LDB, %VAL( C ), LDC, %VAL( D ), LDD,
     $                %VAL( A ), LDA, %VAL( E ), LDE, %VAL( B ), LDB,
     $                %VAL( C ), LDC, %VAL( D ), LDD, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      IP = 2
      IF ( INFO.EQ.0 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( NM, NM, 0 )
         CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS(1) ),
     $                          NM*NM )
         IF ( ITASK.EQ.3 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( NM, NM, 0 )
            CALL mxCopyReal8ToPtr( %VAL( E ), mxGetPr( PLHS(IP) ),
     $                             NM*NM )
            IP = IP + 1
         END IF
         PLHS(IP) = mxCreateDoubleMatrix( NM, P, 0 )
         CALL mxCopyReal8ToPtr( %VAL( B ), mxGetPr( PLHS(IP) ),
     $                          NM*P )
         IP = IP + 1
         PLHS(IP) = mxCreateDoubleMatrix( M, NM, 0 )
         CALL mxCopyReal8ToPtr( %VAL( C ), mxGetPr( PLHS(IP) ),
     $                          M*NM )
         IP = IP + 1
         IF ( NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( M, P, 0 )
            CALL mxCopyReal8ToPtr( %VAL( D ), mxGetPr( PLHS(IP) ),
     $                             M*P )
            IP = IP + 1
            IF ( ITASK.EQ.2 .AND. NLHS.GE.IP ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
               CALL mxCopyReal8ToPtr( RCOND, mxGetPr( PLHS(IP) ), 1 )
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      CALL mxFree( A )
      CALL mxFree( B )
      CALL mxFree( C )
      CALL mxFree( D )
      IF ( ITASK.EQ.2 ) THEN
         DEALLOCATE ( IWORK )
         CALL mxFree( DWORK )
      ELSE IF ( ITASK.EQ.3 ) THEN
         CALL mxFree( E )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( ITASK.EQ.2 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB07ND'')' )
     $         INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of INVERT ***
      END
