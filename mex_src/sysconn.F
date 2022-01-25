#include "fintrf.h"
C SYSCONN.F  - Gateway function to compute a state-space model
C              (A,B,C,D) for various inter-connections of two systems,
C              using SLICOT routines AB05MD, AB05ND, AB05OD, AB05PD,
C              and AB05QD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [A,B,C,D] = SYSCONN(task,A1,B1,C1,D1,A2,B2,C2,D2(,uplo)(,alpha))
C
C   [A,B,C,D] = SYSCONN(1,A1,B1,C1,D1,A2,B2,C2,D2(,uplo))
C   [A,B,C,D] = SYSCONN(2,A1,B1,C1,D1,A2,B2,C2,D2(,alpha))
C   [A,B,C,D] = SYSCONN(3,A1,B1,C1,D1,A2,B2,C2,D2(,alpha))
C   [A,B,C,D] = SYSCONN(4,A1,B1,C1,D1,A2,B2,C2,D2(,alpha))
C   [A,B,C,D] = SYSCONN(5,A1,B1,C1,D1,A2,B2,C2,D2)
C
C Purpose:
C   To compute a state-space model (A,B,C,D) for various inter-
C   connections (depending on the value of task) of two systems,
C   given in state-space form, (A1,B1,C1,D1), and (A2,B2,C2,D2).
C   Specifically,
C
C   task = 1: compute cascaded inter-connection;
C   task = 2: compute feedback inter-connection;
C   task = 3: compute rowwise concatenation (parallel inter-connection
C             on outputs, with separate inputs);
C   task = 4: compute the state-space model (A,B,C,D) corresponding to
C             the sum G = G1 + alpha*G2, where G, G1, and G2 are the
C             transfer-function matrices of the corresponding state-
C             space models (A,B,C,D), (A1,B1,C1,D1), and (A2,B2,C2,D2),
C             respectively;
C   task = 5: append two systems in state-space form, (A1,B1,C1,D1) and
C             (A2,B2,C2,D2), with the transfer-function matrices G1
C             and G2, respectively, and obtain the state-space model
C             (A,B,C,D) corresponding to the transfer-function matrix
C
C                           ( G1 0  )
C                       G = (       ) .                              (1)
C                           ( 0  G2 )
C
C Input parameters:
C   task   - integer specifying the computations to be performed.
C            task = 1 :  compute the cascaded inter-connection;
C            task = 2 :  compute the feedback inter-connection;
C            task = 3 :  compute the rowwise  inter-connection;
C            task = 4 :  compute the parallel inter-connection;
C            task = 5 :  compute the compound system in (1).
C   A1     - the n1-by-n1 state dynamics matrix A1.
C   B1     - the n1-by-m1 input/state matrix B1.
C   C1     - the p1-by-n1 state/output matrix C1.
C   D1     - the p1-by-m1 input/output matrix D1.
C   A2     - the n2-by-n2 state dynamics matrix A2.
C   B2     - the n2-by-m2 input/state matrix B2.  If task = 1, m2 = p1.
C                                                 If task = 2, m2 = p1.
C                                                 If task = 4, m2 = m1.
C   C2     - the p2-by-n2 state/output matrix C2. If task = 2, p2 = m1.
C                                                 If task = 3, p2 = p1.
C                                                 If task = 4, p2 = p1.
C   D2     - the p2-by-m2 input/output matrix D2.
C   uplo   - (optional) integer indicating whether the matrix A should
C            be obtained in the upper or lower block diagonal form:
C            uplo = 1 :  Obtain A in the lower block diagonal form;
C            uplo = 2 :  Obtain A in the upper block diagonal form.
C            Default:  uplo = 1.
C   alpha  - (optional) real coefficient multiplying the transfer-
C            function matrix (or the output equation) of the second
C            system. For task = 2,
C            alpha = +1 corresponds to positive feedback, and
C            alpha = -1 corresponds to negative feedback.
C            For task = 3 or task = 4, alpha is not constrained, but
C            alpha = 0 is not dealt with as a special case.
C            Default:  alpha = -1, for task = 2;
C                      alpha =  1, for task = 3 or task = 4.
C
C Output parameters:
C   A      - the n-by-n state dynamics matrix A of the obtained system,
C            where n = n1+n2.
C   B      - the n-by-m input/state matrix B, where m = m1, if task <= 2
C            or task = 4, and m = m1 + m2, if task = 3 or task = 5.
C   C      - the p-by-n state/output matrix C, where p = p2 if task = 1,
C            p = p1 if task = 2, 3, or 4, and p = p1 + p2, if task = 5.
C   D      - the p-by-m input/output matrix D.
C
C Method:
C
C   task = 1: After cascaded inter-connection of the two systems
C
C   X1'     = A1*X1 + B1*U,
C   V       = C1*X1 + D1*U,  of order n1,
C
C   X2'     = A2*X2 + B2*V,
C   Y       = C2*X2 + D2*V,  of order n2,
C
C   where  '  denotes differentiation with respect to time,
C   the following state-space model is obtained:
C
C   X'      = A*X + B*U,
C   Y       = C*X + D*U,                                             (2)
C
C   where
C
C   A = ( A1     0  ),   B = (  B1   ),
C       ( B2*C1  A2 )        ( B2*D1 )
C
C   C = ( D2*C1  C2 ),   D = ( D2*D1 ).
C
C   This form is returned when uplo = 1. When A1 and A2 are block lower
C   triangular, the resulting state matrix is also block lower
C   triangular.
C
C   By applying a similarity transformation to the system above,
C   using the matrix  ( 0  I ),  where I is the identity matrix of
C                     ( J  0 )
C   order  n2,  and  J  is the identity matrix of order n1, the
C   system matrices become
C
C   A = ( A2  B2*C1 ),   B = ( B2*D1 ),
C       ( 0     A1  )        (  B1   )
C
C   C = ( C2  D2*C1 ),   D = ( D2*D1 ).
C
C   This form is returned when uplo = 2. When A1 and A2 are block upper
C   triangular (for instance, in the real Schur form), the resulting
C   state matrix is also block upper triangular.
C
C   task = 2: After feedback inter-connection of the two systems,
C
C   X1'     = A1*X1 + B1*U1,
C   Y1      = C1*X1 + D1*U1,  of order n1,
C
C   X2'     = A2*X2 + B2*U2,
C   Y2      = C2*X2 + D2*U2,  of order n2,
C
C   the state-space model (2) will be obtained, where
C
C   U = U1 + alpha*Y2,    X  =  ( X1 ),
C   Y = Y1 = U2,                ( X2 )
C
C   A = ( A1  -  alpha*B1*E12*D2*C1       -  alpha*B1*E12*C2    ),
C       (        B2*E21*C1            A2  -  alpha*B2*E21*D1*C2 )
C
C   B = (  B1*E12    ),
C       (  B2*E21*D1 )
C
C   C = (  E21*C1     -  alpha*E21*D1*C2 ),
C
C   D = (  E21*D1 ),
C
C   E21  =  inv( I + alpha*D1*D2 ) and
C   E12  =  inv( I + alpha*D2*D1 ) = I - alpha*D2*E21*D1.
C
C   Taking n1 = 0 and/or n2 = 0 on the function call will solve the
C   constant plant and/or constant feedback cases.
C
C   task = 3: After rowwise concatenation (parallel inter-connection
C   with separate inputs) of the two systems,
C
C   X1'     = A1*X1 + B1*U
C   Y1      = C1*X1 + D1*U
C
C   X2'     = A2*X2 + B2*V
C   Y2      = C2*X2 + D2*V
C
C   (where  '  denotes differentiation with respect to time),
C
C   with the output equation for the second system multiplied by a
C   scalar alpha, the following state-space model will be obtained:
C
C   X'      = A*X + B*(U),
C                     (V)
C
C   Y       = C*X + D*(U),
C                     (V)
C
C   where
C
C   A = ( A1   0  ),         B = ( B1   0  ),
C       ( 0    A2 )              ( 0    B2 )
C
C   C = ( C1   alpha*C2 ),   D = ( D1   alpha*D2 ).
C
C   task = 4: The matrices of the resulting systems are determined as:
C
C       ( A1   0  )             ( B1 )
C   A = (         ) ,       B = (    ) ,
C       ( 0    A2 )             ( B2 )
C
C   C = ( C1  alpha*C2 ) ,  D = D1 + alpha*D2 .
C
C   task = 5: The matrices of the resulting systems are determined as:
C
C       ( A1   0  )         ( B1  0  )
C   A = (         ) ,   B = (        ) ,
C       ( 0    A2 )         ( 0   B2 )
C
C       ( C1   0  )         ( D1  0  )
C   C = (         ) ,   D = (        ) .
C       ( 0    C2 )         ( 0   D2 )
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, July 2003.
C
C Revisions:
C   V. Sima, Feb. 2004, Apr. 2009, Dec. 2010, Sep. 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
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
      CHARACTER         OVER, UPLO
      INTEGER           INFO, LDA, LDA1, LDA2, LDB, LDB1, LDB2, LDC,
     $                  LDC1, LDC2, LDD, LDD1, LDD2, LDWORK, M, M1, M2,
     $                  N, N1, N2, P, P1, P2
      DOUBLE PRECISION  ALPHA
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      mwPointer         A1, A2, B1, B2, C1, C2, D1, D2, DW, DWORK
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           LOVER, USE1
      INTEGER           IP, ITASK, ITMP, J, MB1, MB2, MD1, MD2, MX, NA1,
     $                  NA2, NC1, NC2
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          AB05MD, AB05ND, AB05OD, AB05PD, AB05QD, DLACPY,
     $                  MA02KV
C .. Intrinsic functions ..
      INTRINSIC         ABS, MAX
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.9 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSCONN requires at least 9 input arguments' )
      ELSE IF ( NLHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSCONN requires at least 4 output arguments' )
      END IF
C
C   task, A1(n1xn1), B1(n1,m1), C1(p1,n1), D1(p1,m1),
C         A2(n2xn2), B2(n2,m2), C2(p2,n2), D2(p2,m2)(, uplo)(, alpha)
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
      IF ( ITASK.LT.1 .OR. ITASK.GT.5 ) THEN
         CALL mexErrMsgTxt
     $     ( 'TASK has 1, 2, 3, 4, or 5 the only admissible values' )
      END IF
C
      OVER  = 'O'
      LOVER = .TRUE.
C
      N1 = mxGetM( PRHS(2) )
      M1 = mxGetN( PRHS(3) )
      P1 = mxGetM( PRHS(4) )
      N2 = mxGetM( PRHS(6) )
      M2 = mxGetN( PRHS(7) )
      P2 = mxGetM( PRHS(8) )
C
      IF ( ITASK.LE.2 ) THEN
         IF ( M2.NE.P1 ) THEN
            CALL mexErrMsgTxt(
     $         'The columnsize of B2 must equal the rowsize of C1' )
         ELSE IF ( ITASK.EQ.2 .AND. P2.NE.M1 ) THEN
            CALL mexErrMsgTxt(
     $         'The rowsize of C2 must equal the columnsize of B1' )
         END IF
C
      ELSE IF ( ITASK.LE.4 ) THEN
         IF ( P2.NE.P1 ) THEN
            CALL mexErrMsgTxt(
     $         'The rowsize of C2 must equal the rowsize of C1' )
         ELSE IF ( ITASK.EQ.4 .AND. M2.NE.M1 ) THEN
            CALL mexErrMsgTxt(
     $         'The columnsize of B2 must equal the columnsize of B1' )
         END IF
      END IF
C
      IF ( mxGetN( PRHS(2) ).NE.N1 ) THEN
         CALL mexErrMsgTxt( 'A1 must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A1 must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(3) ).NE.N1 ) THEN
         CALL mexErrMsgTxt( 'B1 must have the same rowsize as A1' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B1 must be a real matrix' )
      END IF
C
      IF ( mxGetN( PRHS(4) ).NE.N1 ) THEN
         CALL mexErrMsgTxt( 'C1 must have the same columnsize as A1')
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C1 must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(5) ).NE.P1 ) THEN
         CALL mexErrMsgTxt( 'D1 must have the same rowsize as C1' )
      END IF
      IF ( mxGetN( PRHS(5) ).NE.M1 ) THEN
         CALL mexErrMsgTxt( 'D1 must have the same columnsize as B1')
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D1 must be a real matrix' )
      END IF
C
      IF ( mxGetN( PRHS(6) ).NE.N2 ) THEN
         CALL mexErrMsgTxt( 'A2 must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(6) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A2 must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(7) ).NE.N2 ) THEN
         CALL mexErrMsgTxt( 'B2 must have the same rowsize as A2' )
      END IF
      IF ( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(7) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B2 must be a real matrix' )
      END IF
C
      IF ( mxGetN( PRHS(8) ).NE.N2 ) THEN
         CALL mexErrMsgTxt( 'C2 must have the same columnsize as A2')
      END IF
      IF ( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(8) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C2 must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(9) ).NE.P2 ) THEN
         CALL mexErrMsgTxt( 'D2 must have the same rowsize as C2' )
      END IF
      IF ( mxGetN( PRHS(9) ).NE.M2 ) THEN
         CALL mexErrMsgTxt( 'D2 must have the same columnsize as B2')
      END IF
      IF ( mxIsNumeric( PRHS(9) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(9) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D2 must be a real matrix' )
      END IF
C
C   uplo/alpha
C
      IP = 10
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            IF ( ITASK.EQ.1 ) THEN
               CALL mexErrMsgTxt( 'UPLO must be a scalar' )
            ELSE IF ( ITASK.LT.5 ) THEN
               CALL mexErrMsgTxt( 'ALPHA must be a scalar' )
            END IF
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            IF ( ITASK.EQ.1 ) THEN
               CALL mexErrMsgTxt( 'UPLO must be an integer scalar' )
            ELSE IF ( ITASK.LT.5 ) THEN
               CALL mexErrMsgTxt( 'ALPHA must be a real scalar' )
            END IF
         END IF
         IF ( ITASK.EQ.1 ) THEN
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            ITMP = TEMP
            IF ( ITMP.LT.1 .OR. ITMP.GT.2 ) THEN
               CALL mexErrMsgTxt
     $              ( 'UPLO has 1 or 2 the only admissible values' )
            END IF
            USE1 = ITMP.EQ.1
            IF ( USE1 ) THEN
               UPLO = 'L'
            ELSE
               UPLO = 'U'
            END IF
         ELSE IF ( ITASK.LT.5 ) THEN
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), ALPHA, 1 )
            IF ( ITASK.EQ.2 .AND. ABS( ALPHA ).NE.ONE ) THEN
               CALL mexErrMsgTxt( 'ALPHA must be +1 or -1' )
            END IF
         END IF
      ELSE
         IF ( ITASK.EQ.1 ) THEN
            USE1 = .TRUE.
            UPLO = 'L'
         ELSE IF ( ITASK.EQ.2 ) THEN
            ALPHA = -ONE
         ELSE IF ( ITASK.LT.5 ) THEN
            ALPHA = ONE
         END IF
      END IF
C
      IF ( ITASK.GT.1 )
     $   USE1 = .TRUE.
C
C Determine the lenghts of working arrays.
C
      N = N1 + N2
      IF ( ITASK.LE.2 .OR. ITASK.EQ.4 ) THEN
         M = M1
      ELSE
         M = M1 + M2
      END IF
      IF ( ITASK.EQ.1 ) THEN
         P = P2
      ELSE IF ( ITASK.NE.5 ) THEN
         P = P1
      ELSE
         P = P1 + P2
      END IF
C
      LDA = MAX( 1, N )
      LDB = LDA
      IF ( ITASK.EQ.1 ) THEN
         IF ( USE1 ) THEN
            LDD = MAX( 1, P1, P2 )
         ELSE
            LDD = MAX( 1, P2 )
         END IF
      ELSE
         LDD = MAX( 1, P )
      END IF
      IF ( N.GT.0 ) THEN
         LDC = LDD
      ELSE
         LDC = 1
      END IF
C
      IF ( USE1 ) THEN
         LDA1 = LDA
         LDB1 = LDB
         LDC1 = LDC
         LDD1 = LDD
         LDA2 = MAX( 1, N2 )
         LDB2 = LDA2
         LDD2 = MAX( 1, P2 )
         IF ( N2.GT.0 ) THEN
            LDC2 = LDD2
         ELSE
            LDC2 = 1
         END IF
         NA1 = N
         NA2 = N2
         NC1 = N
         NC2 = N2
         MB1 = M
         MD1 = M
         MB2 = M2
         MD2 = M2
      ELSE
         LDA2 = LDA
         LDB2 = LDB
         LDC2 = LDC
         LDD2 = LDD
         LDA1 = MAX( 1, N1 )
         LDB1 = LDA1
         LDD1 = MAX( 1, P1 )
         IF ( N1.GT.0 ) THEN
            LDC1 = LDD1
         ELSE
            LDC1 = 1
         END IF
         NA1 = N1
         NA2 = N
         NC1 = N1
         NC2 = N
         MB1 = M1
         MD1 = M1
         MB2 = MAX( M, M2 )
         MD2 = MB2
      END IF
C
C   ldwork
C
      IF ( ITASK.EQ.1 ) THEN
         LDWORK = MAX( 1, P1*MAX( N1, M1, N2, P2 ) )
      ELSE IF ( ITASK.EQ.2 ) THEN
         LDWORK = MAX( P1*P1, M1*M1, N1*P1 )
         IF( LOVER ) THEN
            IF( M1.GT.N*N2 )
     $         LDWORK = MAX( LDWORK, M1*( M1 + 1 ) )
            LDWORK = N1*P1 + LDWORK
         END IF
         LDWORK = MAX( 1, LDWORK )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      A1 = mxCalloc( LDA1*NA1, 8 )
      A2 = mxCalloc( LDA2*NA2, 8 )
      B1 = mxCalloc( LDB1*MB1, 8 )
      B2 = mxCalloc( LDB2*MB2, 8 )
      C1 = mxCalloc( LDC1*NC1, 8 )
      C2 = mxCalloc( LDC2*NC2, 8 )
      D1 = mxCalloc( LDD1*MD1, 8 )
      D2 = mxCalloc( LDD2*MD2, 8 )
      IF ( ITASK.EQ.1 ) THEN
         DWORK = mxCalloc( LDWORK, 8 )
      ELSE IF ( ITASK.EQ.2 ) THEN
         DWORK = mxCalloc( LDWORK, 8 )
         ALLOCATE( IWORK( P1 ) )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( A1 ), N1*N1 )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL( B1 ), N1*M1 )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL( C1 ), P1*N1 )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), %VAL( D1 ), P1*M1 )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), %VAL( A2 ), N2*N2 )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), %VAL( B2 ), N2*M2 )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(8) ), %VAL( C2 ), P2*N2 )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(9) ), %VAL( D2 ), P2*M2 )
C
C Reformat the arrays for overwritting with results, if needed.
C
      IF ( USE1 ) THEN
         IF ( N1.LT.LDA ) THEN
            MX = MAX( 1, N1 )
            DW = mxCalloc( N1*MAX( N1, M1 ), 8 )
            CALL MA02KV( N1*N1, 1, 1, %VAL( A1 ), %VAL( DW ) )
            CALL DLACPY( 'Full', N1, N1, %VAL( DW ), MX, %VAL( A1 ),
     $                   LDA )
            CALL MA02KV( N1*M1, 1, 1, %VAL( B1 ), %VAL( DW ) )
            CALL DLACPY( 'Full', N1, M1, %VAL( DW ), MX, %VAL( B1 ),
     $                   LDB )
            CALL mxFree( DW )
         END IF
C
         IF ( P1.LT.LDC .AND. P1.GT.0 ) THEN
            DW = mxCalloc( P1*N1, 8 )
            CALL MA02KV( P1*N1, 1, 1, %VAL( C1 ), %VAL( DW ) )
            CALL DLACPY( 'Full', P1, N1, %VAL( DW ), MAX( 1, P1 ),
     $                   %VAL( C1 ), LDC )
            CALL mxFree( DW )
         END IF
         IF ( P1.LT.LDD .AND. P1.GT.0 ) THEN
            DW = mxCalloc( P1*M1, 8 )
            CALL MA02KV( P1*M1, 1, 1, %VAL( D1 ), %VAL( DW ) )
            CALL DLACPY( 'Full', P1, M1, %VAL( DW ), MAX( 1, P1 ),
     $                   %VAL( D1 ), LDD )
            CALL mxFree( DW )
         END IF
C
      ELSE
         IF ( N2.LT.LDA ) THEN
            MX = MAX( 1, N2 )
            DW = mxCalloc( N2*MAX( N2, M2 ), 8 )
            CALL MA02KV( N2*N2, 1, 1, %VAL( A2 ), %VAL( DW ) )
            CALL DLACPY( 'Full', N2, N2, %VAL( DW ), MX, %VAL( A2 ),
     $                   LDA )
            CALL MA02KV( N2*M2, 1, 1, %VAL( B2 ), %VAL( DW ) )
            CALL DLACPY( 'Full', N2, M2, %VAL( DW ), MX, %VAL( B2 ),
     $                   LDB )
            CALL mxFree( DW )
         END IF
C
         IF ( P2.LT.LDC .AND. P2.GT.0 ) THEN
            DW = mxCalloc( P2*N2, 8 )
            CALL MA02KV( P2*N2, 1, 1, %VAL( C2 ), %VAL( DW ) )
            CALL DLACPY( 'Full', P2, N2, %VAL( DW ), MAX( 1, P2 ),
     $                   %VAL( C2 ), LDC )
            CALL mxFree( DW )
         END IF
         IF ( P2.LT.LDD .AND. P2.GT.0 ) THEN
            DW = mxCalloc( P2*M2, 8 )
            CALL MA02KV( P2*M2, 1, 1, %VAL( D2 ), %VAL( DW ) )
            CALL DLACPY( 'Full', P2, M2, %VAL( DW ), MAX( 1, P2 ),
     $                   %VAL( D2 ), LDD )
            CALL mxFree( DW )
         END IF
      END IF
C
C Do the actual computations.
C
      IF ( ITASK.EQ.1 ) THEN
         IF ( USE1 ) THEN
            CALL AB05MD( UPLO, OVER, N1, M1, P1, N2, P2, %VAL( A1 ),
     $                   LDA1, %VAL( B1 ), LDB1, %VAL( C1 ), LDC1,
     $                   %VAL( D1 ), LDD1, %VAL( A2 ), LDA2, %VAL( B2 ),
     $                   LDB2, %VAL( C2 ), LDC2, %VAL( D2 ), LDD2, N,
     $                   %VAL( A1 ), LDA, %VAL( B1 ), LDB, %VAL( C1 ),
     $                   LDC, %VAL( D1 ), LDD, %VAL( DWORK ), LDWORK,
     $                   INFO )
            IF ( LDD.GT.P2 .AND. P2.GT.0 ) THEN
               CALL DLACPY( 'Full', P, N, %VAL( C1 ), LDC, %VAL( C1 ),
     $                      P2 )
               CALL DLACPY( 'Full', P, M, %VAL( D1 ), LDD, %VAL( D1 ),
     $                      P2 )
            END IF
         ELSE
            CALL AB05MD( UPLO, OVER, N1, M1, P1, N2, P2, %VAL( A1 ),
     $                   LDA1, %VAL( B1 ), LDB1, %VAL( C1 ), LDC1,
     $                   %VAL( D1 ), LDD1, %VAL( A2 ), LDA2, %VAL( B2 ),
     $                   LDB2, %VAL( C2 ), LDC2, %VAL( D2 ), LDD2, N,
     $                   %VAL( A2 ), LDA, %VAL( B2 ), LDB, %VAL( C2 ),
     $                   LDC, %VAL( D2 ), LDD, %VAL( DWORK ), LDWORK,
     $                   INFO )
         END IF
      ELSE IF ( ITASK.EQ.2 ) THEN
         CALL AB05ND( OVER, N1, M1, P1, N2, ALPHA, %VAL( A1 ), LDA1,
     $                %VAL( B1 ), LDB1, %VAL( C1 ), LDC1, %VAL( D1 ),
     $                LDD1, %VAL( A2 ), LDA2, %VAL( B2 ), LDB2,
     $                %VAL( C2 ), LDC2, %VAL( D2 ), LDD2, N, %VAL( A1 ),
     $                LDA, %VAL( B1 ), LDB, %VAL( C1 ), LDC, %VAL( D1 ),
     $                LDD, IWORK, %VAL( DWORK ), LDWORK, INFO )
      ELSE IF ( ITASK.EQ.3 ) THEN
         CALL AB05OD( OVER, N1, M1, P1, N2, M2, ALPHA, %VAL( A1 ), LDA1,
     $                %VAL( B1 ), LDB1, %VAL( C1 ), LDC1, %VAL( D1 ),
     $                LDD1, %VAL( A2 ), LDA2, %VAL( B2 ), LDB2,
     $                %VAL( C2 ), LDC2, %VAL( D2 ), LDD2, N, M,
     $                %VAL( A1 ), LDA, %VAL( B1 ), LDB, %VAL( C1 ), LDC,
     $                %VAL( D1 ), LDD, INFO )
      ELSE IF ( ITASK.EQ.4 ) THEN
         CALL AB05PD( OVER, N1, M, P, N2, ALPHA, %VAL( A1 ), LDA1,
     $                %VAL( B1 ), LDB1, %VAL( C1 ), LDC1, %VAL( D1 ),
     $                LDD1, %VAL( A2 ), LDA2, %VAL( B2 ), LDB2,
     $                %VAL( C2 ), LDC2, %VAL( D2 ), LDD2, N, %VAL( A1 ),
     $                LDA, %VAL( B1 ), LDB, %VAL( C1 ), LDC, %VAL( D1 ),
     $                LDD, INFO )
      ELSE IF ( ITASK.EQ.5 ) THEN
         CALL AB05QD( OVER, N1, M1, P1, N2, M2, P2, %VAL( A1 ), LDA1,
     $                %VAL( B1 ), LDB1, %VAL( C1 ), LDC1, %VAL( D1 ),
     $                LDD1, %VAL( A2 ), LDA2, %VAL( B2 ), LDB2,
     $                %VAL( C2 ), LDC2, %VAL( D2 ), LDD2, N, M, P,
     $                %VAL( A1 ), LDA, %VAL( B1 ), LDB, %VAL( C1 ), LDC,
     $                %VAL( D1 ), LDD, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
         PLHS(2) = mxCreateDoubleMatrix( N, M, 0 )
         PLHS(3) = mxCreateDoubleMatrix( P, N, 0 )
         PLHS(4) = mxCreateDoubleMatrix( P, M, 0 )
         IF ( USE1 ) THEN
            CALL mxCopyReal8ToPtr( %VAL( A1 ), mxGetPr( PLHS(1) ),
     $                             N*N )
            CALL mxCopyReal8ToPtr( %VAL( B1 ), mxGetPr( PLHS(2) ),
     $                             N*M )
            CALL mxCopyReal8ToPtr( %VAL( C1 ), mxGetPr( PLHS(3) ),
     $                             P*N )
            CALL mxCopyReal8ToPtr( %VAL( D1 ), mxGetPr( PLHS(4) ),
     $                             P*M )
         ELSE
            CALL mxCopyReal8ToPtr( %VAL( A2 ), mxGetPr( PLHS(1) ),
     $                             N*N )
            CALL mxCopyReal8ToPtr( %VAL( B2 ), mxGetPr( PLHS(2) ),
     $                             N*M )
            CALL mxCopyReal8ToPtr( %VAL( C2 ), mxGetPr( PLHS(3) ),
     $                             P*N )
            CALL mxCopyReal8ToPtr( %VAL( D2 ), mxGetPr( PLHS(4) ),
     $                             P*M )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      CALL mxFree( A1 )
      CALL mxFree( A2 )
      CALL mxFree( B1 )
      CALL mxFree( B2 )
      CALL mxFree( C1 )
      CALL mxFree( C2 )
      CALL mxFree( D1 )
      CALL mxFree( D2 )
      IF ( ITASK.EQ.1 ) THEN
         CALL mxFree( DWORK )
      ELSE IF ( ITASK.EQ.2 ) THEN
         CALL mxFree( DWORK )
         DEALLOCATE ( IWORK )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( ITASK.EQ.2 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB05ND'')' )
     $         INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SYSCONN ***
      END
