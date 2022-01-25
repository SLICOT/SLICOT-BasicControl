#include "fintrf.h"
C FSTOEP.F - Gateway function for factoring symmetric positive definite
C            (block) Toeplitz matrices and/or solving associated linear
C            systems using SLICOT routines MB02CD, MB02DD, MB02ED.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C                  [R(,X)] = fstoep(task,T(,B))
C             [G(,Li,R,X)] = fstoep(task,T(,B))
C          [G(,Li),R,H,CS] = fstoep(task,T)
C                [Ru(,Xa)] = fstoep(task,Ta,H,CS,G,R(,Ba))
C  [Gu(,Liu,Ru,Hu,CSu,Xa)] = fstoep(task,Ta,H,CS,G,R(,Li,Ba))
C                      [X] = fstoep(task,T,B)
C
C   task =  1 :                 [R] = fstoep( 1,T)
C                             [R,X] = fstoep( 1,T,B)
C   task =  2 :                 [G] = fstoep( 2,T)
C   task =  3 :               [G,R] = fstoep( 3,T)
C                           [G,R,X] = fstoep( 3,T,B)
C   task =  4 :              [G,Li] = fstoep( 4,T)
C                          [G,Li,X] = fstoep( 4,T,B)
C   task =  5 :            [G,Li,R] = fstoep( 5,T)
C                        [G,Li,R,X] = fstoep( 5,T,B)
C   task =  6 :          [G,R,H,CS] = fstoep( 6,T)
C   task =  7 :       [G,Li,R,H,CS] = fstoep( 7,T)
C   task =  8 :                [Ru] = fstoep( 8,Ta,H,CS,G,R)
C                           [Ru,Xa] = fstoep( 8,Ta,H,CS,G,R,Ba)
C   task =  9 :             [Gu,Ru] = fstoep( 9,Ta,H,CS,G,R)
C                        [Gu,Ru,Xa] = fstoep( 9,Ta,H,CS,G,R,Ba)
C   task = 10 :      [Gu,Ru,Hu,CSu] = fstoep(10,Ta,H,CS,G,R)
C                [Gu,Liu,Ru,Hu,CSu] = fstoep(10,Ta,H,CS,G,R,Li)
C             [Gu,Liu,Ru,Hu,CSu,Xa] = fstoep(10,Ta,H,CS,G,R,Li,Ba)
C   task = 11 :                 [X] = fstoep(11,T,B)
C
C Purpose:
C   To factor a symmetric positive definite block Toeplitz matrix BT
C   and/or its inverse, and/or to compute/update the generator of its
C   inverse, given the first block row / column T of BT, and/or to solve
C   the associated linear systems X*BT = B / BT*X = B.
C   For task <= 5 or task = 11, the first block row / column of BT
C   contains T and has the dimension K-by-N*K / N*K-by-K, while for
C   task >= 6 and task <= 10, the first block row / column of the block
C   Toeplitz matrix, denoted BTA, contains [ T Ta ] / [ T  ] and has the
C   dimension K-by-(N+M)*K / (N+M)*K-by-K.            [ Ta ]
C
C   task =  1: Compute R, the Cholesky factor of BT, i.e., R is
C              upper / lower triangular, so that R'*R = BT / R*R' = BT;
C
C   task =  2: Compute the generator G of inv(BT), of dimension
C              2*K-by-N*K / N*K-by-2*K;
C
C   task =  3: Compute both G and R;
C
C   task =  4: Compute both G and Li, where Li is the Cholesky factor
C              of inv(BT), i.e., Li is lower / upper triangular, so that
C              Li'*Li = inv(BT) / Li*Li' = inv(BT);
C
C   task =  5: Compute both G, Li, and R;
C
C   task =  6: Compute G and R, and deliver the information needed to
C              update a Cholesky factorization.
C
C   task =  7: Compute G, Li, and R, and deliver the information needed
C              to update a Cholesky factorization.
C
C   task =  8: Compute updated R, Ru, given additional blocks of
C              data, Ta, and previous factorization results.
C
C   task =  9: Compute updated G and R, Gu and Ru, given additional
C              blocks of data, Ta, and previous factorization results.
C
C   task = 10: Compute updated G, (Li,) and R, as well as details of the
C              transformations used (needed for subsequent updating),
C              given additional blocks of data, Ta, and previous
C              factorization results.
C
C   task = 11: Solve BT*X = B or X*BT = B.
C
C   Note:      The linear systems BT*X = B or X*BT = B can also be
C              solved when task <= 5, task <> 2, by specifying an
C              additional input and output parameter.
C              Moreover, when task = 8, 9, or 10, one can solve the
C              linear systems BTA*Xa = Ba or Xa*BTA = Ba, by specifying
C              an additional input and output parameter.
C
C Input parameters:
C   task  - integer option to determine the computation to perform:
C           =  1 : compute the Cholesky factor R of the matrix BT;
C           =  2 : compute the generator G of inv(BT);
C           =  3 : compute both G and R;
C           =  4 : compute both G and Li, where Li is the Cholesky
C                  factor of inv(BT);
C           =  5 : compute both G, Li and R;
C           =  6 : compute G and R, and deliver the information needed
C                  to update the factorization;
C           =  7 : compute G, Li, and R, and deliver the information
C                  needed to update the factorization.
C           =  8 : compute updated R, given additional blocks of data;
C           =  9 : compute updated G and R, given additional data;
C           = 10 : compute updated G, (Li,) and R, as well as details of
C                  the transformations used, given additional data;
C           = 11 : solve X*BT = B or BT*X = B.
C   T     - real K-by-N*K / N*K-by-K matrix, containing (part) of the
C           first block row / column of the symmetric positive definite
C           block Toeplitz matrix BT (BTA), with K-by-K blocks.
C   B     - real right hand side NRHS-by-N*K / N*K-by-NRHS matrix.
C   Ta    - real K-by-M*K / M*K-by-K matrix containing additional data
C           of the first block row / column of the symmetric positive
C           definite block Toeplitz matrix BTA.
C   H     - real K-by-N*K / N*K-by-K matrix, containing part of the
C           information about the Householder transformations used.
C   CS    - real 3*(N-1)*K vector containing further details on the
C           Householder transformations and hyperbolic rotations used.
C   G     - real 2*K-by-N*K / N*K-by-2*K generator G of inv(BT).
C   R     - real N*K-by-N*K Cholesky factor of BT.
C   Li    - real N*K-by-N*K Cholesky factor of inv(BT).
C   Ba    - real right hand side NRHS-by-(N+M)*K / (N+M)*K-by-NRHS
C           matrix.
C
C Output parameters:
C   R     - the N*K-by-N*K Cholesky factor of BT.
C   G     - the 2*K-by-N*K / N*K-by-2*K generator G of inv(BT).
C   Li    - the N*K-by-N*K Cholesky factor of inv(BT).
C   X     - the NRHS-by-N*K / N*K-by-NRHS solution matrix of the system
C           X*BT = B / BT*X = B.
C   H     - the K-by-N*K / N*K-by-K matrix, containing part of the
C           information about the Householder transformations used.
C   CS    - the 3*(N-1)*K vector containing further details on the
C           Householder transformations and hyperbolic rotations used.
C   Ru    - the (N+M)*K-by-(N+M)*K Cholesky factor of BTA.
C   Gu    - the 2*K-by-(N+M)*K / (N+M)*K-by-2*K generator of inv(BTA).
C   Liu   - the (N+M)*K-by-(N+M)*K Cholesky factor of inv(BTA).
C   Hu    - the K-by-(N+M)*K / (N+M)*K-by-K matrix containing part of
C           the information about the Householder transformations used.
C   CSu   - the 3*(N+M-1)*K vector containing further details on the
C           Householder transformations and hyperbolic rotations used.
C   Xa    - the NRHS-by-(N+M)*K / (N+M)*K-by-NRHS solution matrix of the
C           system Xa*BTA = Ba / BTA*Xa = Ba.
C
C Comments:
C   1. If the second dimension of T (Ta) is larger than the first one,
C      but it is not a multiple of its first dimension, then the number
C      of columns considered is the largest possible multiple of K.
C      If the first dimension of T (Ta) is larger than the second one,
C      but it is not a multiple of its second dimension, then the number
C      of rows considered is the largest possible multiple of K.
C   2. If one dimension of T is zero, then both K and N are set to 0.
C   3. If task >= 6 and task <= 10, T and Ta must have the same shape.
C   4. If task = 6, 7, or 10, on output, G or Gu contain information
C      needed to perform (further) updates. To get the true generator,
C      the submatrix (K+1:2*K,1:K) / (1:K,K+1:2*K) of G or Gu must be
C      set to zeros(K).
C   5. If N = 1, then the returned results satisfy R'*R = BT
C      (or Ru'*Ru = BTA), Li'*Li = inv(BT) (or Liu'*Liu = inv(BTA)),
C      and G has the size 2*K-by-K (or Gu has the size 2*K-by-(1+M)*K).
C   6. If only the solution of the system X*BT = B / BT*X = B is
C      desired, then task = 11 should be used, for efficiency.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2000.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000,
C   Nov. 2005, Apr. 2009, June 2012, July 2012, Aug. 2012.
C
C **********************************************************************
C
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
      CHARACTER         JOB, TYPET
      INTEGER           INFO, K, LDB, LCS, LDG, LDL, LDR, LDT, LDTA,
     $                  LDWORK, M, N, NRH
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      mwPointer         B, CS, DWORK, G, L, R, T, TA, WL, WR
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           GETX, ROWWIS, UPDATE, UPL
      INTEGER           IP, IPB, ISIZE, K1, KA, MK, NB, NC, NCA, NCB,
     $                  NCG, NCL, NCR, NG, NK, NK1, NKA, NLM, NR, NRA,
     $                  NRG, TASK
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          DLACPY, DTRMM, DTRSM, MA02KD, MA02LD, MB02CD,
     $                  MB02DD, MB02ED
C
C ..Intrinsic functions..
      INTRINSIC         MAX, MIN, MOD
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'FSTOEP requires at least 2 input arguments' )
      ELSE IF ( NLHS.GT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'FSTOEP requires at most 6 output arguments' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   task
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
      TASK = TEMP
      IF ( TASK.LT.1 .OR. TASK.GT.11 ) THEN
         CALL mexErrMsgTxt
     $        ( 'TASK has 1, 2, ..., 11 the only admissible values')
      END IF
C
      UPL = .FALSE.
      IF ( TASK.EQ.1 ) THEN
         JOB = 'O'
         NLM = 2
      ELSE IF ( TASK.EQ.2 ) THEN
         JOB = 'G'
         NLM = 1
      ELSE IF ( TASK.EQ.3 ) THEN
         JOB = 'R'
         NLM = 3
      ELSE IF ( TASK.EQ.4 ) THEN
         JOB = 'L'
         NLM = 3
      ELSE IF ( TASK.EQ.5 ) THEN
         JOB = 'A'
         NLM = 4
      ELSE IF ( TASK.EQ.6 ) THEN
         JOB = 'R'
         NLM = 4
      ELSE IF ( TASK.EQ.7 ) THEN
         JOB = 'A'
         NLM = 5
      ELSE IF ( TASK.EQ.8 ) THEN
         JOB = 'O'
         NLM = 2
      ELSE IF ( TASK.EQ.9 ) THEN
         JOB = 'R'
         NLM = 3
      ELSE IF ( TASK.EQ.10 ) THEN
         IF ( NRHS.EQ.6 ) THEN
            JOB = 'R'
            NLM = 4
         ELSE IF ( NRHS.EQ.7 ) THEN
            JOB = 'A'
            NLM = 5
            UPL = .TRUE.
         ELSE IF ( NRHS.EQ.8 ) THEN
            JOB = 'A'
            NLM = 6
            UPL = .TRUE.
         ELSE
            CALL mexErrMsgTxt
     $           ( 'FSTOEP requires at least 6 input arguments' )
         END IF
      ELSE IF ( TASK.EQ.11 ) THEN
         NLM = 1
      END IF
C
      IF ( NLHS.GT.NLM ) THEN
         WRITE( TEXT, '('' FSTOEP requires at most '',I4,
     $            '' output arguments'')' ) NLM
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( ( TASK.EQ.2 .OR. TASK.EQ.6 .OR. TASK.EQ.7 ) .AND. NRHS.GT.2 )
     $   THEN
            CALL mexErrMsgTxt
     $           ( 'FSTOEP requires at most 2 input arguments' )
      ELSE IF ( ( TASK.LE.5 .OR. TASK.EQ.11 ) .AND. NRHS.GT.3 ) THEN
            CALL mexErrMsgTxt
     $           ( 'FSTOEP requires at most 3 input arguments' )
      ELSE IF ( ( TASK.EQ.8 .OR. TASK.EQ.9 ) .AND. NRHS.GT.7 ) THEN
            CALL mexErrMsgTxt
     $           ( 'FSTOEP requires at most 7 input arguments' )
      ELSE IF ( NRHS.GT.8 ) THEN
            CALL mexErrMsgTxt
     $           ( 'FSTOEP requires at most 8 input arguments' )
      END IF
C
      IF ( TASK.LE.7 .OR. TASK.EQ.11 ) THEN
         IP = 2
         UPDATE = .FALSE.
      ELSE
         IP = 3
         UPDATE = .TRUE.
         IF ( NRHS.LT.6 ) THEN
            CALL mexErrMsgTxt
     $           ( 'FSTOEP requires at least 6 input arguments' )
         END IF
      END IF
C
C   T (or H) (K-by-NxK or NxK-by-K)
C
      NR = mxGetM( PRHS(IP) )
      NC = mxGetN( PRHS(IP) )
      K  = MIN( NR, NC )
      NK = MAX( NR, NC )
      K1 = MAX( 1, K )
C
      IF ( K.LT.0 ) THEN
         IF ( UPDATE ) THEN
            CALL mexErrMsgTxt( 'H must be a matrix' )
         ELSE
            CALL mexErrMsgTxt( 'T must be a matrix' )
         END IF
      END IF
      IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
         IF ( UPDATE ) THEN
            CALL mexErrMsgTxt( 'H must be a real matrix' )
         ELSE
            CALL mexErrMsgTxt( 'T must be a real matrix' )
         END IF
      END IF
C
      IF ( K.EQ.0 ) THEN
         N = 0
      ELSE
         N = NK / K
      END IF
      NK  = N*K
      NK1 = MAX( 1, NK )
      ROWWIS = K.EQ.NR .OR. K.EQ.0
C
      GETX = ( ( ( TASK.EQ.8 .OR.  TASK.EQ.9 ) .AND. NRHS.EQ.7  ) .OR.
     $       ( ( ( TASK.LE.5 .AND. TASK.NE.2 ) .OR.  TASK.EQ.11 ) .AND.
     $                                               NRHS.EQ.3  ) .OR.
     $        ( TASK.EQ.10 .AND. NRHS.EQ.8 ) ) .AND. NLHS.EQ.NLM
      IF ( GETX ) THEN
C
C   B(NRH-by-N*K or N*K-by-NRH) or Ba(NRH-by-(N+M)*K / (N+M)*K-by-NRH)
C
         IPB = NRHS
         IF ( ROWWIS ) THEN
            NRH = mxGetM( PRHS(IPB) )
            NB  = mxGetN( PRHS(IPB) )
         ELSE
            NB  = mxGetM( PRHS(IPB) )
            NRH = mxGetN( PRHS(IPB) )
         END IF
C
         IF ( TASK.LE.5 ) THEN
            IF ( NB.NE.NK ) THEN
               IF ( ROWWIS ) THEN
                  CALL mexErrMsgTxt
     $                 ( 'B must have the same column size as T' )
               ELSE
                  CALL mexErrMsgTxt
     $                 ( 'B must have the same row size as T' )
               END IF
            END IF
            IF ( NRH.LT.0 ) THEN
               CALL mexErrMsgTxt( 'B must be a matrix' )
            END IF
            IF ( mxIsNumeric( PRHS(IPB) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IPB) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'B must be a real matrix' )
            END IF
         END IF
      END IF
C
      IF ( UPDATE ) THEN
C
C   Ta(K-by-M*K or M*K-by-K)
C
         NRA = mxGetM( PRHS(2) )
         NCA = mxGetN( PRHS(2) )
         KA  = MIN( NRA, NCA )
         NKA = MAX( NRA, NCA )
C
         IF ( KA.LT.0 ) THEN
            CALL mexErrMsgTxt( 'Ta must be a matrix' )
         END IF
         IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(2) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'Ta must be a real matrix' )
         END IF
C
         IF ( KA.EQ.0 ) THEN
            M = 0
         ELSE
            M = NKA / KA
         END IF
         MK = M*K
C
         IF ( K.GT.0 .AND. ( ( ROWWIS .AND. NRA.NE.NR ) .OR.
     $                  ( .NOT.ROWWIS .AND. NCA.NE.NC ) ) ) THEN
            CALL mexErrMsgTxt( 'Ta must have the same shape as H' )
         END IF
C
         IF ( GETX ) THEN
            IF ( NB.NE.NK+MK ) THEN
               IF ( ROWWIS ) THEN
                  CALL mexErrMsgTxt
     $                 ( 'Ba must have the same column size as [T Ta]' )
               ELSE
                  CALL mexErrMsgTxt
     $                 ( 'Ba must have the same row size as [T; Ta]' )
               END IF
            END IF
            IF ( NRH.LT.0 ) THEN
               CALL mexErrMsgTxt( 'Ba must be a matrix' )
            END IF
            IF ( mxIsNumeric( PRHS(IPB) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IPB) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'Ba must be a real matrix' )
            END IF
         END IF
C
C   CS(3*(N-1)*K)
C
         ISIZE = mxGetM( PRHS(4) )*mxGetN( PRHS(4) )
C
         IF ( ISIZE.NE.MAX( 0, 3*(N-1)*K ) ) THEN
            WRITE( TEXT, '('' CS must have '',I7,'' elements'')' )
     $             MAX( 0, 3*( N - 1 )*K )
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(4) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'CS must be a real vector' )
         END IF
C
C   G(2*K-by-N*K or N*K-by-2*K)
C
         IF ( ROWWIS ) THEN
            NRG = 2*K
            NCG = NK
         ELSE
            NRG = NK
            NCG = 2*K
         END IF
C
         IF ( mxGetM( PRHS(5) ).NE.NRG ) THEN
            WRITE( TEXT, '('' G must have '',I7,'' rows'')' ) NRG
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxGetN( PRHS(5) ).NE.NCG ) THEN
            WRITE( TEXT, '('' G must have '',I7,'' columns'')' ) NCG
            CALL mexErrMsgTxt( TEXT )
         END IF
C
         IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(5) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'G must be a real matrix' )
         END IF
C
C   R(N*K-by-N*K)
C
         IF ( mxGetM( PRHS(6) ).NE.NK .OR.
     $        mxGetN( PRHS(6) ).NE.NK ) THEN
            WRITE( TEXT,
     $             '('' R must have '',I7,'' rows and columns'')' ) NK
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'R must be a real matrix' )
         END IF
C
         IF ( NRHS.EQ.7 .AND. TASK.EQ.10 ) THEN
C
C   Li(N*K-by-N*K)
C
            IF ( mxGetM( PRHS(7) ).NE.NK .OR.
     $           mxGetN( PRHS(7) ).NE.NK ) THEN
               WRITE( TEXT,
     $              '('' Li must have '',I7,'' rows and columns'')' ) NK
               CALL mexErrMsgTxt( TEXT )
            END IF
            IF ( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(7) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'Li must be a real matrix' )
            END IF
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms
C in DGELQF, DGEQRF, DORMLQ, and DORMQR.
C
      LDWORK = 1
      IF ( ROWWIS ) THEN
         TYPET = 'R'
         LDT = K1
         IF ( UPDATE )
     $      LDTA = LDT
         IF ( TASK.GT.1 .AND. TASK.LE.10 ) THEN
            NRG = 2*K
            LDG = MAX( 1, NRG )
            IF ( TASK.GE.9 ) THEN
               NCG = NK + MK
            ELSE
               NCG = NK
            END IF
            NG = 2*K
         END IF
      ELSE
         TYPET = 'C'
         LDT = MAX( 1, NR )
         IF ( UPDATE )
     $      LDTA = MAX( 1, MK )
         IF ( TASK.GT.1 .AND. TASK.LE.10 ) THEN
            IF ( TASK.GE.9 ) THEN
               NRG = NK + MK
            ELSE
               NRG = NK
            END IF
            LDG = MAX( 1, NRG )
            NCG = 2*K
            NG  = NK
         END IF
      END IF
C
      IF ( TASK.LE.7 ) THEN
         IF ( MOD( TASK, 2 ).EQ.1 .OR. TASK.GE.6 ) THEN
            IF ( TASK.EQ.1 ) THEN
               LDG = 1
               NCG = 0
            END IF
            LDR = MAX( 1, NK )
            NCR = NK
         ELSE
            LDR = 1
            NCR = 0
         END IF
         IF ( TASK.GE.4 .AND. TASK.NE.6 ) THEN
            LDL = MAX( 1, NK )
            NCL = NK
         ELSE
            LDL = 1
            NCL = 0
         END IF
         LCS    = 3*( N - 1 )*K
         LDWORK = MAX( 1, ( N - 1 )*K )
C
      ELSE IF ( UPDATE ) THEN
         NCR = NK + MK
         LDR = MAX( 1, NCR )
         IF ( TASK.EQ.10 ) THEN
            LDL = LDR
            NCL = NCR
         ELSE
            LDL = 1
            NCL = 0
         END IF
         LCS    = 3*( N + M - 1 )*K
         LDWORK = MAX( 1, ( N + M - 1 )*K )
      END IF
C
      IF ( GETX ) THEN
         IF ( ROWWIS ) THEN
            LDB = MAX( 1, NRH )
            NCB = NB
         ELSE
            LDB = MAX( 1, NB )
            NCB = NRH
         END IF
         LDWORK = MAX( LDWORK, N*K*K + ( N + 2 )*K )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      DWORK = mxCalloc( LDWORK,  8 )
      T     = mxCalloc( LDT*NC,  8 )
      IF ( TASK.LE.10 ) THEN
         CS = mxCalloc(     LCS, 8 )
         G  = mxCalloc( LDG*NCG, 8 )
         L  = mxCalloc( LDL*NCL, 8 )
         R  = mxCalloc( LDR*NCR, 8 )
      END IF
      IF ( UPDATE )
     $   TA = mxCalloc( LDTA*NCA, 8 )
      IF ( GETX )
     $   B  = mxCalloc(  LDB*NCB, 8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(T), NR*NC )
      IF ( GETX )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(IPB) ), %VAL(B),
     $                          NRH*NB )
      IF ( UPDATE ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL(TA),
     $                          NRA*NCA )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL(CS),
     $                          ISIZE )
         WR = mxCalloc( MAX( 2*K*NK, NK*NK ), 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), %VAL(WR),
     $                          2*K*NK )
         IF ( NG.GT.0 )
     $      CALL DLACPY( 'Full', NG, 2*K*NK/NG, %VAL(WR), NG, %VAL(G),
     $                   LDG )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), %VAL(WR),
     $                          NK*NK )
         CALL DLACPY( 'Full', NK, NK, %VAL(WR), NK1, %VAL(R), LDR )
         IF ( TASK.EQ.10 .AND. NRHS.GE.7 ) THEN
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), %VAL(WR),
     $                            NK*NK )
            CALL DLACPY( 'Full', NK, NK, %VAL(WR), NK1, %VAL(L), LDL )
         END IF
         CALL mxFree( WR )
      END IF
C
C Do the actual computations.
C
      IF ( TASK.LE.7 ) THEN
         CALL MB02CD( JOB, TYPET, K, N, %VAL(T), LDT, %VAL(G), LDG,
     $                %VAL(R), LDR, %VAL(L), LDL, %VAL(CS), LCS,
     $                %VAL(DWORK), LDWORK, INFO )
C
      ELSE IF ( UPDATE .AND. M.GT.0 ) THEN
         IF ( ROWWIS ) THEN
            WL = mxCalloc( MK*NCL, 8 )
            WR = mxCalloc( LDR*MK, 8 )
            CALL MA02KD( 'Full', NK, K, 1, (N-1)*K+1, K+1, 1, %VAL(R),
     $                   LDR, %VAL(WR), LDR )
            CALL MB02DD( JOB, TYPET, K, M, N, %VAL(TA), LDTA, %VAL(T),
     $                   LDT, %VAL(G), LDG, %VAL(WR), LDR, %VAL(WL), MK,
     $                   %VAL(CS), LCS, %VAL(DWORK), LDWORK, INFO )
            CALL MA02KD( 'Full', NK+MK, MK, 1, 1, 1, NK+1, %VAL(WR),
     $                   LDR, %VAL(R), LDR )
            IF( UPL )
     $         CALL MA02KD( 'Full', MK, NK+MK, 1, 1, NK+1, 1, %VAL(WL),
     $                      MK, %VAL(L), LDL )
         ELSE
            WL = mxCalloc( LDL*MK, 8 )
            WR = mxCalloc( MK*NCR, 8 )
            CALL MA02KD( 'Full', K, NK, (N-1)*K+1, 1, 1, K+1, %VAL(R),
     $                   LDR, %VAL(WR), MK )
            CALL MB02DD( JOB, TYPET, K, M, N, %VAL(TA), LDTA, %VAL(T),
     $                   LDT, %VAL(G), LDG, %VAL(WR), MK, %VAL(WL), LDL,
     $                   %VAL(CS), LCS, %VAL(DWORK), LDWORK, INFO )
            CALL MA02KD( 'Full', MK, NK+MK, 1, 1, NK+1, 1, %VAL(WR), MK,
     $                   %VAL(R), LDR )
            IF( UPL )
     $         CALL MA02KD( 'Full', NK+MK, MK, 1, 1, 1, NK+1, %VAL(WL),
     $                      LDL, %VAL(L), LDL )
         END IF
         CALL mxFree( WL )
         CALL mxFree( WR )
      ELSE IF ( M.EQ.0 ) THEN
         INFO = 0
      END IF
C
      IF ( TASK.EQ.11 ) THEN
         CALL MB02ED( TYPET, K, N, NRH, %VAL(T), LDT, %VAL(B), LDB,
     $                %VAL(DWORK), LDWORK, INFO )
C
      ELSE IF ( GETX ) THEN
         IF ( MOD( TASK, 2 ).EQ.1 .OR. TASK.EQ.8 ) THEN
            IF ( ROWWIS ) THEN
               CALL DTRSM( 'Right', 'Upper', 'NonTranspose',
     $                     'NonUnit', NRH, NB, ONE, %VAL(R), LDR,
     $                     %VAL(B), LDB )
               CALL DTRSM( 'Right', 'Upper', 'Transpose', 'NonUnit',
     $                     NRH, NB, ONE, %VAL(R), LDR, %VAL(B), LDB )
            ELSE
               CALL DTRSM( 'Left', 'Lower', 'NonTranspose',
     $                     'NonUnit', NB, NRH, ONE, %VAL(R), LDR,
     $                     %VAL(B), LDB )
               CALL DTRSM( 'Left', 'Lower', 'Transpose',
     $                     'NonUnit', NB, NRH, ONE, %VAL(R), LDR,
     $                     %VAL(B), LDB )
            END IF
         ELSE IF ( TASK.EQ.4 .OR. TASK.EQ.10 ) THEN
            IF ( ROWWIS ) THEN
               CALL DTRMM( 'Right', 'Lower', 'Transpose',
     $                     'NonUnit', NRH, NB, ONE, %VAL(L), LDL,
     $                     %VAL(B), LDB )
               CALL DTRMM( 'Right', 'Lower', 'NonTranspose',
     $                     'NonUnit', NRH, NB, ONE, %VAL(L), LDL,
     $                     %VAL(B), LDB )
            ELSE
               CALL DTRMM( 'Left', 'Upper', 'Transpose', 'NonUnit',
     $                     NB, NRH, ONE, %VAL(L), LDL, %VAL(B), LDB )
               CALL DTRMM( 'Left', 'Upper', 'NonTranspose',
     $                     'NonUnit', NB, NRH, ONE, %VAL(L), LDL,
     $                     %VAL(B), LDB )
            END IF
         END IF
      END IF
C
C Copy output to MATLAB workspace.
C
      IP = 1
      IF ( UPDATE ) THEN
         NK  = NK + MK
         NK1 = MAX( 1, NK )
      END IF
      IF ( INFO.EQ.0 ) THEN
         IF ( NLHS.GE.1 ) THEN
            IF ( TASK.EQ.1 .OR. TASK.EQ.8 ) THEN
               IP = 2
            ELSE IF ( TASK.LE.10 ) THEN
               IF ( NK.GE.1 .AND. TASK.NE.6 .AND. TASK.NE.7 .AND.
     $                            TASK.NE.10 ) THEN
                  IF ( ROWWIS ) THEN
                     CALL MA02LD( 'Full', K, K, K+1, 1, ZERO, ZERO,
     $                            %VAL(G), LDG )
                  ELSE
                     CALL MA02LD( 'Full', K, K, 1, K+1, ZERO, ZERO,
     $                            %VAL(G), LDG )
                  END IF
               END IF
               PLHS(1) = mxCreateDoubleMatrix( NRG, NCG, 0 )
               CALL mxCopyReal8ToPtr( %VAL(G), mxGetPr( PLHS(1) ),
     $                                NRG*NCG )
               IF ( ( TASK.EQ.3 .OR. TASK.EQ.6 .OR. TASK.EQ.9 .OR.
     $              ( TASK.EQ.10 .AND. NRHS.EQ.6 ) ) .AND.
     $              NLHS.GE.2 ) THEN
                  IP = 3
               ELSE IF ( TASK.GE.4 .AND. NLHS.GE.2 ) THEN
                  IF ( ( ROWWIS .OR. N.EQ.1 ) .AND. NK.GT.1 ) THEN
                     CALL MA02LD( 'Upper', NK-1, NK-1, 1, 2, ZERO, ZERO,
     $                            %VAL(L), LDL )
                  ELSE IF ( NK.GT.1 ) THEN
                     CALL MA02LD( 'Lower', NK-1, NK-1, 2, 1, ZERO, ZERO,
     $                            %VAL(L), LDL )
                  END IF
                  PLHS(2) = mxCreateDoubleMatrix( NK, NK, 0 )
                  CALL mxCopyReal8ToPtr( %VAL(L),
     $                                   mxGetPr( PLHS(2) ), NK*NK )
                  IF ( TASK.EQ.4 ) THEN
                     IP = 3
                  ELSE
                     IP = 4
                  END IF
               END IF
            END IF
            IF ( TASK.NE.2 .AND. TASK.NE.4 .AND. TASK.LE.10 .AND.
     $           NLHS.GE.IP-1 ) THEN
               IF ( ( ROWWIS .OR. N.EQ.1 ) .AND. NK.GT.1 ) THEN
                  CALL MA02LD( 'Lower', NK-1, NK-1, 2, 1, ZERO, ZERO,
     $                         %VAL(R), LDR )
               ELSE IF ( NK.GT.1 ) THEN
                  CALL MA02LD( 'Upper', NK-1, NK-1, 1, 2, ZERO, ZERO,
     $                         %VAL(R), LDR )
               END IF
               PLHS(IP-1) = mxCreateDoubleMatrix( NK, NK, 0 )
               CALL mxCopyReal8ToPtr( %VAL(R),
     $                                mxGetPr( PLHS(IP-1) ), NK*NK )
            END IF
            IF ( ( TASK.EQ.6 .OR. TASK.EQ.7 .OR. TASK.EQ.10 )
     $             .AND. NLHS.GE.IP ) THEN
               IF ( TASK.LT.10 ) THEN
                  PLHS(IP) = mxCreateDoubleMatrix( NR, NC, 0 )
                  CALL mxCopyReal8ToPtr( %VAL(T),
     $                                   mxGetPr( PLHS(IP) ),
     $                                   NR*NC )
               ELSE
                  WR = mxCalloc( K*NK, 8 )
                  IF ( ROWWIS ) THEN
                     CALL DLACPY( 'Full', K, N*K, %VAL(T), LDT,
     $                            %VAL(WR), K1 )
                     CALL MA02KD( 'Full', K, MK, 1, 1, 1, N*K+1,
     $                            %VAL(TA), LDTA, %VAL(WR), K1 )
                     PLHS(IP) = mxCreateDoubleMatrix( K, NK, 0 )
                  ELSE
                     CALL DLACPY( 'Full', N*K, K, %VAL(T), LDT,
     $                            %VAL(WR), NK1 )
                     CALL MA02KD( 'Full', MK, K, 1, 1, N*K+1, 1,
     $                            %VAL(TA), LDTA, %VAL(WR), NK1 )
                     PLHS(IP) = mxCreateDoubleMatrix( NK, K, 0 )
                  END IF
                  CALL mxCopyReal8ToPtr ( %VAL(WR),
     $                                    mxGetPr( PLHS(IP) ),
     $                                    NK*K )
                  CALL mxFree( WR )
               END IF
               IP = IP + 1
               IF ( NLHS.GE.IP ) THEN
                  ISIZE = MAX( 0, 3*NK - 3*K )
                  PLHS(IP) = mxCreateDoubleMatrix( ISIZE, 1, 0 )
                  CALL mxCopyReal8ToPtr( %VAL(CS),
     $                                   mxGetPr( PLHS(IP) ),
     $                                   ISIZE )
                  IP = IP + 1
               END IF
            END IF
            IF ( GETX ) THEN
               IF ( ROWWIS ) THEN
                  PLHS(IP) = mxCreateDoubleMatrix( NRH, NB, 0 )
               ELSE
                  PLHS(IP) = mxCreateDoubleMatrix( NB, NRH, 0 )
               END IF
               CALL mxCopyReal8ToPtr( %VAL(B),
     $                                mxGetPr( PLHS(IP) ), NRH*NB )
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      CALL mxFree( DWORK )
      CALL mxFree(     T )
      IF ( TASK.LE.10 ) THEN
         CALL mxFree( CS )
         CALL mxFree(  G )
         CALL mxFree(  L )
         CALL mxFree(  R )
      END IF
      IF ( UPDATE )
     $   CALL mxFree( TA )
      IF ( GETX )
     $   CALL mxFree( B )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( TASK.LE.7 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB02CD'')' )
     $             INFO
         ELSE IF ( UPDATE ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB02DD'')' )
     $             INFO
         ELSE IF ( TASK.EQ.11 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB02ED'')' )
     $             INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of FSTOEP ***
      END
