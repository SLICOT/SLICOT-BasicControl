#include "fintrf.h"
C HAEIGZ.F - Gateway function for computing a balancing transformation
C            or the eigenvalues of a complex Hamiltonian matrix, using
C            SLICOT Library routines MB04DZ and MB03XZ.
C            The gateway accepts real or complex input matrices.
C
C RELEASE 2.0 of SLICOT Basic Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [W(,Ao,Go(,U1,U2(,l,scal)))]   = HaeigZ(A,QG(,job,jobu,balanc))
C   [W(,Se,De(,Ue1,Ue2(,l,scal)))] = HaeigZ(A,QG(,job,jobu,balanc))
C   [l,scal(,Ab,QGb)]              = HaeigZ(A,QG,job,balanc)
C
C   [l,scal]                 = HaeigZ(A,QG,-1,balanc)
C   [l,scal,Ab,QGb]          = HaeigZ(A,QG,-1,balanc)
C   [W]                      = HaeigZ(A,QG)
C   [W,Se,De]                = HaeigZ(A,QG)
C   [W,Se,De,Ue1,Ue2]        = HaeigZ(A,QG,0,1)
C   [W,Se,De,Ue1,Ue2,l,scal] = HaeigZ(A,QG,0,1,balanc)
C   [W,Ao]                   = HaeigZ(A,QG,1)
C   [W,Ao,U1,U2]             = HaeigZ(A,QG,1,1)
C   [W,Ao,Go]                = HaeigZ(A,QG,2)
C   [W,Ao,Go,U1,U2]          = HaeigZ(A,QG,2,1)
C   [W,Ao,Go,U1,U2,l,scal]   = HaeigZ(A,QG,2,1,balanc)
C
C Purpose:
C   To compute the eigenvalues of a complex n-by-n Hamiltonian matrix H,
C   with
C
C                 [  A   G  ]         H        H
C           H  =  [       H ],   G = G ,  Q = Q ,                 (1)
C                 [  Q  -A  ]
C
C   where A, G and Q are complex n-by-n matrices.
C
C   Due to the structure of H, if lambda is an eigenvalue, then
C   -conjugate(lambda) is also an eigenvalue. This does not mean that
C   purely imaginary eigenvalues are necessarily multiple. The function
C   computes the eigenvalues of H using an embedding to a skew-
C   Hamiltonian matrix He,
C
C                  [  Ae   Ge  ]            T            T
C           He  =  [         T ],   Ge = -Ge ,   Qe = -Qe ,        (2)
C                  [  Qe   Ae  ]
C
C   where Ae, Ge, and Qe are real 2*n-by-2*n matrices. Then, an
C   orthogonal symplectic matrix Ue is used to reduce He to the
C   structured real Schur form,
C
C         T          [  Se   De ]            T
C        Ue He Ue =  [        T ],   De = -De ,                    (3)
C                    [  0    Se ]
C
C   where Ue is a 4n-by-4n real symplectic matrix, and Se is upper
C   quasi-triangular (real Schur form).
C
C   Optionally, if job > 0, the matrix i*He is further transformed to
C   the structured complex Schur form
C
C         H            [  Ao  Go ]          H
C        U (i*He) U =  [       H ],  Go = Go ,                     (4)
C                      [  0  -Ao ]
C
C   where U is a 4n-by-4n unitary symplectic matrix, and Ao is upper
C   triangular (Schur form). Optionally, if jobu = 1, the unitary
C   symplectic transformation matrix
C
C         (  U1  U2  )
C     U = (          )
C         ( -U2  U1  )
C
C   is computed.
C
C   If job = -1, an accurate similarity transformation T such that
C   Hb = T\H*T has, as nearly as possible, approximately equal row and
C   column norms. T is a permutation of a diagonal matrix and
C   symplectic. T is stored in an n-vector scal as described in MB04DZ.
C
C Input parameters:
C   A      - the n-by-n matrix A.
C   QG     - an  n-by-(n+1) matrix containing the triangles of the
C            Hermitian matrices G and Q, as follows:
C            the leading n-by-n lower triangular part contains the
C            lower triangle of the matrix Q, and the n-by-n upper
C            triangular part of the submatrix in the columns 2 to n+1
C            contains the upper triangle of the matrix G of H in (1).
C            So, if i >= j, then Q(i,j) = conj(Q(j,i)) is stored in
C            QG(i,j) and G(j,i) = conj(G(i,j)) is stored in QG(j,i+1).
C            QG is an empty matrix if n = 0.
C   job    - (optional) scalar indicating the computation to be
C            performed, as follows:
C            = -1 :  compute a balancing transformation only;
C            =  0 :  compute the eigenvalues only (default);
C            =  1 :  compute the eigenvalues and the matrix Ao in (4);
C            =  2 :  compute the eigenvalues and the matrices Ao and Go
C                    in (4).
C   jobu   - (optional) if job > 0, scalar indicating whether the
C            unitary transformation matrix U is returned, as follows:
C            = 0 :  U is not required (default);
C            = 1 :  on exit, U contains the unitary transformation
C                   matrix.
C   balanc - determines whether H should be permuted (balanc = 1),
C            scaled (balanc = 2), or permuted and scaled (balanc = 3)
C            prior to eigenvalue computations. Otherwise balanc = 0
C            (default). This parameter is optional if job >= 0, but
C            compulsory if job < 0.
C
C Output parameters:
C   l      - if job = -1 or balanc > 0, an integer determined when H was
C            balanced. The balanced A(I,J) = 0 if I > J and J = 1:l-1.
C            The balanced Q(I,J) = 0 if J = 1:l-1 or I = 1:l-1.
C   scal   - if job = -1 or balanc > 0, an n-vector containing details
C            of the permutation and/or scaling factors applied when
C            balancing. See MB04DZ, for details.
C   Ab     - if job = -1, the matrix A of the balanced Hamiltonian Hb.
C            The lower triangular part of the first l-1 columns of A is
C            zero.
C   QGb    - if job = -1, the matrices Q and G of the balanced
C            Hamiltonian Hb, stored compactly. The lower triangular and
C            diagonal part of the first l-1 columns of QGb is zero.
C   W      - the 2n-vector of the eigenvalues of the matrix H.
C   Se     - if job = 0, the computed 2n-by-2n upper real Schur
C            submatrix Se in (3).
C   De     - if job = 0, the computed 2n-by-2n skew-symmetric submatrix
C            De in (3).
C   Ue1    - if job = 0, the computed 2n-by-2n (1,1) block of the matrix
C            Ue in (3).
C   Ue2    - if job = 0, the computed 2n-by-2n (2,1) block of the matrix
C            Ue in (3).
C   Ao     - if job > 0, the computed 2n-by-2n upper triangular
C            submatrix Ao in (4).
C   Go     - if job = 2, the computed 2n-by-2n Hermitian matrix Go
C            in (4).
C   U1     - if jobu = 1, a 2n-by-2n matrix containing the computed
C            matrix U1.
C   U2     - if jobu = 1, an 2n-by-2n matrix containing the computed
C            matrix U2.
C
C Contributors:
C   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2011.
C
C Revisions:
C   V. Sima, July 2012, Sep. 2012, Oct. 2012, Nov. 2012.
C
C     ******************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C ..  Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D+0 )
      COMPLEX*16        CZERO
      PARAMETER         ( CZERO = ( 0.0D+0, 0.0D+0 ) )
C
C     .. Mex-file interface parameters ..
      mwPointer         PLHS( * ), PRHS( * )
      INTEGER*4         NLHS, NRHS
C
C     .. Mex-file integer functions ..
      mwPointer         mxCalloc, mxCreateDoubleMatrix,
     $                  mxGetPr, mxGetPi
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C     .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         BALANC, JOB, JOBU
      INTEGER           ILO, INFO, LDA, LDQG, LDU1, LDU2, LDWORK,
     $                  LZWORK, N
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: BWORK( : )
      mwPointer         A, DW, DWORK, QG, SCALE, U1, U2, WI, WR, ZWORK
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           IA, IBAL, IJOB, IJOBU, IP, IQG, IU1, IU2,
     $                  LBWORK, LZWTMP, K, L, M1, MN, N1, N2, NN, NN2
      DOUBLE PRECISION  DUM(  2 ), TEMP
      COMPLEX*16        ZDUM( 1 )
C
C     .. External Subroutines ..
      EXTERNAL          MA02ES, MA02EZ, MA02KV, MA02KW, MA02LD, MA02LZ,
     $                  MB03XZ, MB04DZ, ZLACPY
C
C     ..Intrinsic Functions..
      INTRINSIC         INT, MAX, MIN
C
C     Check for proper number of arguments.
C 
      IF( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HAEIGZ requires at least 2 input arguments.' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HAEIGZ requires at least 1 output arguments.' )
      END IF
C
C   A(nxn), QG(nx(n+1)) (, job, jobu, balanc).
C
      N  = mxGetM( PRHS(1) )
      N1 = mxGetN( PRHS(1) )
      IF ( MIN( N, N1 ).GT.0 ) THEN
         IF ( mxIsNumeric( PRHS(1) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'A must be a real or complex matrix' )
         END IF
      END IF
      IF ( N.NE.N1 ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      N2  = 2*N
      NN  = N*N
      NN2 = N2*N2
C
      M1 = mxGetM( PRHS(2) )
      N1 = mxGetN( PRHS(2) )
      IF ( MIN( M1, N1 ).GT.1 ) THEN
         IF ( mxIsNumeric( PRHS(2) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'QG must be a real or complex matrix' )
         END IF
      END IF
      IF ( M1.NE.N ) THEN
         CALL mexErrMsgTxt( 'QG must have the same number of rows as A'
     $                    )
      END IF
      IF ( N1.NE.N+1 ) THEN
         CALL mexErrMsgTxt( 'QG must have one more columns than rows' )
      END IF
C
      IP = 3
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'JOB must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IJOB = TEMP
         IF ( IJOB.LT.-1 .OR. IJOB.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'JOB has -1, 0, 1, or 2 the only admissible values' )
         END IF
         IP = IP + 1
      ELSE
         IJOB = 0
      END IF
C
      IF ( IJOB.EQ.0 ) THEN
         JOB = 'E'
      ELSE IF ( IJOB.EQ.1 ) THEN
         JOB = 'S'
      ELSE IF ( IJOB.EQ.2 ) THEN
         JOB = 'G'
      END IF
C
      IF ( IJOB.GE.0 .AND. NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'JOBU must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'JOBU must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IJOBU = TEMP
         IF ( IJOBU.LT.0 .OR. IJOBU.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'JOBU has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
C
      ELSE
         IJOBU = 0
      END IF
C
      IF ( IJOB.GE.0 ) THEN
         IF ( IJOBU.EQ.0 ) THEN
            JOBU = 'N'
         ELSE
            JOBU = 'U'
         END IF
      ELSE IF( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HAEIGZ requires at least 4 input arguments.' )
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'BALANC must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'BALANC must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IBAL = TEMP
         IF ( IBAL.LT.0 .OR. IBAL.GT.3 ) THEN
            CALL mexErrMsgTxt
     $         ( 'BALANC has 0, 1, 2, or 3 the only admissible values' )
         END IF
C
      ELSE
         IBAL = 0
      END IF
C
      IF ( IBAL.EQ.0 ) THEN
         BALANC = 'N'
      ELSE IF ( IBAL.EQ.1 ) THEN
         BALANC = 'P'
      ELSE IF ( IBAL.EQ.2 ) THEN
         BALANC = 'S'
      ELSE
         BALANC = 'B'
      END IF
C
C Determine dimensions of the arrays and the workspace.
C
      IF( IJOB.LT.0 ) THEN
         K  = N
         IA = 1
      ELSE IF( IJOB.GT.0 ) THEN
         K  = N2
         IA = 1
         LBWORK = N2
      ELSE
         K  = N
         IA = 3
         LBWORK = 0
      END IF
C
      IF( IJOB.LT.0 ) THEN
         IQG = IA + NN
         L   = K  + 1
      ELSE
         IQG = IA + NN2
         L   = MIN( K+1, N2 )
         LZWTMP = MAX( 1, NN + N )
      END IF
C
      IF ( IJOB.EQ.0 .AND. IJOBU.GT.0 ) THEN
         IU1 = IQG + NN2 + N2
         IU2 = IU1 + NN2
      END IF
C
      LDA  = MAX( 1, K )
      LDQG = LDA
      IF ( IJOBU.EQ.0 ) THEN
         N1   = 0
         LDU1 = 1
         LDU2 = 1
      ELSE
         N1   = N2
         LDU1 = MAX( 1, N2 )
         LDU2 = MAX( 1, N2 )
      END IF
C
C Allocate part of the variable dimension local arrays.
C !Fortran 90/95
C
      A  = mxCalloc(  LDA*K, 16 )
      QG = mxCalloc( LDQG*L, 16 )
      IF( IJOB.GE.0 ) THEN
         ALLOCATE( BWORK( LBWORK ) )
         ZWORK = mxCalloc( LZWTMP, 16 )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      MN = MAX( 1, N )
      IF ( mxIsComplex( PRHS(1) ).EQ.0 ) THEN
         DW = mxCalloc( NN, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), %VAL( DW ),
     $                          NN )
         IF( IJOB.GE.0 ) THEN
            CALL MA02KW( NN, 1, 1, %VAL( DW ), %VAL( ZWORK ) )
            CALL ZLACPY( 'Full', N, N, %VAL( ZWORK ), MN, %VAL( A ),
     $                   LDA )
         ELSE
            CALL MA02KW( NN, 1, 1, %VAL( DW ), %VAL( A ) )
         END IF
         CALL mxFree( DW )
      ELSE IF( IJOB.GE.0 ) THEN
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS(1) ),
     $                              mxGetPi( PRHS(1) ),
     $                              %VAL( ZWORK ), NN )
         CALL ZLACPY( 'Full', N, N, %VAL( ZWORK ), MN, %VAL( A ), LDA )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS(1) ),
     $                              mxGetPi( PRHS(1) ),
     $                              %VAL( A ), NN )
      END IF
C
      IF ( mxIsComplex( PRHS(2) ).EQ.0 ) THEN
         DW = mxCalloc( NN+N, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( DW ),
     $                          NN+N )
         IF( IJOB.GE.0 ) THEN
            CALL MA02KW( NN+N, 1, 1, %VAL( DW ), %VAL( ZWORK ) )
            CALL ZLACPY( 'Full', N, N+1, %VAL( ZWORK ), MN, %VAL( QG ),
     $                   LDQG )
         ELSE
            CALL MA02KW( NN+N, 1, 1, %VAL( DW ), %VAL( QG ) )
         END IF
         CALL mxFree( DW )
      ELSE IF( IJOB.GE.0 ) THEN
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS(2) ),
     $                              mxGetPi( PRHS(2) ),
     $                              %VAL( ZWORK ), NN+N )
         CALL ZLACPY( 'Full', N, N+1, %VAL( ZWORK ), MN, %VAL( QG ),
     $                LDQG )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS(2) ),
     $                              mxGetPi( PRHS(2) ),
     $                              %VAL( QG ), NN+N )
      END IF
      IF ( IJOB.GE.0 )
     $   CALL mxFree( ZWORK )
C
      IF ( IJOB.GE.0 ) THEN
C
C Compute optimal ldwork, lzwork.
C
         CALL MB03XZ( BALANC, JOB, JOBU, N, ZDUM, LDA, ZDUM, LDQG, ZDUM,
     $                LDU1, ZDUM, LDU2, DUM, DUM, ILO, DUM, DUM, -1,
     $                ZDUM, -1, BWORK, INFO )
         LDWORK = INT(  DUM( 1 ) )
         LZWORK = INT( ZDUM( 1 ) )
      END IF
C
C Allocate the remaining variable dimension local arrays.
C !Fortran 90/95
C
      SCALE = mxCalloc( N, 8 )
      IF ( IJOB.GE.0 ) THEN
         WI    = mxCalloc(      N2,  8 )
         WR    = mxCalloc(      N2,  8 )
         DWORK = mxCalloc(  LDWORK,  8 )
         U1    = mxCalloc( LDU1*N1, 16 )
         U2    = mxCalloc( LDU2*N1, 16 )
         ZWORK = mxCalloc(  LZWORK, 16 )
      END IF
C
C Go the actual computations.
C
      IF ( IJOB.LT.0 ) THEN
         CALL MB04DZ( BALANC, N, %VAL( A ), LDA, %VAL( QG ), LDQG, ILO,
     $                %VAL( SCALE ), INFO )
      ELSE
         CALL MB03XZ( BALANC, JOB, JOBU, N, %VAL( A ), LDA, %VAL( QG ),
     $                LDQG, %VAL( U1 ), LDU1, %VAL( U2 ), LDU2,
     $                %VAL( WR ), %VAL( WI ), ILO, %VAL( SCALE ),
     $                %VAL( DWORK ), LDWORK, %VAL( ZWORK ), LZWORK,
     $                BWORK, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         IF ( IJOB.LT.0 ) THEN
            PLHS( 1 ) = mxCreateDoubleMatrix( 1, 1, 0 )
            TEMP = ILO
            CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS( 1 ) ), 1 )
C
            IF ( NLHS.GE.2 ) THEN
               PLHS( 2 ) = mxCreateDoubleMatrix( N, 1, 0 )
               CALL mxCopyReal8ToPtr( %VAL( SCALE ),
     $                                mxGetPr( PLHS( 2 ) ), N )
            END IF
            IP = 3
         ELSE
            PLHS( 1 ) = mxCreateDoubleMatrix( N2, 1, 1 )
            CALL mxCopyReal8ToPtr( %VAL( WR ),
     $                             mxGetPr( PLHS( 1 ) ), N2 )
            CALL mxCopyReal8ToPtr( %VAL( WI ),
     $                             mxGetPi( PLHS( 1 ) ), N2 )
            IP = 2
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            IF ( IJOB.LT.0 ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( N, N, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( A ),
     $                                    mxGetPr( PLHS( IP ) ),
     $                                    mxGetPi( PLHS( IP ) ),
     $                                    NN )
               IP = IP + 1
C
               IF ( NLHS.GE.IP ) THEN
                  PLHS( IP ) = mxCreateDoubleMatrix( N, N+1, 1 )
                  CALL mxCopyComplex16ToPtr( %VAL( QG ),
     $                                       mxGetPr( PLHS( IP ) ),
     $                                       mxGetPi( PLHS( IP ) ),
     $                                       NN+N )
                  IP = IP + 1
               END IF
            ELSE IF ( IJOB.GT.0 ) THEN
               CALL MA02LZ( 'Lower', N2-1, N2-1, 2, 1, CZERO, CZERO,
     $                      %VAL( A ), LDA )
               PLHS( IP ) = mxCreateDoubleMatrix( N2, N2, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( A ),
     $                                    mxGetPr( PLHS( IP ) ),
     $                                    mxGetPi( PLHS( IP ) ),
     $                                    NN2 )
               IP = IP + 1
C
               IF ( NLHS.GE.IP .AND. IJOB.NE.1 ) THEN
                  CALL MA02EZ( 'Upper', 'Conjugate', 'Not skew', N2, 
     $                         %VAL( QG ), LDQG )
                  PLHS( IP ) = mxCreateDoubleMatrix( N2, N2, 1 )
                  CALL mxCopyComplex16ToPtr( %VAL( QG ),
     $                                       mxGetPr( PLHS( IP ) ),
     $                                       mxGetPi( PLHS( IP ) ),
     $                                       NN2 )
                  IP = IP + 1
               END IF
            ELSE
               DW = mxCalloc( NN2, 8 )
               CALL MA02KV( NN2, IA, 1, %VAL( DWORK ), %VAL( DW ) )
               IF ( N2.GT.2 )
     $            CALL MA02LD( 'Lower', N2-2, N2-2, 3, 1, ZERO, ZERO,
     $                         %VAL( DW ), N2 )
               PLHS( IP ) = mxCreateDoubleMatrix( N2, N2, 0 )
               CALL mxCopyReal8ToPtr( %VAL( DW ),
     $                                mxGetPr( PLHS( IP ) ), NN2 )
               IP = IP + 1
C
               IF ( NLHS.GE.IP .AND. IJOB.NE.1 ) THEN
                  CALL MA02KV( NN2, IQG+N2, 1, %VAL( DWORK ),
     $                         %VAL( DW ) )
                  CALL MA02ES( 'Upper', N2, %VAL( DW ), N2 )
                  PLHS( IP ) = mxCreateDoubleMatrix( N2, N2, 0 )
                  CALL mxCopyReal8ToPtr( %VAL( DW ),
     $                                   mxGetPr( PLHS( IP ) ),
     $                                   NN2 )
                  IP = IP + 1
               END IF
               IF ( IJOBU.EQ.0 )
     $            CALL mxFree( DW )
            END IF
         END IF
C
         IF ( NLHS.GE.IP .AND. IJOBU.GT.0 ) THEN
            IF ( IJOB.GT.0 ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( N2, N2, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( U1 ),
     $                                    mxGetPr( PLHS( IP ) ),
     $                                    mxGetPi( PLHS( IP ) ),
     $                                    NN2 )
               IP = IP + 1
C
               IF ( NLHS.GE.IP ) THEN
                  PLHS( IP ) = mxCreateDoubleMatrix( N2, N2, 1 )
                  CALL mxCopyComplex16ToPtr( %VAL( U2 ),
     $                                       mxGetPr( PLHS( IP ) ),
     $                                       mxGetPi( PLHS( IP ) ),
     $                                       NN2 )
                  IP = IP + 1
               END IF
            ELSE
               PLHS( IP ) = mxCreateDoubleMatrix( N2, N2, 0 )
               CALL MA02KV( NN2, IU1, 1, %VAL( DWORK ), %VAL( DW ) )
               CALL mxCopyReal8ToPtr( %VAL( DW ),
     $                                mxGetPr( PLHS( IP ) ), NN2 )
               IP = IP + 1
C
               IF ( NLHS.GE.IP ) THEN
                  PLHS( IP ) = mxCreateDoubleMatrix( N2, N2, 0 )
                  CALL MA02KV( NN2, IU2, 1, %VAL( DWORK ),
     $                         %VAL( DW ) )
                  CALL mxCopyReal8ToPtr( %VAL( DW ),
     $                                   mxGetPr( PLHS( IP ) ), NN2 )
                  IP = IP + 1
               END IF
               CALL mxFree( DW )
            END IF
         END IF
C
         IF ( IJOB.GE.0 .AND. NLHS.GE.IP ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( 1, 1, 0 )
            TEMP = ILO
            CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS( IP ) ), 1 )
            IP = IP + 1
C
            IF ( NLHS.GE.IP ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( N, 1, 0 )
               CALL mxCopyReal8ToPtr( %VAL( SCALE ),
     $                                mxGetPr( PLHS( IP ) ), N )
            END IF
         END IF
C
      END IF
C
C Deallocate variable dimension arrays.
C !Fortran 90/95
C
      CALL mxFree( A     )
      CALL mxFree( QG    )
      CALL mxFree( SCALE )
      IF ( IJOB.GE.0 ) THEN
         DEALLOCATE(  BWORK )
         CALL mxFree( DWORK )
         CALL mxFree( U1    )
         CALL mxFree( U2    )
         CALL mxFree( WI    )
         CALL mxFree( WR    )
         CALL mxFree( ZWORK )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( IJOB.LT.0 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB04DZ'')'
     $           ) INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB03XZ'')'
     $           ) INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of HaeigZ ***
      END
