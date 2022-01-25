#include "fintrf.h"
C SKEWHAMILEIGZ.F - Gateway function for computing the eigenvalues of a
C                   complex skew-Hamiltonian/Hamiltonian pencil, using
C                   SLICOT routine MB04BZ.
C                   The gateway accepts real or complex input matrices.
C
C RELEASE 2.0 of SLICOT Basic Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [ALPHA,BETA(,Ao,Do,Bo,Fo(,Q)(,norms,unrel,ublk,publk))] =
C                            skewHamileigZ(A,DE,B,FG(,job,compq,addinf))
C
C   [ALPHA,BETA]                 = skewHamileigZ(A,DE,B,FG)
C   [ALPHA,BETA,Qo]              = skewHamileigZ(A,DE,B,FG,0,1)
C   [ALPHA,BETA,Ao,Do,Bo,Fo]     = skewHamileigZ(A,DE,B,FG,1)
C   [ALPHA,BETA,...,norms,unrel] = skewHamileigZ(A,DE,B,FG,...,addinf)
C   [ALPHA,BETA,...,norms,unrel,ublk,publk]
C                                = skewHamileigZ(A,DE,B,FG,...,addinf)
C
C Purpose:
C   To compute the eigenvalues of a complex n-by-n skew-Hamiltonian/
C   Hamiltonian pencil aS - bH, with
C
C         (  A  D  )         (  B  F  )
C     S = (      H ) and H = (      H ).                           (1)
C         (  E  A  )         (  G -B  )
C
C   The structured Schur form of the embedded real skew-Hamiltonian/
C   skew-Hamiltonian pencil aB_S - bB_T, defined as
C
C           (  Re(A)  -Im(A)  |  Re(D)  -Im(D)  )
C           (                 |                 )
C           (  Im(A)   Re(A)  |  Im(D)   Re(D)  )
C           (                 |                 )
C     B_S = (-----------------+-----------------) , and
C           (                 |      T       T  )
C           (  Re(E)  -Im(E)  |  Re(A )  Im(A ) )
C           (                 |      T       T  )
C           (  Im(E)   Re(E)  | -Im(A )  Re(A ) )
C                                                                  (2)
C           ( -Im(B)  -Re(B)  | -Im(F)  -Re(F)  )
C           (                 |                 )
C           (  Re(B)  -Im(B)  |  Re(F)  -Im(F)  )
C           (                 |                 )
C     B_T = (-----------------+-----------------) ,  T = i*H,
C           (                 |      T       T  )
C           ( -Im(G)  -Re(G)  | -Im(B )  Re(B ) )
C           (                 |      T       T  )
C           (  Re(G)  -Im(G)  | -Re(B ) -Im(B ) )
C
C   is determined and used to compute the eigenvalues. Optionally,
C   if job = 1, the pencil aB_S - bB_H (B_H = -i*B_T) is transformed
C   by a unitary matrix Q to the structured Schur form
C
C              ( Ao  Do  )              ( Bo  Fo  )
C     B_Sout = (       H ) and B_Hout = (       H ),               (3)
C              (  0  Ao  )              (  0 -Bo  )
C
C   where Ao and Bo are upper triangular, Do is skew-Hermitian, and
C   Fo is Hermitian. The embedding doubles the multiplicities of the
C   eigenvalues of the pencil aS - bH. Optionally, if compq = 1, the
C   unitary matrix Q is computed.
C
C Input parameters:
C   A      - the m-by-m matrix A, with m = n/2.
C   DE     - an  m-by-(m+1) matrix containing the triangles of the
C            skew-Hermitian matrices D and E, as follows:
C            the leading m-by-m lower triangular part contains the
C            lower triangle of the matrix E, and the m-by-m upper
C            triangular part of the submatrix in the columns 2 to m+1
C            contains the upper triangle of the matrix D of S in (1).
C            So, if i >= j, then E(i,j) = -conj(E(j,i)) is stored in
C            DE(i,j) and D(j,i) = -conj(D(i,j)) is stored in DE(j,i+1).
C            DE is an empty matrix if m = 0.
C   B      - the m-by-m matrix B.
C   FG     - an  m-by-(m+1) matrix containing the triangles of the
C            Hermitian matrices F and G, as follows:
C            the leading m-by-m lower triangular part contains the lower
C            triangle of the matrix G, and the m-by-m upper triangular
C            part of the submatrix in the columns 2 to m+1 contains the
C            upper triangle of the matrix F of H in (1).
C            So, if i >= j, then G(i,j) = conj(G(j,i)) is stored in
C            FG(i,j) and F(j,i) = conj(F(i,j)) is stored in FG(j,i+1).
C            FG is an empty matrix if m = 0.
C   job    - (optional) scalar indicating the computation to be
C            performed, as follows:
C            = 0 :  compute the eigenvalues only (default);
C            = 1 :  compute the eigenvalues and the matrices of the
C                   transformed pencil in (3).
C   compq  - (optional) scalar indicating whether the unitary
C            transformation matrix Q is required, as follows:
C            = 0 :  Q is not required (default);
C            = 1 :  on exit, Q contains the unitary transformation
C                   matrix.
C   addinf - (optional) scalar indicating whether additional information
C            is required for possibly inaccurate eigenvalues as follows:
C            = 0 :  additional information is not required (default);
C            = 1 :  additional information is required.
C
C Output parameters:
C   ALPHA  - the n-vector of the numerators alpha defining the
C            eigenvalues of the pencil aS - bH.
C   BETA     the n-vector of the denominators beta defining the
C            eigenvalues of the pencil aS - bH.
C            Together, the quantities alpha = ALPHA(j) and
C            beta = BETA(j) represent the j-th eigenvalue of the pencil
C            aS - bH, in the form lambda = alpha/beta. Since lambda may
C            overflow, the ratios should not, in general, be computed.
C   Ao     - if job = 1, the computed n-by-n submatrix Ao in (3).
C            The strictly lower triangular part is not zeroed.
C   Do     - if job = 1, the computed n-by-n upper triangular part of
C            the skew-Hermitian matrix Do in (3). The strictly lower
C            triangular part is not set.
C   Bo     - if job = 1, the computed n-by-n submatrix Bo in (3).
C            The strictly lower triangular part is not zeroed.
C   Fo     - if job = 1, the computed n-by-n upper triangular part of
C            the Hermitian matrix Fo in (3). The strictly lower
C            triangular part is not set.
C   Q      - if compq = 1, an 2*n-by-2*n matrix containing the computed
C            unitary matrix.
C   norms  - (optional) if addinf = 1, vector of length 2 containing the
C            Frobenius norms of the matrices B_S and B_T in (2). These
C            norms are used in the tests to decide that some eigenvalues
C            are considered as unreliable (possibly inaccurate).
C   unrel  - (optional) if addinf = 1, vector containing the indices of
C            unreliable eigenvalues.
C   ublk   - (optional) if addinf = 1, vector containing the s double
C            values corresponding to the 1-by-1 blocks, with real or
C            purely imaginary eigenvalues, followed by t groups of
C            two 2-by-2 matrices corresponding to the 2-by-2 blocks.
C            The first one or four values in each group, respectively,
C            correspond to the T matrix and the following one or four
C            values correspond to the S matrix. The eigenvalues of a
C            1-by-1 group (T11, S11) are given by -1i*T11/S11, while
C            the eigenvalues of a 2-by-2 group are similarly given by
C            -1i*eig(T22,S22). If a 2-by-2 group starts at location i,
C            then T22 = [ublk(i) ublk(i+2); ublk(i+1) ublk(i+3)],
C            followed similarly by S22. Note that all 1-by-1 groups
C            appear before any of 2-by-2 groups, in the order of the
C            eigenvalues.
C   publk  - (optional) if addinf = 1, vector of pointers to the
C            starting location of each 1-by-1 or 2-by-2 pairs in ublk
C            corresponding to unrealible eigenvalues.
C
C Contributors:
C   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2011.
C
C Revisions:
C   V. Sima, July 2012, Sep. 2012, Oct. 2012, Sep. 2016, Jan. 2017,
C   Apr. 2020.
C
C     ******************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C     .. Parameters ..
      COMPLEX*16         CIMAG
      PARAMETER          ( CIMAG = ( 0.0D+0, 1.0D+0 ) )
C
C     .. Mex-file interface parameters ..
      mwPointer         PLHS( * ), PRHS( * )
      INTEGER*4         NLHS, NRHS
C
C     .. Mex-file integer functions ..
      mwPointer         mxCalloc, mxCreateDoubleMatrix,
     $                  mxGetPi, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C     .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         COMPQ, JOB
      INTEGER           INFO, LDA, LDB, LDDE, LDFG, LDQ, LDWORK,
     $                  LZWORK, N 
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: BWORK( : )
      INTEGER,          ALLOCATABLE :: IWORK( : )
      DOUBLE PRECISION, ALLOCATABLE :: IND( : )
      mwPointer         A, ALPHAI, ALPHAR, B, BETA, DE, DW, DWA, DWB,
     $                  DWORK, FG, Q, TMP, ZWORK
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           REALA, REALB, REALD, REALF
      INTEGER           ADDINF, I, ICMPQ, IJOB, IP, IS, J, K, L, LBWORK,
     $                  M, M1, MP1, N1, P
      DOUBLE PRECISION  DUM(  3 ), TEMP
      COMPLEX*16        ZDUM( 1 )
C
C     .. External Subroutines ..
      EXTERNAL          MA02KI, MA02KV, MA02KW, MB04BZ, ZLACPY
C
C     ..Intrinsic Functions..
      INTRINSIC         INT, MAX, MIN
C
C     Check for proper number of arguments.
C
      IF( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMILEIGZ requires at least 4 input arguments.' )
      ELSE IF ( NLHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMILEIGZ requires at least 2 output arguments.' )
      END IF
C
C   A(mxm), DE(mx(m+1)), B(mxm), FG(mx(m+1))(, job, compq, addinf).
C
      M  = mxGetM( PRHS(1) )
      N1 = mxGetN( PRHS(1) )
      IF ( MIN( M, N1 ).GT.0 ) THEN
         IF ( mxIsNumeric( PRHS(1) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'A must be a real or complex matrix' )
         END IF
      END IF
      IF ( M.NE.N1 ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      N   = 2*M
      MP1 = M + 1
C
      M1 = mxGetM( PRHS(2) )
      N1 = mxGetN( PRHS(2) )
      IF ( MIN( M1, N1 ).GT.1 ) THEN
         IF ( mxIsNumeric( PRHS(2) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'DE must be a real or complex matrix' )
         END IF
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'DE must have the same number of rows as A'
     $                    )
      END IF
      IF ( N1.NE.MP1 ) THEN
         CALL mexErrMsgTxt( 'DE must have one more columns than rows' )
      END IF
C
      M1 = mxGetM( PRHS(3) )
      N1 = mxGetN( PRHS(3) )
      IF ( MIN( M1, N1 ).GT.0 ) THEN
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'B must be a real or complex matrix' )
         END IF
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'B must have the same number of rows as A' )
      END IF
      IF ( N1.NE.M ) THEN
         CALL mexErrMsgTxt( 'B must be a square matrix' )
      END IF
C
      M1 = mxGetM( PRHS(4) )
      N1 = mxGetN( PRHS(4) )
      IF ( MIN( M1, N1 ).GT.1 ) THEN
         IF ( mxIsNumeric( PRHS(4) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'FG must be a real or complex matrix' )
         END IF
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'FG must have the same number of rows as A'
     $                    )
      END IF
      IF ( N1.NE.MP1 ) THEN
         CALL mexErrMsgTxt( 'FG must have one more columns than rows' )
      END IF
C
      IP = 5
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
         IF ( IJOB.LT.0 .OR. IJOB.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'JOB has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
      ELSE
         IJOB = 0
      END IF
C
      IF ( IJOB.EQ.0 ) THEN
         JOB = 'E'
      ELSE
         JOB = 'T'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'COMPQ must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'COMPQ must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ICMPQ = TEMP
         IF ( ICMPQ.LT.0 .OR. ICMPQ.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPQ has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
      ELSE
         ICMPQ = 0
      END IF
C
      IF ( ICMPQ.EQ.0 ) THEN
         COMPQ = 'N'
      ELSE
         COMPQ = 'C'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ADDINF must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ADDINF must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ADDINF = TEMP
         IF ( ADDINF.LT.0 .OR. ADDINF.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'ADDINF has 0 or 1 the only admissible values' )
         END IF
      ELSE
         ADDINF = 0
      END IF
C
C Determine dimensions of the arrays and the workspace.
C
      IF ( IJOB.EQ.0 ) THEN
         K      = M
         LBWORK = 0
      ELSE
         K      = N
         LBWORK = N
      END IF
      L    = MIN( K+1, N )
      LDA  = MAX( 1, K )
      LDDE = LDA
      LDB  = LDA
      LDFG = LDA
      IF ( ICMPQ.EQ.0 ) THEN
         N1  = 0
         LDQ = 1
      ELSE
         N1  = 2*N
         LDQ = MAX( 1, N1 )
      END IF
C
C Allocate part of the variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE( BWORK( LBWORK ), IND( N ), IWORK( 2*N+3 ) )
      A  = mxCalloc(  LDA*K, 16 )
      B  = mxCalloc(  LDB*K, 16 )
      DE = mxCalloc( LDDE*L, 16 )
      FG = mxCalloc( LDFG*L, 16 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      REALA = mxIsComplex( PRHS( 1 ) ).EQ.0
      REALD = mxIsComplex( PRHS( 2 ) ).EQ.0
      REALB = mxIsComplex( PRHS( 3 ) ).EQ.0
      REALF = mxIsComplex( PRHS( 4 ) ).EQ.0
      IF ( IJOB.NE.0 .OR. REALA .OR. REALD .OR. REALB .OR. REALF )
     $   TMP = mxCalloc( M*MP1, 16 )
      M1 = MAX( 1, M )
C
      IF ( REALA ) THEN
         DW = mxCalloc( M*M, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 1 ) ), %VAL( DW ),
     $                          M*M )
         CALL MA02KW( M*M, 1, 1, %VAL( DW ), %VAL( TMP ) )
         CALL mxFree( DW )
         CALL ZLACPY( 'Full', M, M, %VAL( TMP ), M1, %VAL( A ), LDA )
      ELSE IF ( IJOB.EQ.0 ) THEN
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 1 ) ),
     $                              mxGetPi( PRHS( 1 ) ),
     $                              %VAL( A ), M*M )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 1 ) ),
     $                              mxGetPi( PRHS( 1 ) ),
     $                              %VAL( TMP ), M*M )
         CALL ZLACPY( 'Full', M, M, %VAL( TMP ), M1, %VAL( A ), LDA )
      END IF
C
      IF ( REALD ) THEN
         DW = mxCalloc( M*MP1, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 2 ) ), %VAL( DW ),
     $                          M*MP1 )
         CALL MA02KW( M*MP1, 1, 1, %VAL( DW ), %VAL( TMP ) )
         CALL mxFree( DW )
         CALL ZLACPY( 'Full', M, MP1, %VAL( TMP ), M1, %VAL( DE ), LDDE
     $              )
      ELSE IF ( IJOB.EQ.0 ) THEN
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 2 ) ),
     $                              mxGetPi( PRHS( 2 ) ),
     $                              %VAL( DE ), M*MP1 )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 2 ) ),
     $                              mxGetPi( PRHS( 2 ) ),
     $                              %VAL( TMP ), M*MP1 )
         CALL ZLACPY( 'Full', M, MP1, %VAL( TMP ), M1, %VAL( DE ), LDDE
     $              )
      END IF
C
      IF ( REALB ) THEN
         DW = mxCalloc( M*M, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 3 ) ), %VAL( DW ),
     $                          M*M )
         CALL MA02KW( M*M, 1, 1, %VAL( DW ), %VAL( TMP ) )
         CALL mxFree( DW )
         CALL ZLACPY( 'Full', M, M, %VAL( TMP ), M1, %VAL( B ), LDB )
      ELSE IF ( IJOB.EQ.0 ) THEN
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 3 ) ),
     $                              mxGetPi( PRHS( 3 ) ),
     $                              %VAL( B ), M*M )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 3 ) ),
     $                              mxGetPi( PRHS( 3 ) ),
     $                              %VAL( TMP ),  M*M )
         CALL ZLACPY( 'Full', M, M, %VAL( TMP ), M1, %VAL( B ), LDB )
      END IF
C
      IF ( REALF ) THEN
         DW = mxCalloc( M*MP1, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 4 ) ), %VAL( DW ),
     $                          M*MP1 )
         CALL MA02KW( M*MP1, 1, 1, %VAL( DW ), %VAL( TMP ) )
         CALL mxFree( DW )
         CALL ZLACPY( 'Full', M, MP1, %VAL( TMP ), M1, %VAL( FG ), LDFG
     $              )
      ELSE IF ( IJOB.EQ.0 ) THEN
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 4 ) ),
     $                              mxGetPi( PRHS( 4 ) ),
     $                              %VAL( FG ), M*MP1 )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 4 ) ),
     $                              mxGetPi( PRHS( 4 ) ),
     $                              %VAL( TMP ), M*MP1 )
         CALL ZLACPY( 'Full', M, MP1, %VAL( TMP ), M1, %VAL( FG ), LDFG
     $              )
      END IF
      IF ( IJOB.NE.0 .OR. REALA .OR. REALD .OR. REALB .OR. REALF )
     $   CALL mxFree( TMP )
C
C   ldwork, lzwork.
C
      CALL MB04BZ( JOB, COMPQ, N, ZDUM, LDA, ZDUM, LDDE, ZDUM, LDB,
     $             ZDUM, LDFG, ZDUM, LDQ, DUM, DUM, DUM, IWORK, DUM, -1,
     $             ZDUM, -1, BWORK, INFO )
      LDWORK = INT(  DUM( 1 ) )
      LZWORK = INT( ZDUM( 1 ) )
C
C Allocate the remaining variable dimension local arrays.
C !Fortran 90/95
C
      ALPHAI = mxCalloc(      N,  8 )
      ALPHAR = mxCalloc(      N,  8 )
      BETA   = mxCalloc(      N,  8 )
      DWORK  = mxCalloc( LDWORK,  8 )
      Q      = mxCalloc( LDQ*N1, 16 )
      ZWORK  = mxCalloc( LZWORK, 16 )
C
C Do the actual computations.
C
      CALL MB04BZ( JOB, COMPQ, N, %VAL( A ), LDA, %VAL( DE ), LDDE,
     $             %VAL( B ), LDB, %VAL( FG ), LDFG, %VAL( Q ), LDQ,
     $             %VAL( ALPHAR ), %VAL( ALPHAI ), %VAL( BETA ), IWORK,
     $             %VAL( DWORK ), LDWORK, %VAL( ZWORK ), LZWORK, BWORK,
     $             INFO )
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 .OR. INFO.EQ.3 ) THEN
         PLHS( 1 ) = mxCreateDoubleMatrix( N, 1, 1 )
         CALL mxCopyReal8ToPtr( %VAL( ALPHAR ),
     $                          mxGetPr( PLHS( 1 ) ), N )
         CALL mxCopyReal8ToPtr( %VAL( ALPHAI ),
     $                          mxGetPi( PLHS( 1 ) ), N )
C
         PLHS( 2 ) = mxCreateDoubleMatrix( N, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( BETA ),
     $                          mxGetPr( PLHS( 2 ) ), N )
C
         IP = 3
         L  = IWORK( 1 )
         IF ( IJOB.EQ.1 ) THEN
            IF ( NLHS.GE.IP ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( N, N, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( A ),
     $                                    mxGetPr( PLHS( IP ) ),
     $                                    mxGetPi( PLHS( IP ) ),
     $                                    N*N )
               IP = IP + 1
            END IF
C
            IF ( NLHS.GE.IP ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( N, N, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( DE ),
     $                                    mxGetPr( PLHS( IP ) ),
     $                                    mxGetPi( PLHS( IP ) ),
     $                                    N*N )
               IP = IP + 1
            END IF
C
            IF ( NLHS.GE.IP ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( N, N, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( B ),
     $                                    mxGetPr( PLHS( IP ) ),
     $                                    mxGetPi( PLHS( IP ) ),
     $                                    N*N )
               IP = IP + 1
            END IF
C
            IF ( NLHS.GE.IP ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( N, N, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( FG ),
     $                                    mxGetPr( PLHS( IP ) ),
     $                                    mxGetPi( PLHS( IP ) ),
     $                                    N*N )
               IP = IP + 1
            END IF
C
         END IF
C
         IF ( NLHS.GE.IP .AND. ICMPQ.GT.0 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N1, N1, 1 )
            CALL mxCopyComplex16ToPtr( %VAL( Q ),
     $                                 mxGetPr( PLHS( IP ) ),
     $                                 mxGetPi( PLHS( IP ) ),
     $                                 N1*N1 )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. ADDINF.GT.0 ) THEN
            CALL MA02KV( 2, 2, 1, %VAL( DWORK ), DUM )
            PLHS( IP ) = mxCreateDoubleMatrix( 1, 2, 0 )
            CALL mxCopyReal8ToPtr( DUM, mxGetPr( PLHS( IP ) ), 2 )
            IP = IP + 1
C
            IF ( NLHS.GE.IP ) THEN
               DO 10 K = 1, L
                  IND( K ) = IWORK( K+1 )
   10          CONTINUE
               PLHS( IP ) = mxCreateDoubleMatrix( MIN( 1, L ), L,
     $                                            0 )
               CALL mxCopyReal8ToPtr( IND, mxGetPr( PLHS( IP ) ),
     $                                L )
               IP = IP + 1
            END IF
C
            IF ( NLHS.GE.IP ) THEN
               J = IWORK( 2*L+4 )
               J = 2*( N - 2*J ) + 8*J
               DWA = mxCalloc( J, 8 )
               CALL MA02KV( J, 4, 1, %VAL( DWORK ), %VAL( DWA ) )
               PLHS( IP ) = mxCreateDoubleMatrix( MIN( 1, J ), J,
     $                                            0 )
               CALL mxCopyReal8ToPtr( %VAL( DWA ),
     $                                mxGetPr( PLHS( IP ) ), J )
               CALL mxFree( DWA  )
               IP = IP + 1
            END IF
C
            IF ( NLHS.GE.IP ) THEN
               DO 20 K = 1, L
                  IND( K ) = IWORK( K+L+1 ) - 3
   20          CONTINUE
               PLHS( IP ) = mxCreateDoubleMatrix( MIN( 1, L ), L,
     $                                            0 )
               CALL mxCopyReal8ToPtr( IND, mxGetPr( PLHS( IP ) ),
     $                                L )
            END IF
C
         END IF
C
      END IF
C
C Deallocate variable dimension arrays.
C !Fortran 90/95
C
      DEALLOCATE( BWORK, IND, IWORK )
      CALL mxFree( A      )
      CALL mxFree( ALPHAI )
      CALL mxFree( ALPHAR )
      CALL mxFree( B      )
      CALL mxFree( BETA   )
      CALL mxFree( DE     )
      CALL mxFree( DWORK  )
      CALL mxFree( FG     )
      CALL mxFree( Q      )
      CALL mxFree( ZWORK  )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB04BZ'')'
     $        ) INFO
         IF ( INFO.EQ.3 ) THEN
            CALL mexPrintf( TEXT )
         ELSE
            CALL mexErrMsgTxt( TEXT )
         END IF
      END IF
C
      RETURN
C *** Last line of skewHamileigZ ***
      END
