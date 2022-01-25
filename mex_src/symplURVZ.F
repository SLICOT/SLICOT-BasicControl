#include "fintrf.h"
C SYMPLURVZ.F - Gateway function for computing the eigenvalues of a
C               complex skew-Hamiltonian/Hamiltonian pencil in factored
C               form, using SLICOT routine MB04AZ.
C               The gateway accepts real or complex input matrices.
C
C RELEASE 2.0 of SLICOT Basic Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C                                   
C   [ALPHA,BETA(,Ao,Do,Co,Bo,Fo)(,Q,U)(,norms,unrel,ublk,publk)]
C                      = symplURVZ(Z,B,FG(,job,compq,compu,addinf))
C   [ALPHA,BETA(,Q,U)(,norms,unrel,ublk,publk)]
C                      = symplURVZ(Z,B,FG(,0,compq,compu,addinf))
C
C Purpose:
C   To compute the eigenvalues of a complex N-by-N skew-Hamiltonian/
C   Hamiltonian pencil aS - bH, with
C
C            H  T           (  B  F  )       (  Z11  Z12  )
C     S = J Z  J  Z and H = (      H ), Z =: (            ).         (1)
C                           (  G -B  )       (  Z21  Z22  )
C
C   The structured Schur form of the embedded real skew-Hamiltonian/
C                                                         H  T
C   skew-Hamiltonian pencil, aB_S - bB_T, with B_S = J B_Z  J  B_Z,
C
C           (  Re(Z11)  -Im(Z11)  |  Re(Z12)  -Im(Z12)  )
C           (                     |                     )
C           (  Im(Z11)   Re(Z11)  |  Im(Z12)   Re(Z12)  )
C           (                     |                     )
C     B_Z = (---------------------+---------------------) ,
C           (                     |                     )
C           (  Re(Z21)  -Im(Z21)  |  Re(Z22)  -Im(Z22)  )
C           (                     |                     )
C           (  Im(Z21)   Re(Z21)  |  Im(Z22)   Re(Z22)  )
C                                                                    (2)
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
C   is determined and used to compute the eigenvalues. Optionally, if
C   job = 1, the pencil aB_S - bB_H is transformed by a unitary matrix Q
C   and a unitary symplectic matrix U to the structured Schur form
C                                            H  T
C   aB_Sout - bB_Hout, with B_Sout = J B_Zout  J  B_Zout,
C
C              ( Ao  Do  )              ( Bo  Fo  )
C     B_Zout = (         ) and B_Hout = (       H ),                 (3)
C              (  0  Co  )              (  0 -Bo  )
C
C   where Ao and Bo are upper triangular, Co is lower triangular, and Fo
C   is Hermitian. B_H above is defined as B_H = -i*B_T. The embedding
C   doubles the multiplicities of the eigenvalues of the pencil aS - bH.
C   Optionally, if compq > 0, the unitary matrix Q is computed.
C   Optionally, if compu > 0, the unitary symplectic matrix U is
C   computed.
C
C Input parameters:
C   Z      - the n-by-n matrix Z, n = 2*m.
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
C            transformation matrix Q is returned, or if Q is not
C            required, as follows:
C            = 0 :  Q is not required (default);
C            = 1 :  on exit, Q contains the unitary matrix Q.
C   compu  - (optional) scalar indicating whether the unitary symplectic
C            transformation matrix U is returned, or if U is not
C            required, as follows:
C            = 0 :  U is not required (default);
C            = 1 :  on exit, U contains the relevant part of the unitary 
C                   symplectic matrix U (see below).
C   addinf - (optional) scalar indicating whether additional information
C            is required for numerically infinite or unreliable
C            eigenvalues, as follows:
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
C   Ao     - (optional) if job = 1, the n-by-n matrix Ao in (3).
C            The strictly lower triangular part is not zeroed.
C   Do     - (optional) if job = 1, the n-by-n matrix Do in (3).
C   Co     - (optional) if job = 1, the n-by-n matrix Co in (3).
C            The part over the first superdiagonal is not zeroed.
C   Bo     - (optional) if job = 1, the n-by-n matrix Bo in (3).
C            The strictly lower triangular part is not zeroed.
C   Fo     - (optional) if job = 1, the n-by-n matrix Fo in (3).
C            The strictly lower triangular part is not set.
C   Q      - (optional) if compq > 0, a 2*n-by-2*n matrix containing the
C            unitary matrix Q.
C   U      - (optional) if compu > 0, the n-by-2*n leading rows of the
C            unitary symplectic matrix U.
C            If job = 0 and compq > 0 or compu > 0, Q or U contain the
C            corresponding matrices which transformed B_Z and B_T
C            in (2). Their imaginary part is zero.
C   norms  - (optional) if addinf = 1, vector of length 3 containing the
C            Frobenius norms of the three factors of the formal matrix
C            product used by the algorithm. These norms can be used to
C            decide that some eigenvalues should considered as
C            numerically infinite or unreliable.
C   unrel  - (optional) if addinf = 1, the absolute values in unrel are
C            indices of the possibly inaccurate eigenvalues. Positive
C            values correspond to 1-by-1 diagonal blocks (with either
C            real or purely imaginary eigenvalues), and negative values
C            correspond to 2-by-2 diagonal blocks.
C   ublk   - (optional) if addinf = 1, vector containing the s triplet
C            values corresponding to the 1-by-1 blocks, with real or
C            purely imaginary eigenvalues, followed by t groups of
C            triplet 2-by-2 matrices corresponding to the 2-by-2
C            blocks. These eigenvalues are given by by
C            imag(ev) - real(ev)*i, where i is the purely imaginary
C            unit, and ev is the spectrum of the matrix product
C            A1*inv(A2)*inv(A3), with A1, ..., A3 defining the
C            corresponding 1-by-1 or 2-by-2 matrix triplet.
C   publk  - (optional) if addinf = 1, vector of pointers to the
C            starting location in ublk of each triplet of blocks
C            corresponding to the unreliable eigenvalues.
C            If addinf = 1, but there are no inaccurate eigenvalues,
C            then unrel, ublk, and publk are void.
C
C Contributors:
C   V. Sima, Research Institute for Informatics, Bucharest, Oct. 2009.
C
C Revisions:
C   V. Sima, Feb. 2010, Nov. 2010, Dec. 2010, July 2011, July 2012,
C   Sep. 2012, July 2013, May 2020.
C
C     ******************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
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
      CHARACTER         COMPQ, COMPU, JOB
      INTEGER           INFO, LDB, LDC, LDD, LDFG, LDQ, LDU, LDWORK,
     $                  LDZ, LIWORK, LZWORK, N
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: BWORK( : )
      INTEGER,          ALLOCATABLE :: IWORK( : )
      DOUBLE PRECISION, ALLOCATABLE :: IND( : )
      mwPointer         ALPHAI, ALPHAR, B, BETA, C, D, DW, DWORK, FG, Q,
     $                  TMP, U, Z, ZWORK
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           ADDINF, I, ICMPQ, ICMPU, IJOB, IP, J, L, LBWORK,
     $                  M, M1, M1M, N1, N2, NB, NC, NF, NN, P
      DOUBLE PRECISION  DDUM( 4 ), TEMP
      COMPLEX*16        DUM(  1 )
C
C     .. External Subroutines ..
      EXTERNAL          MA02KW, MB04AZ, ZLACPY
C
C     ..Intrinsic Functions..
      INTRINSIC         INT, MAX, MOD
C
C     Check for proper number of arguments.
C
      IF( NRHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYMPLURVZ requires at least 3 input arguments.' )
      ELSE IF ( NLHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYMPLURVZ requires at least 2 output arguments.' )
      END IF
C
C   Z(nxn), B(mxm), FG(mx(m+1))(, job, compq, compu(, addinf)).
C
      N  = mxGetM( PRHS(1) )
      N1 = mxGetN( PRHS(1) )
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 ) THEN
         CALL mexErrMsgTxt( 'Z must be a real or complex matrix' )
      END IF
      IF ( N.NE.N1 ) THEN
         CALL mexErrMsgTxt( 'Z must be a square matrix' )
      END IF
      IF ( MOD( N, 2 ).NE.0 ) THEN
         CALL mexErrMsgTxt( 'Z must have an even order' )
      END IF
      M = N/2
C
      M1 = mxGetM( PRHS(2) )
      N1 = mxGetN( PRHS(2) )
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 ) THEN
         CALL mexErrMsgTxt( 'B must be a real or complex matrix' )
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'B must have half the number of rows of Z' )
      END IF
      IF ( N1.NE.M ) THEN
         CALL mexErrMsgTxt( 'B must be a square matrix' )
      END IF
C
      M1 = mxGetM( PRHS(3) )
      N1 = mxGetN( PRHS(3) )
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 ) THEN
         CALL mexErrMsgTxt( 'FG must be a real or complex matrix' )
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'FG must have the same number of rows as B')
      END IF
      M1 = M + 1
      IF ( N1.NE.M1 ) THEN
         CALL mexErrMsgTxt( 'FG must have one more columns than rows' )
      END IF
C
      IP = 4
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
C
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
            CALL mexErrMsgTxt( 'COMPU must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'COMPU must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ICMPU = TEMP
         IF ( ICMPU.LT.0 .OR. ICMPU.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPU has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
C
      ELSE
         ICMPU = 0
      END IF
C
      IF ( ICMPU.EQ.0 ) THEN
         COMPU = 'N'
      ELSE
         COMPU = 'C'
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
      NN  = 2*N
      LDZ = MAX( 1, N )
      M1M = MAX( 1, M )
      IF ( IJOB.EQ.0 ) THEN
         LDB = M1M
         NB  = M
         NF  = M1
         LDC = 1
         NC  = 0
      ELSE
         LDB = LDZ
         NB  = N
         NF  = MAX( N, M1 )
         LDC = LDZ
         NC  = N
      END IF
      LDD  = LDC
      LDFG = LDB
      IF ( ICMPQ.EQ.0 ) THEN
         LDQ = 1
         N1  = 0
      ELSE
         LDQ = MAX( 1, NN )
         N1  = NN
      END IF
      IF ( ICMPU.EQ.0 ) THEN
         LDU = 1
         N2  = 0
      ELSE
         LDU = LDZ
         N2  = NN
      END IF
      LIWORK = NN + 9
      IF ( IJOB.EQ.0 ) THEN
         LBWORK = 0
      ELSE
         LBWORK = N
      END IF
C
      ALLOCATE( BWORK( LBWORK ), IND( M ), IWORK( LIWORK ) )
C
      CALL MB04AZ( JOB, COMPQ, COMPU, N, DUM, LDZ, DUM, LDB, DUM, LDFG,
     $             DUM, LDD, DUM, LDC, DUM, LDQ, DUM, LDU, DDUM, DDUM,
     $             DDUM, IWORK, LIWORK, DDUM, -1, DUM, -1, BWORK, INFO )
      LDWORK = INT( DDUM( 1 ) )
      LZWORK = INT(  DUM( 1 ) )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALPHAI = mxCalloc(       N,  8 )
      ALPHAR = mxCalloc(       N,  8 )
      B      = mxCalloc(  LDB*NB, 16 )
      BETA   = mxCalloc(       N,  8 )
      C      = mxCalloc(  LDC*NC, 16 )
      D      = mxCalloc(  LDD*NC, 16 )
      DWORK  = mxCalloc(  LDWORK,  8 )
      FG     = mxCalloc( LDFG*NF, 16 )
      Q      = mxCalloc(  LDQ*N1, 16 )
      U      = mxCalloc(  LDU*N2, 16 )
      Z      = mxCalloc(  LDZ*N,  16 )
      ZWORK  = mxCalloc(  LZWORK, 16 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      IF ( mxIsComplex( PRHS( 1 ) ).EQ.0 ) THEN
         DW = mxCalloc( N*N, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 1 ) ), %VAL( DW ),
     $                          N*N )
         CALL MA02KW( N*N, 1, 1, %VAL( DW ), %VAL( Z ) )
         CALL mxFree( DW )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 1 ) ),
     $                              mxGetPi( PRHS( 1 ) ),
     $                              %VAL( Z ), N*N )
      END IF
C
      TMP = mxCalloc( M*M1, 16 )
      IF ( mxIsComplex( PRHS( 2 ) ).EQ.0 ) THEN
         DW = mxCalloc( M*M, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 2 ) ), %VAL( DW ),
     $                          M*M )
         CALL MA02KW( M*M, 1, 1, %VAL( DW ), %VAL( TMP ) )
         CALL mxFree( DW )
         CALL ZLACPY( 'Full', M, M, %VAL( TMP ), M1M, %VAL( B ), LDB )
      ELSE IF( IJOB.EQ.0 ) THEN
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 2 ) ),
     $                              mxGetPi( PRHS( 2 ) ), %VAL( B ),
     $                              M*M )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 2 ) ),
     $                              mxGetPi( PRHS( 2 ) ),
     $                              %VAL( TMP ), M*M )
         CALL ZLACPY( 'Full', M, M, %VAL( TMP ), M1M, %VAL( B ), LDB )
      END IF
C
      IF ( mxIsComplex( PRHS( 3 ) ).EQ.0 ) THEN
         DW = mxCalloc( M*M1, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 3 ) ), %VAL( DW ),
     $                          M*M1 )
         CALL MA02KW( M*M1, 1, 1, %VAL( DW ), %VAL( TMP ) )
         CALL mxFree( DW )
         CALL ZLACPY( 'Full', M, M1, %VAL( TMP ), M1M, %VAL( FG ),
     $                LDFG )
      ELSE IF( IJOB.EQ.0 ) THEN
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 3 ) ),
     $                              mxGetPi( PRHS( 3 ) ),
     $                              %VAL( FG ), M*M1 )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 3 ) ),
     $                              mxGetPi( PRHS( 3 ) ),
     $                              %VAL( TMP ), M*M1 )
         CALL ZLACPY( 'Full', M, M1, %VAL( TMP ), M1M, %VAL( FG ),
     $                LDFG )
      END IF
      CALL mxFree( TMP )
C
C Do the actual computations.
C
      CALL MB04AZ( JOB, COMPQ, COMPU, N, %VAL( Z ), LDZ, %VAL( B ),
     $             LDB, %VAL( FG ), LDFG, %VAL( D ), LDD, %VAL( C ),
     $             LDC, %VAL( Q ), LDQ, %VAL( U ), LDU, %VAL( ALPHAR ),
     $             %VAL( ALPHAI ), %VAL( BETA ), IWORK, LIWORK,
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
         CALL mxCopyReal8ToPtr( %VAL( BETA ), mxGetPr( PLHS( 2 ) ),
     $                          N )
C
         IP = 3
         IF ( IJOB.GT.0 ) THEN
            IF ( NLHS.GE.IP ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( N, N, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( Z ),
     $                                    mxGetPr( PLHS( IP ) ),
     $                                    mxGetPi( PLHS( IP ) ),
     $                                    N*N )
               IP = IP + 1
            END IF
C
            IF ( NLHS.GE.IP ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( N, N, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( D ),
     $                                    mxGetPr( PLHS( IP ) ),
     $                                    mxGetPi( PLHS( IP ) ),
     $                                    N*N )
               IP = IP + 1
            END IF
C
            IF ( NLHS.GE.IP ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( N, N, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( C ),
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
         END IF
C
         IF ( NLHS.GE.IP .AND. ICMPQ.GT.0 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( NN, NN, 1 )
            CALL mxCopyComplex16ToPtr( %VAL( Q ),
     $                                 mxGetPr( PLHS( IP ) ),
     $                                 mxGetPi( PLHS( IP ) ),
     $                                 NN*NN )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. ICMPU.GT.0 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N, NN, 1 )
            CALL mxCopyComplex16ToPtr( %VAL( U ),
     $                                 mxGetPr( PLHS( IP ) ),
     $                                 mxGetPi( PLHS( IP ) ),
     $                                 N*NN )
         END IF
C
         IF ( NLHS.GE.IP .AND. ADDINF.GT.0 ) THEN
            CALL MA02KV( 3, 2, 1, %VAL( DWORK ), DDUM )
            PLHS( IP ) = mxCreateDoubleMatrix( 1, 3, 0 )
            CALL mxCopyReal8ToPtr( DDUM, mxGetPr( PLHS( IP ) ), 3 )
            IP = IP + 1
C
            IF ( NLHS.GE.IP ) THEN
               L = IWORK( 1 )
               DO 10 I = 1, L
                  IND( I ) = IWORK( I+1 )
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
               P = 3*( N - 2*J ) + 12*J
               DW = mxCalloc( P, 8 )
               CALL MA02KV( P, 5, 1, %VAL( DWORK ), %VAL( DW ) )
               PLHS( IP ) = mxCreateDoubleMatrix( MIN( 1, P ),
     $                                            P, 0 )
               CALL mxCopyReal8ToPtr( %VAL( DW ),
     $                                mxGetPr( PLHS( IP ) ), P )
               CALL mxFree( DW )
               IP = IP + 1
            END IF
C
            IF ( NLHS.GE.IP ) THEN
               DO 20 I = 1, L
                  IND( I ) = IWORK( L+1+I )
   20          CONTINUE
               PLHS( IP ) = mxCreateDoubleMatrix( MIN( 1, L ), L,
     $                                            0 )
               CALL mxCopyReal8ToPtr( IND, mxGetPr( PLHS( IP ) ),
     $                                L )
            END IF
         END IF
      END IF
C
C Deallocate variable dimension arrays.
C !Fortran 90/95
C
      DEALLOCATE( BWORK, IND, IWORK  )
      CALL mxFree( ALPHAI )
      CALL mxFree( ALPHAR )
      CALL mxFree( B      )
      CALL mxFree( BETA   )
      CALL mxFree( C      )
      CALL mxFree( D      )
      CALL mxFree( DWORK  )
      CALL mxFree( FG     )
      CALL mxFree( Q      )
      CALL mxFree( U      )
      CALL mxFree( Z      )
      CALL mxFree( ZWORK  )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 .AND. INFO.NE.3 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB04AZ'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of symplURV ***
      END
 
