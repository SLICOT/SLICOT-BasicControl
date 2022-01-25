#include "fintrf.h"
C SYMPLURV.F - Gateway function for computing the eigenvalues of a real
C              skew-Hamiltonian/Hamiltonian pencil in factored form,
C              using SLICOT routine MB04AD.
C
C RELEASE 2.0 of SLICOT Basic Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [ALPHAR,ALPHAI,BETA(,To,Zo,Ho(,Q1,Q2,U11,U12,U21,U22)
C                                   (,norms,unrel,ublk,publk))] =
C           symplURV(Z,H(,job,compq1(,Q01),compq2(,Q02),
C                       compu1(,U011,U012),compu2(,U021,U022)(,addinf)))
C   [ALPHAR,ALPHAI,BETA]          = symplURV(Z,H)
C   [ALPHAR,ALPHAI,BETA,To,Zo,Ho) = symplURV(Z,H,1)
C   [ALPHAR,ALPHAI,BETA,To,Zo,Ho,Q1,Q2,U11,U12,U21,U22)
C                                 = symplURV(Z,H,1,1,1,1,1)
C
C Purpose:
C   To compute the eigenvalues of a real skew-Hamiltonian/Hamiltonian
C   pencil aS - bH with
C
C                                    (  0  I  )
C     S = T Z = J Z' J' Z, where J = (        ),                   (1)
C                                    ( -I  0  )
C
C   via generalized symplectic URV decomposition. That is, orthogonal
C   matrices Q1 and Q2 and orthogonal symplectic matrices U1 and U2
C   are computed such that
C
C                                 (  T11  T12 )
C     Q1' T U1 = Q1' J Z' J' U1 = (           ) = To,
C                                 (   0   T22 )
C
C                (  Z11  Z12 )
C     U2' Z Q2 = (           ) = Zo,                             (2)
C                (   0   Z22 )
C
C                ( H11  H12 )
C     Q1' H Q2 = (          ) = Ho,
C                (  0   H22 )
C
C   where T11, T22', Z11, Z22', H11 are upper triangular and H22' is
C   upper quasi-triangular.
C   Optionally, if compq1 = 1 (compq2 = 1), the orthogonal matrix Q1
C   (Q2) that fulfills (2) is computed.
C   Optionally, if compu1 = 1 (compu2 = 1), the orthogonal symplectic
C   matrix
C
C          (  U11  U12  )          (  U21  U22  )
C     U1 = (            )   ( U2 = (            ) )
C          ( -U12  U11  )          ( -U22  U21  )
C
C   that fulfills (2) is computed. Only U11 and U12 (U21 and U22 ) are
C   returned.
C
C Input parameters:
C   Z      - the n-by-n matrix Z, n even.
C   H      - the n-by-n Hamiltonian matrix H.
C   job    - (optional) scalar indicating the computation to be
C            performed, as follows:
C            = 0 :  compute the eigenvalues only (default);
C            = 1 :  compute the eigenvalues and the matrices of the
C                   transformed pencil in (2).
C   compq1 - (optional) scalar indicating whether the orthogonal
C            transformation matrix Q1 is returned, or if Q1 is not
C            required, as follows:
C            = 0 :  Q1 is not required (default);
C            = 1 :  on exit, Q1 contains the orthogonal matrix Q1;
C            = 2 :  the orthogonal transformations are accumulated
C                   into Q1;
C                   on input, Q01 must contain an orthogonal matrix Q01;
C                   on exit, Q1 contains Q01*Q1, with Q1 returned for
C                   compq1 = 1.
C   Q01    - if compq1 = 2, the n-by-n orthogonal matrix Q01.
C   compq2 - (optional) scalar indicating whether the orthogonal
C            transformation matrix Q2 is returned, or if Q2 is not
C            required, as follows:
C            = 0 :  Q2 is not required (default);
C            = 1 :  on exit, Q2 contains the orthogonal matrix Q2;
C            = 2 :  the orthogonal transformations are accumulated
C                   into Q2;
C                   on input, Q02 must contain an orthogonal matrix Q02;
C                   on exit, Q2 contains Q02*Q2, with Q2 returned for
C                   compq2 = 1.
C   Q02    - if compq2 = 2, the n-by-n orthogonal matrix Q02.
C   compu1 - (optional) scalar indicating whether the orthogonal
C            symplectic transformation matrix U1 is returned, or if U1
C            is not required, as follows:
C            = 0 :  U1 is not required (default);
C            = 1 :  on exit, U11 and U12 contain the corresponding
C                   submatrices of the orthogonal symplectic matrix U1;
C            = 2 :  the orthogonal transformations are accumulated
C                   into U1;
C                   on input, U011, U012 must contain the corresponding
C                   submatrices of an orthogonal symplectic matrix U01;
C                   on exit, U11 and U12 contain the updated submatrices
C                   U11 and U12 of the matrix product U01*U1, with U1
C                   returned for compu1 = 1.
C   U011,  - if compu1 = 2, the m-by-m submatrices U011 and U012 of U01,
C   U012,    respectively.
C   compu2 - (optional) scalar indicating whether the orthogonal
C            symplectic transformation matrix U2 is returned, or if U2
C            is not required, as follows:
C            = 0 :  U2 is not required (default);
C            = 1 :  on exit, U21 and U22 contain the corresponding
C                   submatrices of the orthogonal symplectic matrix U2;
C            = 2 :  the orthogonal transformations are accumulated
C                   into U2;
C                   on input, U021, U022 must contain the corresponding
C                   submatrices of an orthogonal symplectic matrix U02;
C                   on exit, U21 and U22 contain the updated submatrices
C                   U21 and U22 of the matrix product U02*U2, with U2
C                   returned for compu2 = 1.
C   U021,  - if compu2 = 2, the m-by-m submatrices U021 and U022 of U02,
C   U022,    respectively.
C   addinf - (optional) scalar indicating whether additional information
C            is required for numerically infinite or unreliable
C            eigenvalues, as follows:
C            = 0 :  additional information is not required (default);
C            = 1 :  additional information is required.
C
C Output parameters:
C   ALPHAR,- the m-vectors of real parts and imaginary parts,
C   ALPHAI   respectively, of each scalar alpha defining an eigenvalue
C            of the pencil aS - bH, with m = n/2.
C            If ALPHAI(j) is zero, then the j-th eigenvalue is real.
C   BETA     the m-vector of the scalars beta that define the
C            eigenvalues of the pencil aS - bH.
C            Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
C            beta = BETA(j) represent the j-th eigenvalue of the pencil
C            aS - bH, in the form lambda = alpha/beta. Since lambda may
C            overflow, the ratios should not, in general, be computed.
C            Due to the skew-Hamiltonian/Hamiltonian structure of the
C            pencil, for every eigenvalue lambda, -lambda is also an
C            eigenvalue, and thus it has only to be saved once in
C            ALPHAR, ALPHAI and BETA.
C            Specifically, only eigenvalues with imaginary parts greater
C            than or equal to zero are stored; their conjugate
C            eigenvalues are not stored. If imaginary parts are zero
C            (i.e., for real eigenvalues), only positive eigenvalues
C            are stored.
C   To     - if job = 1, the computed n-by-n matrix To in (2).
C            If job = 0, the computed n-by-n matrix T obtained just
C            before the application of the periodic QZ algorithm.
C   Zo     - if job = 1, the computed n-by-n matrix Zo in (2).
C            Z(m+1:n,1:m) is not set to 0.
C            If job = 0, the computed n-by-n matrix Z obtained just
C            before the application of the periodic QZ algorithm.
C   Ho     - if job = 1, the computed n-by-n matrix Ho in (2).
C            If job = 0, the computed n-by-n matrix H obtained just
C            before the application of the periodic QZ algorithm.
C   Q1     - if compq1 > 0, an n-by-n matrix containing the computed
C            or updated orthogonal matrix Q1.
C   Q2     - if compq2 > 0, an n-by-n matrix containing the computed
C            or updated orthogonal matrix Q2.
C   U11,   - if compu1 > 0, the m-by-m matrices containing the computed
C   U12      or updated submatrices U11 and U12 of the orthogonal
C            symplectic matrix U1.
C   U21,   - if compu2 > 0, the m-by-m matrices containing the computed
C   U22      or updated submatrices U21 and U22 of the orthogonal
C            symplectic matrix U2.
C            If job = 0, To, Zo, Ho, Q1, Q2, U11, U12, U21, and U22
C            contain the corresponding matrices just before the
C            application of the periodic QZ algorithm.
C   norms  - (optional) if addinf = 1, vector of length 6 containing the
C            Frobenius norms of the six factors of the formal matrix
C            product used by the algorithm. These norms can be used to
C            decide that some eigenvalues should considered as
C            numerically infinite or unreliable.
C   unrel  - (optional) if addinf = 1, the absolute values in unrel are
C            indices of the possibly inaccurate eigenvalues. Positive
C            values correspond to 1-by-1 diagonal blocks (with either
C            real or purely imaginary eigenvalues), and negative values
C            correspond to 2-by-2 diagonal blocks.
C   ublk   - (optional) if addinf = 1, vector containing the s sextuple
C            values corresponding to the 1-by-1 blocks, with real or
C            purely imaginary eigenvalues, followed by t groups of
C            sextuple 2-by-2 matrices corresponding to the 2-by-2
C            blocks.
C   publk  - (optional) if addinf = 1, vector of pointers to the
C            starting location in ublk of each sextuple of blocks
C            corresponding to the unreliable eigenvalues.
C            If addinf = 1, but there are no inaccurate eigenvalues,
C            then unrel, ublk, and publk are void.
C
C Contributors:
C   V. Sima, Research Institute for Informatics, Bucharest, Oct. 2009.
C
C Revisions:
C   V. Sima, Feb. 2010, Nov. 2010, Dec. 2010, July 2011, July 2012,
C   Oct. 2012, July 2013, Mar. 2020.
C   M. Voigt, Jan. 2012.
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
     $                  mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C     .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         COMPQ1, COMPQ2, COMPU1, COMPU2, JOB
      INTEGER           INFO, LDH, LDQ1, LDQ2, LDT, LDU11, LDU12, LDU21,
     $                  LDU22, LDWORK, LDZ, LIWORK, N
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK( : )
      DOUBLE PRECISION, ALLOCATABLE :: IND( : )
      mwPointer         ALPHAI, ALPHAR, BETA, DW, DWORK, H, Q1, Q2, T,
     $                  U11, U12, U21, U22, Z
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           ADDINF, I, ICMPQ1, ICMPQ2, ICMPU1, ICMPU2,
     $                  IDUM( 1 ), IJOB, IP, J, L, LDU, M, M1, N1, N2,
     $                  N3, N4, P
      DOUBLE PRECISION  DUM( 7 ), TEMP
C
C     .. External Subroutines ..
      EXTERNAL          MA02KV, MB04AD
C
C     ..Intrinsic Functions..
      INTRINSIC         INT, MAX, MOD
C
C     Check for proper number of arguments.
C
      IF( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYMPLURV requires at least 2 input arguments.' )
      ELSE IF ( NLHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYMPLURV requires at least 3 output arguments.' )
      END IF
C
C   Z(nxn), H(nxn)(, job, compq1(, Q01), compq2(, Q02),
C                         compu1(, U011, U012), compu2(, U021, U022)
C                         (,addinf)).
C
      N  = mxGetM( PRHS(1) )
      N1 = mxGetN( PRHS(1) )
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'Z must be a real matrix' )
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
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'H must be a real matrix' )
      END IF
      IF ( M1.NE.N ) THEN
         CALL mexErrMsgTxt( 'H must have the same number of rows as Z' )
      END IF
      IF ( N1.NE.N ) THEN
         CALL mexErrMsgTxt( 'H must be a square matrix' )
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
            CALL mexErrMsgTxt( 'COMPQ1 must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'COMPQ1 must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ICMPQ1 = TEMP
         IF ( ICMPQ1.LT.0 .OR. ICMPQ1.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPQ1 has 0, 1, or 2 the only admissible values' )
         END IF
         IP = IP + 1
         IF ( ICMPQ1.EQ.2 )
     $      IP = IP + 1
C
      ELSE
         ICMPQ1 = 0
      END IF
C
      IF ( ICMPQ1.EQ.0 ) THEN
         COMPQ1 = 'N'
      ELSE IF ( ICMPQ1.EQ.1 ) THEN
         COMPQ1 = 'I'
      ELSE
         COMPQ1 = 'U'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'COMPQ2 must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'COMPQ2 must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ICMPQ2 = TEMP
         IF ( ICMPQ2.LT.0 .OR. ICMPQ2.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPQ2 has 0, 1, or 2 the only admissible values' )
         END IF
         IP = IP + 1
         IF ( ICMPQ2.EQ.2 )
     $      IP = IP + 1
C
      ELSE
         ICMPQ2 = 0
      END IF
C
      IF ( ICMPQ2.EQ.0 ) THEN
         COMPQ2 = 'N'
      ELSE IF ( ICMPQ2.EQ.1 ) THEN
         COMPQ2 = 'I'
      ELSE
         COMPQ2 = 'U'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'COMPU1 must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'COMPU1 must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ICMPU1 = TEMP
         IF ( ICMPU1.LT.0 .OR. ICMPU1.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPU1 has 0, 1, or 2 the only admissible values' )
         END IF
         IP = IP + 1
         IF ( ICMPU1.EQ.2 )
     $      IP = IP + 2
C
      ELSE
         ICMPU1 = 0
      END IF
C
      IF ( ICMPU1.EQ.0 ) THEN
         COMPU1 = 'N'
      ELSE IF ( ICMPU1.EQ.1 ) THEN
         COMPU1 = 'I'
      ELSE
         COMPU1 = 'U'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'COMPU2 must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'COMPU2 must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ICMPU2 = TEMP
         IF ( ICMPU2.LT.0 .OR. ICMPU2.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPU2 has 0, 1, or 2 the only admissible values' )
         END IF
         IP = IP + 1
         IF ( ICMPU2.EQ.2 )
     $      IP = IP + 2
C
      ELSE
         ICMPU2 = 0
      END IF
C
      IF ( ICMPU2.EQ.0 ) THEN
         COMPU2 = 'N'
      ELSE IF ( ICMPU2.EQ.1 ) THEN
         COMPU2 = 'I'
      ELSE
         COMPU2 = 'U'
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
      LDT = MAX( 1, N )
      LDZ = LDT
      LDH = LDT
      LDU = MAX( 1, M )
      IF ( ICMPQ1.EQ.0 ) THEN
         LDQ1 = 1
         N1   = 0
      ELSE
         LDQ1 = LDT
         N1   = N
      END IF
      IF ( ICMPQ2.EQ.0 ) THEN
         LDQ2 = 1
         N2   = 0
      ELSE
         LDQ2 = LDT
         N2   = N
      END IF
      IF ( ICMPU1.EQ.0 ) THEN
         LDU11 = 1
         LDU12 = 1
         N3    = 0
      ELSE
         LDU11 = LDU
         LDU12 = LDU
         N3    = M
      END IF
      IF ( ICMPU2.EQ.0 ) THEN
         LDU21 = 1
         LDU22 = 1
         N4    = 0
      ELSE
         LDU21 = LDU
         LDU22 = LDU
         N4    = M
      END IF
      LIWORK = N + 18
C
      CALL MB04AD( JOB, COMPQ1, COMPQ2, COMPU1, COMPU2, N, DUM, LDZ,
     $             DUM, LDH, DUM, LDQ1, DUM, LDQ2, DUM, LDU11, DUM,
     $             LDU12, DUM, LDU21, DUM, LDU22, DUM, LDT, DUM, DUM,
     $             DUM, IDUM, LIWORK, DUM, -1, INFO )
      LDWORK = INT( DUM( 1 ) )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE( IND( M ), IWORK( LIWORK ) )
C
      ALPHAI = mxCalloc(       M,  8 )
      ALPHAR = mxCalloc(       M,  8 )
      BETA   = mxCalloc(       M,  8 )
      DWORK  = mxCalloc(  LDWORK,  8 )
      H      = mxCalloc(   LDH*N,  8 )
      Q1     = mxCalloc(  LDQ1*N1, 8 )
      Q2     = mxCalloc(  LDQ2*N2, 8 )
      T      = mxCalloc(   LDT*N,  8 )
      U11    = mxCalloc( LDU11*N3, 8 )
      U12    = mxCalloc( LDU12*N3, 8 )
      U21    = mxCalloc( LDU21*N4, 8 )
      U22    = mxCalloc( LDU22*N4, 8 )
      Z      = mxCalloc(   LDZ*N,  8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 1 ) ), %VAL( Z ), N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 2 ) ), %VAL( H ), N*N )
C
      IP = 4
      IF ( NRHS.GE.IP .AND. ICMPQ1.EQ.2 ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( IP ) ), %VAL( Q1 ),
     $                          N*N )
         IP = IP + 1
      END IF
C
      IF ( NRHS.GE.IP .AND. ICMPQ2.EQ.2 ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( IP ) ), %VAL( Q2 ),
     $                          N*N )
         IP = IP + 1
      END IF
C
      IF ( NRHS.GE.IP .AND. ICMPU1.EQ.2 ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( IP ) ), %VAL( U11 ),
     $                          M*M )
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( IP ) ), %VAL( U12 ),
     $                          M*M )
         IP = IP + 1
      END IF
C
      IF ( NRHS.GE.IP .AND. ICMPU2.EQ.2 ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( IP ) ), %VAL( U21 ),
     $                          M*M )
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( IP ) ), %VAL( U22 ),
     $                          M*M )
      END IF
C
C Do the actual computations.
C
      CALL MB04AD( JOB, COMPQ1, COMPQ2, COMPU1, COMPU2, N, %VAL( Z ),
     $             LDZ, %VAL( H ), LDH, %VAL( Q1 ), LDQ1, %VAL( Q2 ),
     $             LDQ2, %VAL( U11 ), LDU11, %VAL( U12 ), LDU12,
     $             %VAL( U21 ), LDU21, %VAL( U22 ), LDU22, %VAL( T ),
     $             LDT, %VAL( ALPHAR ), %VAL( ALPHAI ), %VAL( BETA ),
     $             IWORK, LIWORK, %VAL( DWORK ), LDWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 .OR. INFO.EQ.3 ) THEN
         PLHS( 1 ) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( ALPHAR ),
     $                          mxGetPr( PLHS( 1 ) ), M )
C
         PLHS( 2 ) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( ALPHAI ),
     $                          mxGetPr( PLHS( 2 ) ), M )
C
         PLHS( 3 ) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( BETA ),
     $                          mxGetPr( PLHS( 3 ) ), M )
C
         IP = 4
         IF ( NLHS.GE.IP ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( T ), mxGetPr( PLHS( IP ) ),
     $                             N*N )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( Z ), mxGetPr( PLHS( IP ) ),
     $                             N*N )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( H ), mxGetPr( PLHS( IP ) ),
     $                             N*N )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. ICMPQ1.GE.1 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( Q1 ),
     $                             mxGetPr( PLHS( IP ) ), N*N )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. ICMPQ2.GE.1 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( Q2 ),
     $                             mxGetPr( PLHS( IP ) ), N*N )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP+1 .AND. ICMPU1.GE.1 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( U11 ),
     $                             mxGetPr( PLHS( IP ) ), M*M )
            IP = IP + 1
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( U12 ),
     $                             mxGetPr( PLHS( IP ) ), M*M )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP+1 .AND. ICMPU2.GE.1 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( U21 ),
     $                             mxGetPr( PLHS( IP ) ), M*M )
            IP = IP + 1
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( U22 ),
     $                             mxGetPr( PLHS( IP ) ), M*M )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. ADDINF.GT.0 ) THEN
            CALL MA02KV( 6, 2, 1, %VAL( DWORK ), DUM )
            PLHS( IP ) = mxCreateDoubleMatrix( 1, 6, 0 )
            CALL mxCopyReal8ToPtr( DUM, mxGetPr( PLHS( IP ) ), 6 )
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
               P = 6*( M - 2*J ) + 24*J
               DW = mxCalloc( P,  8 )
               CALL MA02KV( P, 8, 1, %VAL( DWORK ), %VAL( DW ) )
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
      DEALLOCATE( IND, IWORK   )
      CALL mxFree( ALPHAI )
      CALL mxFree( ALPHAR )
      CALL mxFree( BETA   )
      CALL mxFree( DWORK  )
      CALL mxFree( H      )
      CALL mxFree( Q1     )
      CALL mxFree( Q2     )
      CALL mxFree( T      )
      CALL mxFree( U11    )
      CALL mxFree( U12    )
      CALL mxFree( U21    )
      CALL mxFree( U22    )
      CALL mxFree( Z      )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 .AND. INFO.NE.3 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB04AD'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of symplURV ***
      END
 
