#include "fintrf.h"
C POLEZERO.F - Gateway function for computing the poles and zeros
C              of a standard or descriptor system, and the Kronecker
C              structure of the system pencil, using SLICOT routines
C              AB08MD, AB08ND, and AG08BD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [(rpoles,ipoles(,betap))(,rzeros,izeros,betaz)(,rnk) ...
C                      (,infz,Kronl,Kronr(,infe,niz,Af,Ef))] ...
C                           = polezero(task,A(,B,C,D)(,E)(,tol,bal))
C   [rpoles,ipoles(,betap)] = polezero(0,A(,E))
C   [rnk]                   = polezero(1,A,B,C,D(,tol,bal))
C   [rnk,infz,Kronl,Kronr(,infe,niz,Af,Ef)] = ...
C                           = polezero(2,A,B,C,D(,E)(,tol,bal))
C   [rzeros,izeros,betaz,rnk,infz,Kronl,Kronr(,infe,niz,Af,Ef)] ...
C                           = polezero(3,A,B,C,D(,E)(,tol,bal))
C   [rpoles,ipoles(,betap),rzeros,izeros,betaz,rnk,infz,Kronl, ...
C                                        Kronr(,infe,niz,Af,Ef)] ...
C                           = polezero(4,A,B,C,D(,E)(,tol,bal))
C
C Purpose:
C   To compute the normal rank, poles, zeros, and the Kronecker
C   structure of the system pencil for a standard or descriptor system.
C
C Input parameters:
C   task   - integer specifying the computations to be performed.
C            = 0 :  compute the poles of the system;
C            = 1 :  compute the normal rank of the transfer-function
C                   matrix of a standard system;
C            = 2 :  compute the normal rank and Kronecker structure of
C                   the system;
C            = 3 :  compute the normal rank, Kronecker structure, and
C                   zeros of the system;
C            = 4 :  compute the normal rank, Kronecker structure, poles,
C                   and zeros of the system.
C   A      - the l-by-n state dynamics matrix A.
C            If task = 0, task = 1, or task = 4, it is assumed that
C            l = n.
C   B      - the l-by-m input/state matrix B.
C            If task = 0, B must not be specified as an input parameter.
C   C      - the p-by-n state/output matrix C.
C            If task = 0, C must not be specified as an input parameter.
C   D      - the p-by-m input/output matrix D.
C            If task = 0, D must not be specified as an input parameter.
C   E      - (optional) the l-by-n descriptor matrix E.
C            For a standard system, one must have l = n, and matrix E
C            may not be specified as an input parameter, or it may be
C            set to either an empty matrix or an identity matrix of
C            order n.
C            If task > 1, l = n = 1, and the number of input parameters
C            is 6 or larger, then the 6-th parameter is taken as E,
C            not tol.
C   tol    - (optional) if task > 0, real scalar containing the
C            tolerance to be used in rank decisions to determine the
C            effective rank, which is defined as the order of the
C            largest leading (or trailing) triangular submatrix in the
C            QR (or RQ) factorization with column (or row) pivoting
C            whose estimated condition number is less than 1/tol.
C            If tol <= 0, then default tolerances are used instead,
C            defined in terms of the size of the system matrix and eps,
C            where eps is the machine precision.
C            Default: tol = 0.
C            If task = 0, tol must not be specified as an input
C            parameter.
C   bal    - (optional) if task > 0, integer indicating whether the
C            system should be balanced (scaled).
C            = 0 :  use balancing;
C            = 1 :  do not use balancing.
C            Default: bal = 0.
C            If task = 0, bal must not be specified as an input
C            parameter.
C
C Output parameters:
C   rpoles - (optional) if task = 0 or task = 4, an n-vector containing
C            the real parts of the (alpha part, if E is general, of the)
C            system poles.
C            Otherwise, rpoles must not be specified as an output
C            parameter.
C   ipoles - (optional) if task = 0 or task = 4, an n-vector containing
C            the imaginary parts of the (alpha part, if E is general, of
C            the) system poles.
C            Otherwise, ipoles must not be specified as an output
C            parameter.
C   betap  - (optional) if E is general, and task = 0 or task = 4, an
C            n-vector containing the beta part of the system poles.
C            Otherwise, betap must not be specified as an output
C            parameter.
C   rzeros - (optional) if task >= 3, a vector containing the real
C            parts of the alpha part of the system zeros.
C            Otherwise, rzeros must not be specified as an output
C            parameter.
C   izeros - (optional) if task >= 3, a vector containing the imaginary
C            parts of the alpha part of the system zeros.
C            Otherwise, izeros must not be specified as an output
C            parameter.
C   betaz  - (optional) if task >= 3, a vector containing the beta part
C            of the system zeros.
C            Otherwise, betaz must not be specified as an output
C            parameter.
C   rnk    - (optional) if task > 0, the normal rank of the system
C            pencil.
C            Otherwise, rnk must not be specified as an output
C            parameter.
C   infz   - (optional) if task >= 2, integer vector of length dinfz
C            containing information on the infinite elementary divisors
C            as follows: the system has infz(i) infinite elementary
C            divisors of degree i (in the Smith form), where
C            i = 1,2,...,dinfz.
C            Otherwise, infz must not be specified as an output
C            parameter.
C   Kronl  - (optional) if task >= 2, integer vector containing the left
C            Kronecker (row) indices.
C            Otherwise, Kronl must not be specified as an output
C            parameter.
C   Kronr  - (optional) if task >= 2, integer vector containing the
C            right Kronecker (column) indices.
C            Otherwise, Kronr must not be specified as an output
C            parameter.
C   infe   - (optional) if task >= 2, integer vector containing the
C            multiplicities of infinite eigenvalues.
C            Otherwise, infe must not be specified as an output
C            parameter.
C   niz    - (optional) if task >= 2, the number of infinite zeros.
C            Otherwise, niz must not be specified as an output
C            parameter.
C   Af     - (optional) if task >= 2, the nfz-by-nfz matrix Af of the
C            reduced pencil (see Method). If task >= 3, Af is in a real
C            Schur form.
C            Otherwise, Af must not be specified as an output
C            parameter.
C   Ef     - (optional) if task >= 2, the nfz-by-nfz matrix Ef of the
C            reduced pencil (see Method). If task >= 3, Ef is in an
C            upper triangular form.
C            Otherwise, Ef must not be specified as an output
C            parameter.
C
C Method:
C   If task = 0, the (generalized) eigenvalues of the system (pencil)
C   are computed.
C
C   If task = 1, the (n+p)-by-(m+n) compound matrix (B  A) is reduced
C                                                   (D  C)
C
C   to one with the same invariant zeros and with D of full row rank.
C   The normal rank of the transfer-function matrix is the rank of D.
C
C   If task >= 2 and identity E, the gateway extracts from the system
C   matrix of a state-space system (A,B,C,D) a regular pencil
C   Af - lambda*Ef which has the invariant zeros of the system as
C   generalized eigenvalues, as follows:
C
C        (a) construct the (n+p)-by-(m+n) compound matrix (B  A);
C                                                         (D  C)
C
C        (b) reduce the above system to one with the same invariant
C            zeros and with D of full row rank;
C
C        (c) pertranspose the system;
C
C        (d) reduce the system to one with the same invariant zeros and
C            with D square invertible;
C
C        (e) perform a unitary transformation on the columns of
C            (A - lambda*I  B) in order to reduce it to
C            (      C       D)
C
C            (Af - lambda*Ef  X)
C            (                 ), with Y and Ef square invertible;
C            (      0         Y)
C
C        (f) compute the right and left Kronecker indices of the system
C            (A,B,C,D), which together with the orders of the infinite
C            zeros (determined by steps (a) - (e)) constitute the
C            complete set of structural invariants under strict
C            equivalence transformations of a linear system.
C
C   If task >= 2 and general E, the gateway extracts from the system
C   matrix of a descriptor system (A-lambda*E,B,C,D) a regular pencil
C   Af-lambda*Ef which has the finite zeros of the system as generalized
C   eigenvalues. The procedure has the following computational steps:
C
C        (a) construct the (l+p)-by-(m+n) system pencil
C
C             S(lambda) = ( B  A )-lambda*( 0  E );
C                         ( D  C )        ( 0  0 )
C
C        (b) reduce S(lambda) to S1(lambda) with the same finite
C            zeros and right Kronecker structure but with E upper
C            triangular and nonsingular;
C
C        (c) reduce S1(lambda) to S2(lambda) with the same finite zeros
C            and right Kronecker structure but with D of full row rank;
C
C        (d) reduce S2(lambda) to S3(lambda) with the same finite zeros
C            and with D square invertible;
C
C        (e) perform a unitary transformation on the columns of
C
C            S3(lambda) = (A-lambda*E   B) in order to reduce it to
C                         (     C       D)
C
C            (Af-lambda*Ef   X), with Y and Ef square invertible;
C            (     0         Y)
C
C        (f) compute the right and left Kronecker indices of the system
C            matrix, which together with the multiplicities of the
C            finite and infinite eigenvalues constitute the
C            complete set of structural invariants under strict
C            equivalence transformations of a linear system.
C
C Comments
C   1. If only the normal rank of the system pencil of a standard system
C      is desired, using task = 1 will ensure the maximum efficiency.
C   2. If E is general, the poles and zeros are returned as generalized
C      eigenvalues, with complex numerators alpha and real denominators
C      beta, but the ratios are not computed.
C   3. If Af and Ef are not needed, they should not be specified as
C      output arguments for task >= 3, since then the zeros are computed
C      slightly faster.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Feb. 2005,
C   Jan. 2009, Mar. 2009, Apr. 2009, July 2011, June 2012, Feb. 2015.
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
C .. Scalar parameters used by SLICOT and LAPACK subroutines ..
      CHARACTER         EQUIL
      INTEGER           DINFZ, INFO, INFO1, INFOV, INFOV1, L, LDA, LDAF,
     $                  LDB, LDBF, LDC, LDD, LDE, LDWORK, M, N, NFZ,
     $                  NINFE, NIZ, NKROL, NKROR, NRANK, P, SDIM
      DOUBLE PRECISION  TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER, ALLOCATABLE ::          INFE(:), INFZ(:), IWORK(:),
     $                                 KRONL(:), KRONR(:)
      DOUBLE PRECISION, ALLOCATABLE :: DW(:)
      mwPointer                        A, AF, AS, B, BETAP, BETAZ, BF,
     $                                 C, D, DWORK, E, ES, PI, PR, ZI,
     $                                 ZR
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           BWORK(1), DESCRP, OK, WITHS
      INTEGER           BAL, I, IP, IPE, ITASK, J, LE, LINP, LIWORK,
     $                  LOUT, M1N, N1N, NE
      DOUBLE PRECISION  DUM(1), TEMP
C
C     .. External Functions ..
      LOGICAL           SELECT
      DOUBLE PRECISION  DASUM
      EXTERNAL          DASUM, SELECT
C
C .. External subroutines ..
      EXTERNAL          AB08MD, AB08ND, AG08BD, DAXPY, DCOPY, DGEES,
     $                  DGEEV, DGGES, DGGEV, DLACPY
C
C .. Intrinsic functions ..
      INTRINSIC         INT, MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'POLEZERO requires at least 2 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'POLEZERO requires at least 1 output argument' )
      END IF
C
C   task, A(lxn), B(lxm), C(p,n), D(p,m), E(l,n), tol, bal.
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
      IF ( ITASK.LT.0 .OR. ITASK.GT.4 ) THEN
         CALL mexErrMsgTxt
     $    ( 'ITASK has 0, 1, 2, 3 or 4 the only admissible values')
      END IF
C
      IF ( ITASK.EQ.0 ) THEN
         LINP = 2
         LOUT = 2
      ELSE
         LINP = 5
         IF ( ITASK.EQ.1 ) THEN
            LOUT = 1
         ELSE IF ( ITASK.EQ.2 ) THEN
            LOUT = 2*ITASK
         ELSE
            LOUT = 2*ITASK + 1
         END IF
      END IF
C
      IF ( NRHS.LT.LINP ) THEN
         WRITE( TEXT, '('' POLEZERO requires at least '',I2,
     $                  '' input arguments'')' ) LINP
         CALL mexErrMsgTxt( TEXT )
      ELSE IF ( NLHS.LT.LOUT ) THEN
         WRITE( TEXT, '('' POLEZERO requires at least '',I2,
     $                  '' output arguments'')' ) LOUT
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      IF ( ITASK.EQ.3 ) THEN
         WITHS = NLHS.GE.10
      ELSEIF ( ITASK.EQ.4 ) THEN
         WITHS = NLHS.GE.13
      ELSE
         WITHS = .FALSE.
      END IF
C
      L = mxGetM( PRHS(2) )
      N = mxGetN( PRHS(2) )
C
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF ( ( ITASK.LE.1 .OR. ITASK.EQ.4 ) .AND. N.NE.L ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
C
      IF ( ITASK.GT.0 ) THEN
         M = mxGetN( PRHS(3) )
         P = mxGetM( PRHS(4) )
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(3) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'B must be a real matrix' )
         END IF
         IF ( mxGetM( PRHS(3) ).NE.L ) THEN
            CALL mexErrMsgTxt
     $         ( 'B must have the same number of rows as A' )
         END IF
C
         IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(4) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'C must be a real matrix' )
         END IF
         IF ( mxGetN( PRHS(4) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $         ( 'C must have the same number of columns as A' )
         END IF
C
         IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(5) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'D must be a real matrix' )
         END IF
         IF ( mxGetM( PRHS(5) ).NE.P ) THEN
            CALL mexErrMsgTxt
     $         ( 'D must have the same number of rows as C' )
         END IF
         IF ( mxGetN( PRHS(5) ).NE.M ) THEN
            CALL mexErrMsgTxt
     $         ( 'D must have the same number of columns as B' )
         END IF
         IP = 6
      ELSE
         IP = 3
      END IF
C
      IPE = IP
      LDA = MAX( 1, L )
C
      DESCRP = .FALSE.
      IF ( ITASK.NE.1 .AND. NRHS.GE.IP ) THEN
         LE = mxGetM( PRHS(IP) )
         NE = mxGetN( PRHS(IP) )
         OK = mxIsNumeric( PRHS(IP) ).NE.0 .AND.
     $        mxIsComplex( PRHS(IP) ).NE.1
         IF ( OK ) THEN
            IF ( LE*NE.EQ.0 ) THEN
               IP = IP + 1
            ELSE IF ( LE*NE.EQ.1 ) THEN
               IF ( L*N.EQ.1 ) THEN
                  CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), DUM,
     $                                   1 )
                  DESCRP = DUM(1).NE.ONE
                  IF ( DESCRP ) THEN
                     E = mxCalloc( LDA*N, 8 )
                     CALL DCOPY( 1, DUM, 1, %VAL(E), 1 )
                     IP = IP + 1
                  END IF
               END IF
            ELSE
               IF ( LE.NE.L ) THEN
                  CALL mexErrMsgTxt
     $               ( 'E must have the same number of rows as A' )
               END IF
               IF ( NE.NE.N ) THEN
                  CALL mexErrMsgTxt
     $               ( 'E must have the same number of columns as A' )
               END IF
               IF ( L.EQ.N ) THEN
                  E = mxCalloc( LDA*N, 8 )
                  CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ),
     $                                   %VAL(E), N*N )
                  DUM(1) = ONE
                  CALL DAXPY( N, -ONE, DUM, 0, %VAL(E), N+1 )
                  IF ( DASUM( N*N, %VAL(E), 1 ).NE.ZERO ) THEN
                     DESCRP = .TRUE.
                     CALL DAXPY( N, ONE, DUM, 0, %VAL(E), N+1 )
                  ELSE
                     CALL mxFree( E )
                  END IF
               ELSE
                  DESCRP = .TRUE.
               END IF
               IP = IP + 1
            END IF
         ELSE IF ( LE*NE.NE.1 .OR. L*N.EQ.1 ) THEN
            CALL mexErrMsgTxt( 'E must be a real matrix' )
         ELSE
            CALL mexErrMsgTxt( 'TOL must be a real scalar' )
         END IF
      END IF
C
      IF ( ITASK.GT.0 ) THEN
         IF ( NRHS.GE.IP ) THEN
C
C   tol
C
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'TOL must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'TOL must be a real scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TOL, 1 )
            IP = IP + 1
         ELSE
            TOL = ZERO
         END IF
C
         IF ( NRHS.GE.IP ) THEN
C
C   bal
C
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'BAL must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'BAL must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            BAL = TEMP
            IF ( BAL.LT.0 .OR. BAL.GT.1 ) THEN
               CALL mexErrMsgTxt
     $           ( 'BAL has 0 or 1 the only admissible values' )
            END IF
         ELSE
            BAL = 0
         END IF
C
         IF ( BAL.EQ.0 ) THEN
            EQUIL = 'S'
         ELSE
            EQUIL = 'N'
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C Use a value for LDWORK large enough for enabling calls of block
C algorithms in AB08MD, AB08ND, and AG08BD, depending on task.
C Optimal workspace is also used in DGEES, DGEEV, DGGES, or DGGEV, if
C task = 0 or task >= 3.
C
      IF ( ITASK.GT.0 ) THEN
         LDB = LDA
         LDC = MAX( 1, P )
         LDD = LDC
         IF ( .NOT.DESCRP ) THEN
            LDAF = MAX( 1, N + M )
            LDBF = MAX( 1, N + P )
         END IF
      END IF
      IF ( DESCRP )
     $   LDE = LDA
      IF ( ITASK.EQ.0 .OR. ITASK.EQ.4 ) THEN
         LIWORK = 0
         IF ( DESCRP ) THEN
            CALL DGGEV( 'No vectors', 'No vectors', N, DUM, LDA, DUM,
     $                  LDA, DUM, DUM, DUM, DUM, 1, DUM, 1, DUM, -1,
     $                  INFOV )
         ELSE
            CALL DGEEV( 'No vectors', 'No vectors', N, DUM, LDA, DUM,
     $                  DUM, DUM, 1, DUM, 1, DUM, -1, INFOV1 )
         END IF
         LDWORK = INT( DUM(1) )
         IF ( ITASK.EQ.4 ) THEN
            IF ( WITHS ) THEN
               CALL DGGES( 'No vectors', 'No vectors', 'Not ordered',
     $                     SELECT, N, DUM, LDA, DUM, LDA, SDIM, DUM,
     $                     DUM, DUM, DUM, 1, DUM, 1, DUM, -1, BWORK,
     $                     INFO )
            ELSE
               CALL DGGEV( 'No vectors', 'No vectors', N, DUM, LDA, DUM,
     $                     LDA, DUM, DUM, DUM, DUM, 1, DUM, 1, DUM, -1,
     $                     INFOV )
            END IF
            LDWORK = MAX( LDWORK, INT( DUM(1) ) )
         END IF
      ELSE IF ( ITASK.EQ.1 ) THEN
         LIWORK = 2*N + MAX( M, P ) + 1
         CALL AB08MD( EQUIL, N, M, P, DUM, LDA, DUM, LDB, DUM, LDC, DUM,
     $                LDD, NRANK, TOL, IWORK, DUM, -1, INFO )
         LDWORK = INT( DUM(1) )
      ELSE IF ( ITASK.EQ.3 ) THEN
         IF ( WITHS ) THEN
            CALL DGGES( 'No vectors', 'No vectors', 'Not ordered',
     $                  SELECT, MIN( L, N ), DUM, LDA, DUM, LDA, SDIM,
     $                  DUM, DUM, DUM, DUM, 1, DUM, 1, DUM, -1, BWORK,
     $                  INFO )
         ELSE
            CALL DGGEV( 'No vectors', 'No vectors', MIN( L, N ), DUM,
     $                  LDA, DUM, LDA, DUM, DUM, DUM, DUM, 1, DUM, 1,
     $                  DUM, -1, INFOV )
         END IF
         LDWORK = INT( DUM(1) )
      ELSE
         LDWORK = 0
      END IF
      IF ( ITASK.GE.2 ) THEN
         IF ( DESCRP ) THEN
            LIWORK = N + MAX( 1, M )
            CALL AG08BD( EQUIL, L, N, M, P, DUM, LDA, DUM, LDE, DUM,
     $                   LDB, DUM, LDC, DUM, LDD, NFZ, NRANK, NIZ,
     $                   DINFZ, NKROR, NINFE, NKROL, INFZ, KRONR, INFE,
     $                   KRONL, TOL, IWORK, DUM, -1, INFO )
            IF ( BAL.EQ.0 )
     $         LDWORK = MAX( LDWORK, 4*( L + N ) )
         ELSE
            LIWORK = MAX( M, P )
            CALL AB08ND( EQUIL, N, M, P, DUM, LDA, DUM, LDB, DUM, LDC,
     $                   DUM, LDD, NFZ, NRANK, DINFZ, NKROR, NKROL,
     $                   INFZ, KRONR, KRONL, DUM, LDAF, DUM, LDBF, TOL,
     $                   IWORK, DUM, -1, INFO )
         END IF
         LDWORK = MAX( LDWORK, INT( DUM(1) ), MAX( N, M, P ) + 1 )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( IWORK( LIWORK ) )
      A     = mxCalloc(  LDA*N, 8 )
      DWORK = mxCalloc( LDWORK, 8 )
      IF ( ITASK.GT.0 ) THEN
         B = mxCalloc( LDB*M, 8 )
         C = mxCalloc( LDC*N, 8 )
         D = mxCalloc( LDD*M, 8 )
         IF ( DESCRP ) THEN
            ALLOCATE ( INFE( 1 + MIN( L + P, N + M ) ), INFZ( N + 1 ),
     $                 KRONL( L + P + 1 ), KRONR( N + M + 1 ) )
            ALLOCATE ( DW( 1 + MAX( L + P, N + M ) ) )
         ELSE IF ( ITASK.GT.1 ) THEN
            ALLOCATE ( INFE( 2*N + MIN( P, M ) ), INFZ( N + 1 ),
     $                 KRONL( MAX( N, P ) + 1 ),
     $                 KRONR( MAX( N, M ) + 1 ) )
            ALLOCATE ( DW( 1 + MAX( N, P, M ) ) )
            AF = mxCalloc( LDAF*( N + MIN( P, M ) ), 8 )
            BF = mxCalloc( LDBF*( N + M ),           8 )
         END IF
      END IF
      IF ( DESCRP .AND. L.NE.N )
     $   E = mxCalloc( LDE*N, 8 )
C
      IF ( ITASK.EQ.0 .OR. ITASK.EQ.4 ) THEN
         PI = mxCalloc( N, 8 )
         PR = mxCalloc( N, 8 )
         IF ( DESCRP )
     $      BETAP = mxCalloc( N, 8 )
      END IF
C
      IF ( ITASK.GE.3 ) THEN
         BETAZ = mxCalloc( N, 8 )
         ZI    = mxCalloc( N, 8 )
         ZR    = mxCalloc( N, 8 )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL(A), L*N )
      IF ( ITASK.GT.0 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL(B), L*M )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL(C), P*N )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), %VAL(D), P*M )
      END IF
      IF ( DESCRP .AND. L.NE.N )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(IPE) ), %VAL(E), L*N )
C
C Do the actual computations.
C
      INFO   = 0
      INFO1  = 0
      INFOV  = 0
      INFOV1 = 0
      IF ( ITASK.EQ.0 .OR. ITASK.EQ.4 ) THEN
         IF ( ITASK.EQ.4 ) THEN
            AS = mxCalloc( LDA*N, 8 )
            CALL DLACPY( 'Full', L, N, %VAL(A), LDA, %VAL(AS), LDA )
            IF ( DESCRP )  THEN
               ES = mxCalloc( LDE*N, 8 )
               CALL DLACPY( 'Full', L, N, %VAL(E), LDE, %VAL(ES), LDE )
               CALL DGGEV( 'No vectors', 'No vectors', N, %VAL(AS), LDA,
     $                     %VAL(ES), LDE, %VAL(PR), %VAL(PI),
     $                     %VAL(BETAP), DUM, 1, DUM, 1, %VAL(DWORK),
     $                     LDWORK, INFOV )
               CALL mxFree( ES )
            ELSE
               CALL DGEEV( 'No vectors', 'No vectors', N, %VAL(AS), LDA,
     $                     %VAL(PR), %VAL(PI), DUM, 1, DUM, 1,
     $                     %VAL(DWORK), LDWORK, INFOV1 )
            END IF
            CALL mxFree( AS )
         ELSE
            IF ( DESCRP )  THEN
               CALL DGGEV( 'No vectors', 'No vectors', N, %VAL(A), LDA,
     $                     %VAL(E), LDE, %VAL(PR), %VAL(PI),
     $                     %VAL(BETAP), DUM, 1, DUM, 1, %VAL(DWORK),
     $                     LDWORK, INFOV )
            ELSE
               CALL DGEEV( 'No vectors', 'No vectors', N, %VAL(A), LDA,
     $                     %VAL(PR), %VAL(PI), DUM, 1, DUM, 1,
     $                     %VAL(DWORK), LDWORK, INFOV1 )
            END IF
         END IF
         IF ( INFOV.NE.0 .OR. INFOV1.NE.0 )
     $      GO TO 70
C
      ELSE IF ( ITASK.EQ.1 ) THEN
         CALL AB08MD( EQUIL, N, M, P, %VAL(A), LDA, %VAL(B), LDB,
     $                %VAL(C), LDC, %VAL(D), LDD, NRANK, TOL, IWORK,
     $                %VAL(DWORK), LDWORK, INFO )
         NRANK = N + NRANK
      END IF
C
      IF ( ITASK.GE.2 ) THEN
         IF ( DESCRP )  THEN
            CALL AG08BD( EQUIL, L, N, M, P, %VAL(A), LDA, %VAL(E), LDE,
     $                   %VAL(B), LDB, %VAL(C), LDC, %VAL(D), LDD, NFZ,
     $                   NRANK, NIZ, DINFZ, NKROR, NINFE, NKROL, INFZ,
     $                   KRONR, INFE, KRONL, TOL, IWORK, %VAL(DWORK),
     $                   LDWORK, INFO )
            IF ( ITASK.GE.3 ) THEN
               IF ( WITHS ) THEN
                  CALL DGGES( 'No vectors', 'No vectors', 'Not ordered',
     $                        SELECT, NFZ, %VAL(A), LDA, %VAL(E), LDE,
     $                        SDIM, %VAL(ZR), %VAL(ZI), %VAL(BETAZ),
     $                        DUM, 1, DUM, 1, %VAL(DWORK), LDWORK,
     $                        BWORK, INFO )
               ELSE
                  CALL DGGEV( 'No vectors', 'No vectors', NFZ, %VAL(A),
     $                        LDA, %VAL(E), LDE, %VAL(ZR), %VAL(ZI),
     $                        %VAL(BETAZ), DUM, 1, DUM, 1, %VAL(DWORK),
     $                        LDWORK, INFOV )
               END IF
            END IF
C
         ELSE
            CALL AB08ND( EQUIL, N, M, P, %VAL(A), LDA, %VAL(B), LDB,
     $                   %VAL(C), LDC, %VAL(D), LDD, NFZ, NRANK, DINFZ,
     $                   NKROR, NKROL, INFZ, KRONR, KRONL, %VAL(AF),
     $                   LDAF, %VAL(BF), LDBF, TOL, IWORK, %VAL(DWORK),
     $                   LDWORK, INFO )
C
            IF ( INFO.EQ.0 ) THEN
               NIZ   = 0
               NINFE = 0
C
               DO 5 I = 1, DINFZ
                  IF ( INFZ(I).GT.0 ) THEN
                     NIZ = NIZ + INFZ(I)*I
                     NINFE = NINFE + INFZ(I)
                  END IF
    5          CONTINUE
C
               NINFE = NRANK - NINFE
C
               DO 10 I = 1, NINFE
                  INFE(I) = 1
   10          CONTINUE
C
               DO 20 I = 1, DINFZ
                  DO 15 II = NINFE + 1, NINFE + INFZ(I)
                     INFE(II) = I + 1
   15             CONTINUE
                  NINFE = NINFE + INFZ(I)
   20          CONTINUE
C
            END IF
C
            IF ( INFO.EQ.0 .AND. ITASK.GE.3 ) THEN
               DUM(1) = ONE
               IF ( WITHS ) THEN
                  IF ( NRANK.EQ.0 ) THEN
                     CALL DGEES( 'No vectors', 'Not ordered', SELECT,
     $                           NFZ, %VAL(AF), LDAF, SDIM, %VAL(ZR),
     $                           %VAL(ZI), DUM, 1, %VAL(DWORK), LDWORK,
     $                           BWORK, INFO1 )
                     CALL DCOPY( NFZ, DUM, 0, %VAL(BETAZ), 1 )
                  ELSE
                     CALL DGGES( 'No vectors', 'No vectors',
     $                           'Not ordered', SELECT, NFZ, %VAL(AF),
     $                           LDAF, %VAL(BF), LDBF, SDIM, %VAL(ZR),
     $                           %VAL(ZI), %VAL(BETAZ), DUM, 1, DUM, 1,
     $                           %VAL(DWORK), LDWORK, BWORK, INFO )
                  END IF
               ELSE
                  IF ( NRANK.EQ.0 ) THEN
                     CALL DGEEV( 'No vectors', 'No vectors', NFZ,
     $                           %VAL(AF), LDAF, %VAL(ZR), %VAL(ZI),
     $                           DUM, 1, DUM, 1, %VAL(DWORK), LDWORK,
     $                           INFOV1 )
                     CALL DCOPY( NFZ, DUM, 0, %VAL(BETAZ), 1 )
                  ELSE
                     CALL DGGEV( 'No vectors', 'No vectors', NFZ,
     $                           %VAL(AF), LDAF, %VAL(BF), LDBF,
     $                           %VAL(ZR), %VAL(ZI), %VAL(BETAZ), DUM,
     $                           1, DUM, 1, %VAL(DWORK), LDWORK, INFOV )
                  END IF
               END IF
            END IF
            NRANK = N + NRANK
         END IF
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         IP = 1
         IF ( ITASK.EQ.0 .OR. ITASK.EQ.4 ) THEN
            M1N = MIN( 1, N )
            N1N = N*M1N
            PLHS(IP) = mxCreateDoubleMatrix( N, M1N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(PR), mxGetPr( PLHS(IP) ),
     $                             N1N )
            IP = IP + 1
            PLHS(IP) = mxCreateDoubleMatrix( N, M1N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(PI), mxGetPr( PLHS(IP) ),
     $                             N1N )
            IP = IP + 1
            IF ( DESCRP ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( N, M1N, 0 )
               CALL mxCopyReal8ToPtr( %VAL(BETAP),
     $                                mxGetPr( PLHS(IP) ), N1N )
               IP = IP + 1
            END IF
         END IF
C
         IF ( ITASK.GE.3 ) THEN
            M1N = MIN( 1, NFZ )
            N1N = NFZ*M1N
            PLHS(IP) = mxCreateDoubleMatrix( NFZ, M1N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(ZR), mxGetPr( PLHS(IP) ),
     $                             N1N )
            IP = IP + 1
            PLHS(IP) = mxCreateDoubleMatrix( NFZ, M1N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(ZI), mxGetPr( PLHS(IP) ),
     $                             N1N )
            IP = IP + 1
            PLHS(IP) = mxCreateDoubleMatrix( NFZ, M1N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(BETAZ), mxGetPr( PLHS(IP) ),
     $                             N1N )
            IP = IP + 1
         END IF
C
         IF ( ITASK.GE.1 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
            TEMP = NRANK
            CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS(IP) ), 1 )
            IP = IP + 1
         END IF
C
         IF ( ITASK.GE.2 ) THEN
            M1N = MIN( 1, DINFZ )
            N1N = DINFZ*M1N
C
            DO 30 I = 1, N1N
               DW(I) = INFZ(I)
   30       CONTINUE
C
            PLHS(IP) = mxCreateDoubleMatrix( DINFZ, M1N, 0 )
            CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(IP) ), N1N )
            IP = IP + 1
            M1N = MIN( 1, NKROL )
            N1N = NKROL*M1N
C
            DO 40 I = 1, N1N
               DW(I) = KRONL(I)
   40       CONTINUE
C
            PLHS(IP) = mxCreateDoubleMatrix( NKROL, M1N, 0 )
            CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(IP) ), N1N )
            IP = IP + 1
            M1N = MIN( 1, NKROR )
            N1N = NKROR*M1N
C
            DO 50 I = 1, N1N
               DW(I) = KRONR(I)
   50       CONTINUE
C
            PLHS(IP) = mxCreateDoubleMatrix( NKROR, M1N, 0 )
            CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(IP) ), N1N )
            IP = IP + 1
            IF ( NLHS.GE.IP ) THEN
               M1N = MIN( 1, NINFE )
               N1N = NINFE*M1N
C
               DO 60 I = 1, N1N
                  DW(I) = INFE(I)
   60          CONTINUE
C
               PLHS(IP) = mxCreateDoubleMatrix( NINFE, M1N, 0 )
               CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(IP) ), N1N )
               IP = IP + 1
            END IF
            IF ( NLHS.GE.IP ) THEN
               TEMP = NIZ
               PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
               CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS(IP) ), 1 )
               IP = IP + 1
            END IF
C
            IF ( NLHS.GE.IP ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( NFZ, NFZ, 0 )
               IF ( NLHS.GE.IP + 1 )
     $            PLHS(IP+1) = mxCreateDoubleMatrix( NFZ, NFZ, 0 )
               IF ( DESCRP ) THEN
                  IF( NFZ.GT.0 .AND. NFZ.LT.LDA ) THEN
                     CALL DLACPY( 'Full', NFZ, NFZ, %VAL(A), LDA,
     $                            %VAL(A), NFZ )
                     IF ( NLHS.GE.IP + 1 )
     $                  CALL DLACPY( 'Full', NFZ, NFZ, %VAL(E), LDE,
     $                               %VAL(E), NFZ )
                  END IF
                  CALL mxCopyReal8ToPtr( %VAL(A),
     $                                  mxGetPr( PLHS(IP) ),
     $                                  NFZ*NFZ )
                  IF ( NLHS.GE.IP + 1 )
     $               CALL mxCopyReal8ToPtr( %VAL(E),
     $                                      mxGetPr( PLHS(IP+1) ),
     $                                      NFZ*NFZ )
               ELSE
                  IF( NFZ.GT.0 .AND. NFZ.LT.LDAF ) THEN
                     CALL DLACPY( 'Full', NFZ, NFZ, %VAL(AF), LDAF,
     $                            %VAL(AF), NFZ )
                     IF ( NLHS.GE.IP + 1 ) THEN
                        CALL DLACPY( 'Full', NFZ, NFZ, %VAL(BF), LDBF,
     $                               %VAL(BF), NFZ )
                     END IF
                  END IF
                  CALL mxCopyReal8ToPtr( %VAL(AF),
     $                                   mxGetPr( PLHS(IP) ),
     $                                   NFZ*NFZ )
                  IF ( NLHS.GE.IP + 1 )
     $               CALL mxCopyReal8ToPtr( %VAL(BF),
     $                                      mxGetPr( PLHS(IP+1) ),
     $                                      NFZ*NFZ )
               END IF
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
   70 CONTINUE
      DEALLOCATE ( IWORK )
      CALL mxFree( A     )
      CALL mxFree( DWORK )
      IF ( ITASK.GT.0 ) THEN
         CALL mxFree( B )
         CALL mxFree( C )
         CALL mxFree( D )
         IF ( DESCRP ) THEN
            DEALLOCATE ( DW, INFE, INFZ, KRONL, KRONR )
         ELSE IF ( ITASK.GT.1 ) THEN
            DEALLOCATE ( DW, INFE, INFZ, KRONL, KRONR )
            CALL mxFree( AF )
            CALL mxFree( BF )
         END IF
      END IF
      IF ( DESCRP )
     $   CALL mxFree( E )
C
      IF ( ITASK.EQ.0 .OR. ITASK.EQ.4 ) THEN
         CALL mxFree( PI )
         CALL mxFree( PR )
         IF ( DESCRP )
     $      CALL mxFree( BETAP )
      END IF
C
      IF ( ITASK.GE.3 ) THEN
         CALL mxFree( BETAZ )
         CALL mxFree( ZI    )
         CALL mxFree( ZR    )
      END IF
C
C Error and warning handling.
C
      IF ( MAX( INFO, INFO1, INFOV, INFOV1 ).NE.0 ) THEN
         IF ( INFO.NE.0 )  THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DGGES'')'
     $           ) INFO
         ELSE IF ( INFOV.NE.0 )  THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DGGEV'')'
     $           ) INFOV
         ELSE IF ( INFO1.NE.0 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DGEES'')'
     $           ) INFO1
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DGEEV'')'
     $           ) INFOV1
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of POLEZERO ***
      END
