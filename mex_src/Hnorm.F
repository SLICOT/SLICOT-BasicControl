#include "fintrf.h"
C HNORM.F  - Gateway function for computing various system norms
C            (Hankel norm, H2 norm) and complex stability radius of a
C            continuous-time or discrete-time system in standard form,
C            using SLICOT routines AB13AD, AB13BD, AB13ED and AB13FD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [Hnorm(,ns,Ao,Bo,Co(,Do)(,hsv))]
C                              = Hnorm(job,A,B,C(,D)(,dico,equil
C                                      (,alpha)(,tol)))
C   [low(,high))]              = Hnorm(job,A(,tol))
C
C   [Hknorm(,ns,Ao,Bo,Co,hsv)] = Hnorm(1,A,B,C(,dico,equil,alpha))
C   [Hnorm(,nq,Ao,Bo,Co,Do)]   = Hnorm(2,A,B,C,D(,dico,jobn,tol))
C   [low(,high)]               = Hnorm(3,A(,tol))
C   [beta(,omega)]             = Hnorm(4,A(,tol))
C
C Purpose:
C   1) To compute the Hankel norm, H2 or L2 norm of the continuous-time
C   or discrete-time system in standard form
C                                     -1
C        G(lambda) = C*( lambda*I - A ) *B + D .
C
C   For computing the Hankel norm, or for computing the H2/L2 norm for
C   a continuous-time system, D is assumed 0. For computing the H2 or
C   L2 norm, the transfer-function matrix G must not have poles on the
C   imaginary axis, for a continuous-time system, or on the unit circle,
C   for a discrete-time system. If the H2-norm is computed, the system
C   must be stable.
C
C   2) To estimate beta(A), the 2-norm distance from a real matrix A to
C   the nearest complex matrix with an eigenvalue on the imaginary axis.
C   The number beta(A) is the minimum of the smallest singular value of
C   the matrix (A - jwI), where I is the identity matrix and j**2 = -1,
C   and the minimum is taken over all real w. The estimate is given as
C
C            low <= beta(A) <= high,
C
C   where either
C
C            (1 + tol) * low >= high,
C   or
C            low = 0   and   high = delta,
C
C   and delta is a small number approximately equal to the square root
C   of machine precision times the Frobenius norm of A. If all
C   eigenvalues of A lie in the open left half complex plane, then
C   beta(A) is the distance to the nearest unstable complex matrix,
C   i.e., the complex stability radius.
C
C Input parameters:
C   job    - option parameter indicating the task to be performed:
C            = 1 :  Hankel-norm of the alpha-stable projection of the
C                   transfer-function matrix G (D assumed 0);
C            = 2 :  H2 or L2 norm of a system;
C            = 3 :  complex stability radius, using bisection;
C            = 4 :  complex stability radius, using bisection and SVD.
C   A      - the n-by-n system state matrix A.
C   B      - the n-by-m system input matrix B.
C   C      - the p-by-n system output matrix C.
C   D      - the p-by-m system input/output matrix D.
C   dico   - (optional) specifies the type of the system:
C            = 1 :  continuous-time system;
C            = 2 :  discrete-time system.
C            Default: dico = 1.
C   equil  - (optional) specifies whether the system (A,B,C) should be
C            preliminarily equilibrated:
C            = 1 :  do not perform equilibration;
C            = 2 :  perform equilibration (scaling).
C            Default: equil = 1.
C   alpha  - (optional) constant specifying the alpha-stability boundary
C            for the eigenvalues of the state dynamics matrix A. For a
C            continuous-time system (dico = 1), alpha <= 0 is the
C            boundary value for the real parts of eigenvalues, while for
C            a discrete-time system (dico = 2), 0 <= alpha <= 1 is the
C            boundary value for the moduli of eigenvalues. The alpha-
C            stability domain is defined either as the open half complex
C            plane left to alpha, if dico = 1, or the interior of the
C            alpha-radius circle centered in the origin, if dico = 2.
C            Default: alpha =  -sqrt(epsilon_machine), if dico = 1;
C                     alpha = 1-sqrt(epsilon_machine), if dico = 2,
C            where epsilon_machine is the relative machine precision.
C   jobn   - (optional) constant specifying the norm to be computed:
C            = 1 :  the H2-norm;
C            = 2 :  the L2-norm.
C            Default: jobn = 1.
C   tol    - (optional) if job = 2, absolute tolerance level below
C            which the elements of B are considered zero (used for
C            controllability tests). If tol <= 0, default value is used.
C            Default: n*norm(B,1)*epsilon_machine.
C            If job = 3, tol specifies the accuracy with which low and
C            high approximate beta(A). If tol is less than
C            sqrt(epsilon_machine), then that value is used instead.
C            The recommended value is tol = 9, which gives an estimate
C            of beta(A) correct to within an order of magnitude.
C            Default: sqrt(epsilon_machine).
C            If job = 4, the accuracy with which beta(A) is to be
C            calculated. (See Comment 3 below.) If tol is set less than
C            epsilon_machine, then that value is used instead.
C            Default: epsilon_machine.
C
C Output parameters:
C   Hknorm - the Hankel-norm of the alpha-stable projection of G.
C   Hnorm  - the H2-norm of G, if jobn = 1, or the L2-norm of G,
C            if jobn = 2.
C   ns     - the order of the alpha-stable subsystem.
C   nq     - the order of the resulting numerator Q of the right coprime
C            factorization with inner denominator of G (see Comment 1).
C            Generally, nq = n - nus, where nus is the number of
C            uncontrollable unstable eigenvalues.
C   Ao     - if job = 1, the n-by-n resulted state dynamics matrix in
C            a block diagonal real Schur form with its eigenvalues
C            reordered and separated. Ao has two diagonal blocks: the
C            leading ns-by-ns part has eigenvalues in the alpha-
C            stability domain and the trailing (n-ns)-by-(n-ns) part
C            has eigenvalues outside the alpha-stability domain.
C            If job = 2, the nq-by-nq state dynamics matrix (in a real
C            Schur form) of the numerator factor Q of the right coprime
C            factorization with inner denominator of G (see Comment 1).
C   Bo     - if job = 1, the n-by-m input matrix of the transformed
C            system.
C            If job = 2, the nq-by-m input matrix of the numerator
C            factor Q of the right coprime factorization with inner
C            denominator of G.
C   Co     - if job = 1, the p-by-n output matrix of the transformed
C            system.
C            If job = 2, the p-by-nq output matrix of the numerator
C            factor Q of the right coprime factorization with inner
C            denominator of G.
C   Do     - if job = 2, the p-by-m input/output matrix of the numerator
C            factor Q of the right coprime factorization with inner
C            denominator of G.
C   hsv    - (optional) the ns Hankel singular values (ordered
C            decreasingly) of the alpha-stable part of the system.
C            hsv(1) is the Hankel norm of the alpha-stable subsystem.
C   low    - a lower bound for beta(A).
C   high   - (optional) an upper bound for beta(A).
C   beta   - the computed value of beta(A), which actually is an upper
C            bound.
C   omega  - (optional) the value of w such that the smallest singular
C            value of (A - jwI) equals beta(A).
C
C Comments
C   1) For job = 2, if the transfer-function matrix G is unstable,
C      then a right coprime factorization with inner denominator of G
C      is first computed
C               -1
C        G = Q*R  ,
C
C      where Q and R are stable transfer-function matrices and R is
C      inner. If G is stable, then Q = G and R = I.
C      Let (AQ,BQ,CQ,DQ) be the state-space representation of Q.
C
C      If dico = 1, then the L2-norm of G is computed as
C
C         norm2(G) = norm2(Q) = sqrt(trace(BQ'*X*BQ)),
C
C      where X satisfies the continuous-time Lyapunov equation
C
C         AQ'*X + X*AQ + CQ'*CQ = 0.
C
C      If dico = 2, then the l2-norm of G is computed as
C
C         norm2(G) = norm2(Q) = sqrt(trace(BQ'*X*BQ+DQ'*DQ)),
C
C      where X satisfies the discrete-time Lyapunov equation
C
C         AQ'*X*AQ - X + CQ'*CQ = 0.
C
C   2) For job = 3, the algorithm computes a lower bound low and an
C      upper bound high for beta(A) by a bisection method in the
C      following way. Given a non-negative real number sigma, the
C      Hamiltonian matrix H(sigma) is constructed:
C
C                     |   A      -sigma*I |     | A   G  |
C         H(sigma) =  |                   | :=  |        | .
C                     | sigma*I    -A'    |     | F  -A' |
C
C      It can be shown that H(sigma) has an eigenvalue whose real
C      part is zero if and only if sigma >= beta. Any lower and upper
C      bounds on beta(A) can be improved by choosing a number between
C      them and checking to see if H(sigma) has an eigenvalue with zero
C      real part.  This decision is made by computing the eigenvalues of
C      H(sigma) using the square reduced algorithm of Van Loan.
C
C   3) For job = 4, a reliable, quadratically convergent algorithm is
C      used, which finds (using a bisection strategy) an interval which
C      contains beta(A), and then switches to a modified bisection
C      strategy which converges quadratically to a minimizer. Note that
C      the efficiency of the strategy degrades if there are several
C      local minima that are near or equal the global minimum.
C      The computed function value beta satisfies
C
C         beta(A) <= beta + epsilon,
C
C         beta/(1+tol) - delta <= max(beta(A), sqrt(2*n*eps)*norm(A)),
C
C      where norm(A) is the Frobenius norm of A,
C
C         epsilon = p(n) * eps * norm(A),
C      and
C         delta   = p(n) * sqrt(eps) * norm(A),
C
C      eps = epsilon_machine, and p(n) is a low degree polynomial. It
C      is recommended to choose tol greater than sqrt(eps). Although
C      rounding errors can cause failure for smaller values of tol,
C      nevertheless, the calculation usually succeeds. Regardless of
C      success or failure, the first inequality holds.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2002.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009,
C   June 2012, July 2012, April 2017.
C
C **********************************************************************
C
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..      Print warnings is disabled (PRTWAR = .FALSE.).
      LOGICAL           PRTWAR
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( PRTWAR = .FALSE., ONE = 1.0D0, ZERO = 0.0D0 )
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
      CHARACTER         DICO, EQUIL, JOBN
      INTEGER           INFO, IWARN, LCWORK, LDA, LDB, LDC, LDD, LDWORK,
     $                  M, N, NQ, NS, P
      DOUBLE PRECISION  ALPHA, BETA, HIGH, LOW, OMEGA, TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      mwPointer         A, B, C, CWORK, D, DWORK, HSV
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           J1, J12, J2, J3
      CHARACTER*120     TEXT
      INTEGER           IP, ISYS, ITMP, JOB
      DOUBLE PRECISION  HNRM, TEMP
C
C .. External functions ..
      LOGICAL           MA02HD
      DOUBLE PRECISION  AB13AD, AB13BD, DLAMCH
      EXTERNAL          AB13AD, AB13BD, DLAMCH, MA02HD
C
C .. External subroutines ..
      EXTERNAL          AB13ED, AB13FD, DLACPY, DLASET
C
C ..Intrinsic functions..
      INTRINSIC         MAX, MIN, SQRT
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HNORM requires at least 2 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HNORM requires at least 1 output argument' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   job, A(nxn), (B(nxm), C(pxn), (D(pxm),) (dico,equil/jobn,alpha,)
C   (tol)
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR.
     $     mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      JOB = TEMP
      IF ( JOB.LT.1 .OR. JOB.GT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'JOB has 1, 2, 3, or 4 the only admissible values' )
      END IF
C
      J1  = JOB.EQ.1
      J2  = JOB.EQ.2
      J12 = J1 .OR. J2
      J3  = JOB.EQ.3
C
      IF ( J12 .AND. NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HNORM requires at least 4 input arguments' )
         IF ( J2 .AND. NRHS.LT.5 )
     $      CALL mexErrMsgTxt
     $           ( 'HNORM requires at least 5 input arguments' )
      END IF
C
      N = mxGetM( PRHS(2) )
C
      IF ( mxGetN( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a numeric matrix' )
      END IF
C
      IF ( J12 ) THEN
         M = mxGetN( PRHS(3) )
         P = mxGetM( PRHS(4) )
         IF ( mxGetM( PRHS(3) ).NE.N ) THEN
            CALL mexErrMsgTxt( 'B must have the same row dimension as A'
     $                       )
         END IF
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(3) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'B must be a numeric matrix' )
         END IF
         IF ( mxGetN( PRHS(4) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $          ('C must have the same column dimension as A' )
         END IF
         IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(4) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'C must be a numeric matrix' )
         END IF
         IP = 5
         IF ( J2 ) THEN
            IF ( mxGetM( PRHS(5) ).NE.P ) THEN
               CALL mexErrMsgTxt
     $             ( 'D must have the same row dimension as C' )
            END IF
            IF ( mxGetN( PRHS(5) ).NE.M ) THEN
               CALL mexErrMsgTxt
     $             ( 'D must have the same column dimension as B' )
            END IF
            IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(5) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'D must be a numeric matrix' )
            END IF
            IP = IP + 1
         END IF
C
C   dico
C
         IF ( NRHS.GE.IP ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
                 CALL mexErrMsgTxt( 'DICO must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
                 CALL mexErrMsgTxt( 'DICO must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            ISYS = TEMP
            IF ( ISYS.NE.1 .AND. ISYS.NE.2 ) THEN
                 CALL mexErrMsgTxt( 'DICO must be 1 or 2' )
            END IF
            IF ( ISYS.EQ.1 ) THEN
               DICO = 'C'
            ELSE
               DICO = 'D'
            END IF
            IP = IP + 1
         ELSE
            ISYS = 1
            DICO = 'C'
         END IF
C
C   equil/jobn
C
         IF ( NRHS.GE.IP ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               IF ( J1 ) THEN
                  CALL mexErrMsgTxt( 'EQUIL must be a scalar' )
               ELSE
                  CALL mexErrMsgTxt( 'JOBN must be a scalar' )
               END IF
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               IF ( J1 ) THEN
                  CALL mexErrMsgTxt( 'EQUIL must be an integer scalar' )
               ELSE
                  CALL mexErrMsgTxt( 'JOBN must be an integer scalar' )
               END IF
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            ITMP = TEMP
            IF ( ITMP.NE.1 .AND. ITMP.NE.2 ) THEN
               IF ( J1 ) THEN
                  CALL mexErrMsgTxt( 'EQUIL must be 1 or 2' )
               ELSE
                  CALL mexErrMsgTxt( 'JOBN must be 1 or 2' )
               END IF
            ELSE
               IF ( J1 ) THEN
                  IF ( ITMP.EQ.1 ) THEN
                     EQUIL = 'N'
                  ELSE
                     EQUIL = 'S'
                  END IF
               ELSE
                  IF ( ITMP.EQ.1 ) THEN
                     JOBN = 'H'
                  ELSE
                     JOBN = 'L'
                  END IF
               END IF
            END IF
            IP = IP + 1
         ELSE
            ITMP = 1
            IF ( J1 ) THEN
               EQUIL = 'N'
            ELSE
               JOBN = 'H'
            END IF
         END IF
C
C   alpha
C
         IF ( J1 .AND. NRHS.GE.IP ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
                 CALL mexErrMsgTxt( 'ALPHA must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
                 CALL mexErrMsgTxt( 'ALPHA must be a real scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), ALPHA, 1 )
         ELSE
            ALPHA = -SQRT( DLAMCH( 'Epsilon' ) )
            IF ( ISYS.EQ.2 )
     $         ALPHA = ONE + ALPHA
         END IF
      ELSE
         IP = 3
      END IF
C
C   tol
C
      IF ( .NOT.J1 .AND. NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'TOL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'TOL must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TOL, 1 )
      ELSE
         TOL = ZERO
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms.
C
      LDA = MAX( 1, N )
      IF ( J12 ) THEN
         LDB = LDA
         LDC = MAX( 1, P )
         IF ( J2 )
     $      LDD = LDC
      END IF
      IF ( J1 ) THEN
         LDWORK = N*( MAX( N, M, P ) + 5 ) + N*( N + 1 )/2
      ELSE IF ( J12 ) THEN
         LDWORK = MAX( M*( N + M ) + MAX( N*( N + 5 ), M*( M + 2 ),
     $                 4*P ), N*( MAX( N, P ) + 4 ) + MIN( N, P ) )
      ELSE IF ( J3 ) THEN
         LDWORK = 3*N*( N + 1 )
      ELSE
         LDWORK = 3*N*( N + 2 )
         LCWORK = MAX( 1, N*( N + 3 ) )
      END IF
      LDWORK = MAX( 1, LDWORK )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      A     = mxCalloc(  LDA*N, 8 )
      DWORK = mxCalloc( LDWORK, 8 )
      IF ( J12 ) THEN
         B = mxCalloc( LDB*M, 8 )
         C = mxCalloc( LDC*N, 8 )
         IF ( J1 ) THEN
            HSV = mxCalloc( N, 8 )
         ELSE
            D = mxCalloc( LDD*M, 8 )
         END IF
      ELSE IF ( JOB.EQ.4 ) THEN
         CWORK = mxCalloc( LCWORK, 16 )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL(A), N*N )
      IF ( J12 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL(B), N*M )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL(C), P*N )
         IF ( J2 ) THEN
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), %VAL(D),
     $                             P*M )
C
C Test whether D is the zero matrix.
C
            IF ( ISYS.EQ.1 .AND. MIN( P, M ).GT.0 .AND.
     $           .NOT.MA02HD( 'Full', P, M, ZERO, %VAL(D), LDD ) ) THEN
               WRITE( TEXT, '('' Matrix D is taken as zero. '')' )
               CALL mexPrintf( TEXT )
               CALL DLASET( 'Full', P, M, ZERO, ZERO, %VAL(D), LDD )
            END IF
         END IF
      END IF

C
C Do the actual computations.
C
      IF ( J1 ) THEN
         HNRM = AB13AD( DICO, EQUIL, N, M, P, ALPHA, %VAL(A), LDA,
     $                  %VAL(B), LDB, %VAL(C), LDC, NS, %VAL(HSV),
     $                  %VAL(DWORK), LDWORK, INFO )
      ELSE IF ( J12 ) THEN
         HNRM = AB13BD( DICO, JOBN, N, M, P, %VAL(A), LDA, %VAL(B), LDB,
     $                  %VAL(C), LDC, %VAL(D), LDD, NQ, TOL,
     $                  %VAL(DWORK), LDWORK, IWARN, INFO )
      ELSE IF ( J3 ) THEN
         CALL AB13ED( N, %VAL(A), LDA, LOW, HIGH, TOL, %VAL(DWORK),
     $                LDWORK, INFO )
      ELSE
         CALL AB13FD( N, %VAL(A), LDA, BETA, OMEGA, TOL, %VAL(DWORK),
     $                LDWORK, %VAL(CWORK), LCWORK, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( 1, 1, 0 )
      IF ( J12 ) THEN
         CALL mxCopyReal8ToPtr( HNRM, mxGetPr( PLHS(1) ), 1 )
      ELSE IF ( J3 ) THEN
         CALL mxCopyReal8ToPtr( LOW, mxGetPr( PLHS(1) ), 1 )
      ELSE
         CALL mxCopyReal8ToPtr( BETA, mxGetPr( PLHS(1) ), 1 )
      END IF
C
      IF( NLHS.GT.1 ) THEN
         PLHS(2) = mxCreateDoubleMatrix( 1, 1, 0 )
         IF ( J12 ) THEN
            IF ( J2 )
     $         NS = NQ
            TEMP = NS
            CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS(2) ), 1 )
         ELSE IF ( J3 ) THEN
            CALL mxCopyReal8ToPtr( HIGH, mxGetPr( PLHS(2) ), 1 )
         ELSE
            CALL mxCopyReal8ToPtr( OMEGA, mxGetPr( PLHS(2) ), 1 )
         END IF
C
         IF( J12 .AND. NLHS.GT.2 ) THEN
            IF ( NS.GT.0 .AND. NS.NE.N ) THEN
               CALL DLACPY( 'Full', NS, NS, %VAL(A), LDA, %VAL(A), NS )
               CALL DLACPY( 'Full', NS, M,  %VAL(B), LDB, %VAL(B), NS )
            END IF
            PLHS(3) = mxCreateDoubleMatrix( NS, NS, 0 )
            CALL mxCopyReal8ToPtr( %VAL(A), mxGetPr( PLHS(3) ),
     $                             NS*NS )
            PLHS(4) = mxCreateDoubleMatrix( NS, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL(B), mxGetPr( PLHS(4) ),
     $                             NS*M )
            PLHS(5) = mxCreateDoubleMatrix( P, NS, 0 )
            CALL mxCopyReal8ToPtr( %VAL(C), mxGetPr( PLHS(5) ),
     $                             P*NS )
            IF( NLHS.GT.5 ) THEN
               IF( J1 ) THEN
                  PLHS(6) = mxCreateDoubleMatrix( NS, MIN( NS, 1 ),
     $                                            0 )
                  CALL mxCopyReal8ToPtr( %VAL(HSV),
     $                                   mxGetPr( PLHS(6) ),
     $                                   NS*MIN( NS, 1 ) )
               ELSE
                  PLHS(6) = mxCreateDoubleMatrix( P, M, 0 )
                  CALL mxCopyReal8ToPtr( %VAL(D),
     $                                   mxGetPr( PLHS(6) ), P*M )
               END IF
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      CALL mxFree( A )
      CALL mxFree( DWORK )
      IF ( J12 ) THEN
         CALL mxFree( B )
         CALL mxFree( C )
         IF ( J1 ) THEN
            CALL mxFree( HSV )
         ELSE
            CALL mxFree( D )
         END IF
      ELSE IF ( JOB.EQ.4 ) THEN
         CALL mxFree( CWORK )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( J1 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB13AD'')'
     $           ) INFO
         ELSE IF ( J2 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB13BD'')'
     $           ) INFO
         ELSE IF ( J3 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB13ED'')'
     $           ) INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB13FD'')'
     $           ) INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
C
      ELSE IF ( IWARN.NE.0 .AND. J2 .AND. PRTWAR ) THEN
         WRITE( TEXT, '('' IWARN = '',I4,'' ON EXIT FROM AB13BD'')'
     $        ) IWARN
         CALL mexPrintf( TEXT )
      END IF
C
      RETURN
C *** Last line of HNORM ***
      END
