#include "fintrf.h"
C TOTALLS.F - Gateway function for solving the Total Least Squares
C             (TLS) problem using a singular value decomposition (SVD)
C             approach or a Partial SVD (PSVD) approach, based on
C             SLICOT routines MB02MD and MB02ND.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [X(,V(,S)(,Q,NL),ro(,to),rcnd)] = TotalLS(A,B(,meth(,job),r(,t),
C                                                 tol(,reltol),printw))
C   [X(,V,S,ro,rcnd)]               = TotalLS(A,B,1(,job,r,tol,printw))
C   [X(,V,Q,NL,ro,to,rcnd)]         = TotalLS(A,B(,2,r,t,tol,reltol,
C                                                 printw))
C
C Purpose:
C   To solve the Total Least Squares (TLS) problem using a Singular
C   Value Decomposition (SVD) approach, or a Partial SVD (PSVD)
C   approach. The TLS problem assumes an overdetermined set of linear
C   equations AX = B, where both the data matrix A as well as the
C   observation matrix B are inaccurate. The gateway also solves
C   determined and underdetermined sets of equations by computing the
C   minimum norm solution. It is assumed that all preprocessing measures
C   (scaling, coordinate transformations, whitening, ... ) of the data
C   have been performed in advance.
C
C Input parameters:
C   A      - the m-by-n data matrix A.
C   B      - the m-by-l observation matrix B (right-hand sides).
C            Define the matrix C = [A|B] by concatenating A and B.
C   meth   - (optional) scalar indicating the method to be used for
C            solving the TLS problem, as follows:
C            = 1 :  SVD  method;
C            = 2 :  PSVD method.
C            Default: meth = 2.
C   job    - (optional) if meth = 1, scalar determining whether the
C            values of r and tol are given as input parameters or are
C            computed by the gateway, as follows:
C            = 0 :  Compute neither r nor tol;
C            = 1 :  Compute r only;
C            = 2 :  Compute tol only;
C            = 3 :  Compute both r and tol.
C            Default: job = 3.
C   r      - (optional) if meth = 1, and job = 0 or job = 2, or if
C            meth = 2 and r >= 0, scalar specifying the rank of the
C            TLS approximation [A+DA|B+DB];  r <= min(m,n) and r >= 0
C            for meth = 1.
C            Otherwise, the rank is computed by the gateway.
C            If meth = 1, and job = 1 or job = 3, then r must not be
C            specified as an input parameter.
C            Default for meth = 2: r = -1.
C   t      - (optional) if meth = 2, and r < 0, then the rank of the
C            TLS approximation [A+DA|B+DB] is computed using the real
C            number t as (min(m,n+l) - d), where d is the number of
C            singular values of [A|B] less than or equal to t;  t >= 0.
C            If meth = 2 and r >= 0, t is an initial estimate for
C            computing a lower bound on the r largest singular values
C            of [A|B]. If t < 0 on entry however, then t is computed
C            by the gateway.
C            Default: t = EPS, if r <  0, where EPS is the machine
C                              precision (see LAPACK Library routine
C                              DLAMCH);
C                     t = -1,  if r >= 0.
C            If meth = 1, t must not be specified as an input parameter.
C   tol    - (optional) if meth = 1, real scalar defining the tolerance
C            used to determine the rank of the TLS approximation
C            [A+DA|B+DB] and to check the multiplicity of the singular
C            values of matrix C. Specifically, S(i) and S(j) (i < j) are
C            considered to be equal if SQRT(S(i)**2 - S(j)**2) <= tol,
C            and the TLS approximation [A+DA|B+DB] has rank r if
C            S(i) > tol*S(1) (or S(i) > tol, if tol specifies sdev
C            (see below)), for i = 1,2,...,r. The value tol is also used
C            to check the singularity of the upper triangular matrix F
C            (as defined in Method).
C            If job = 0 or job = 1, then tol must specify the desired
C            tolerance. If tol <= 0, the tolerance is taken as EPS.
C            Otherwise, the tolerance is computed by the gateway and the
C            input value of tol defines a non-negative value sdev, i.e.,
C            the estimated standard deviation of the error on each
C            element of the matrix C.
C            If meth = 2, real scalar defining the multiplicity of
C            singular values by considering all singular values within
C            an interval of length tol as coinciding. The value tol is
C            used in checking how many singular values are less than or
C            equal to t. Also in computing an appropriate upper bound t
C            by a bisection method, tol is used as a stopping criterion
C            defining the minimum (absolute) subinterval width. The
C            value tol is also taken as an absolute tolerance for
C            negligible elements in the QR/QL iterations. If tol <= 0,
C            then the tolerance is taken as specified in SLICOT Library
C            routine MB04YD.
C            Default: tol = 0.
C   reltol - (optional) if meth = 2, real scalar specifying the minimum
C            relative width of an interval. When an interval is narrower
C            than tol, or than reltol times the larger (in magnitude)
C            endpoint, then it is considered to be sufficiently small
C            and bisection has converged. If reltol <= BASE * EPS, where
C            BASE is machine radix, then the tolerance is taken
C            as BASE * EPS.
C            Default: reltol = 0.
C            If meth = 1, reltol must not be specified as an input
C            parameter.
C   printw - (optional) switch for printing the warning messages.
C            = 1:  print warning messages;
C            = 0:  do not print warning messages.
C            Default: printw = 0.
C
C Output parameters:
C   X      - the n-by-l solution matrix to the TLS problem specified
C            by A and B.
C   V      - (optional) if meth = 1, the (n+l)-by-(n+l) matrix of
C            (transformed) right singular vectors, including null space
C            vectors, if any, of C. Specifically, the first ro columns
C            (see the description of ro below) always contain the first
C            ro right singular vectors, corresponding to the largest
C            singular values of C. If l = 0 or ro = 0, the remaining
C            n+l-ro columns contain the remaining right singular
C            vectors. Otherwise, these n+l-ro columns contain the matrix
C            V2 transformed as described in Step 3 of the TLS algorithm
C            (see Method).
C            If meth = 2, the columns of V whose index i corresponds
C            with NL(i) = 1, are the possibly transformed n+l-ro base
C            vectors of the right singular subspace corresponding to the
C            singular values of C = [A|B] which are less than or equal
C            to to (see below the description of the parameters ro and
C            to). Specifically, if l = 0 or ro = 0, these vectors are
C            indeed the base vectors above. Otherwise, these vectors
C            form the matrix V2, transformed as described in Step 4 of
C            the PTLS algorithm (see Method). The TLS solution is
C            computed from these vectors. The other columns contain no
C            useful information.
C   S      - (optional) if meth = 1, the min(m,n+l) singular values of
C            the matrix C, ordered such that S(1) >= S(2) >= ... >=
C            >= S(p-1) >= S(p) >= 0, where p = min(m,n+l).
C            If meth = 2, then S must not be specified as an output
C            parameter.
C   Q      - (optional) if meth = 2, the max(1,2*p-1) vector containing
C            the partially diagonalized bidiagonal matrix J computed
C            from C, at the moment that the desired singular subspace
C            has been found. Specifically, the leading p entries of Q
C            contain the diagonal elements q(1),q(2),...,q(p) and the
C            entries Q(p+1), ..., Q(2*p-1) contain the superdiagonal
C            elements e(1),...,e(p-1) of J.
C            If meth = 1, then Q must not be specified as an output
C            parameter.
C   NL     - (optional) if meth = 2, an (n+l) vector; the indices of
C            the elements of this vector with value 1 indicate the
C            columns in V containing the (transformed) base vectors of
C            the right singular subspace of C from which the TLS
C            solution has been computed.
C            If meth = 1, then NL must not be specified as an output
C            parameter.
C   ro     - (optional) if meth = 1, and job = 1 or job = 3, or if
C            meth = 2, and r < 0, then ro contains the computed
C            (effective) rank of the TLS approximation [A+DA|B+DB].
C            If meth = 1, and job = 0 or job = 2, or if meth = 2, and
C            r >= 0, the input value r may be changed by the gateway if
C            the r-th and the (r+1)-th singular values of C = [A|B] are
C            considered to be equal, or if the upper triangular matrix F
C            (as defined in Method) is (numerically) singular.
C   to     - (optional) if meth = 2 and r >= 0 on entry, then to
C            contains the computed bound such that precisely ro singular
C            values of C = [A|B] are greater than to + tol.
C            If meth = 2, and r < 0, then to = t.
C            If meth = 1, then to must not be specified as an output
C            parameter.
C   rcnd   - (optional) real scalar containing the reciprocal of the
C            condition number of the matrix F.
C
C Method:
C   The methods used are extensions of the classical TLS algorithm.
C
C   Let [A|B] denote the matrix formed by adjoining the columns of B to
C   the columns of A on the right.
C
C   Total Least Squares (TLS) definition:
C
C     Given matrices A and B, find a matrix X satisfying
C
C          (A + DA) X = B + DB,
C
C     where A and DA are m-by-n matrices, B and DB are m-by-l matrices
C     and X is an n-by-l matrix.
C     The solution X must be such that the Frobenius norm of [DA|DB] is
C     a minimum and each column of B + DB is in the range of A + DA.
C     Whenever the solution is not unique, the minimum norm solution X
C     is singled out.
C
C   Define matrix C = [A|B] and s(i) as its i-th singular value for
C   i = 1,...,min(m,n+l). If m < n+l, then s(j) = 0 for j = m+1,...,n+l.
C
C   The Classical TLS algorithm proceeds as follows:
C
C   Step 1: Compute part of the singular value decomposition (SVD)
C           USV' of C = [A|B], namely compute S and V'. (An initial
C           QR factorization of C is used when m is larger enough
C           than n+l.)
C
C   Step 2: If the rank r is not given, compute the rank r0 of the data
C           [A|B] based on tol as follows: if job = 0 or job = 1,
C
C              s(1) >= ... >= s(r0) > tol*s(1) >= ... >= s(n+l).
C
C           Otherwise, tol can be computed from the standard deviation
C           sdev of the errors on [A|B]:
C
C              tol = sqrt(2 * max(m,n+l)) * sdev,
C
C           and the rank r0 is determined (if job = 1 or 3) using
C
C              s(1) >= ... >= s(r0) > tol >= ... >= s(n+l).
C
C           The rank r of the approximation [A+DA|B+DB] is then equal
C           to the minimum of n and r0.
C
C   Step 3: Let V2 be the matrix of the columns of V corresponding to
C           the (n+l-r) smallest singular values of C, i.e., the last
C           (n+l-r) columns of V. Using Householder transformations,
C           compute the orthogonal matrix Q such that:
C
C                     |VH   Y|
C            V2 x Q = |      |,                                      (1)
C                     |0    F|
C
C           where VH is an n-by-(n-r) matrix, Y is an n-by-l matrix and
C           F is an l-by-l upper triangular matrix. If F is singular,
C           then lower the rank r with the multiplicity of s(r) and
C           repeat this step.
C
C   Step 4: If F is nonsingular then the solution X is obtained by
C           solving the following equations by forward elimination:
C
C              X F = -Y.                                             (2)
C
C   Consider now the Partial Total Least Squares (PTLS) approach. Since
C   the TLS solution can be computed from any orthogonal basis of the
C   right singular subspace corresponding to the smallest singular
C   values of C, the Partial Singular Value Decomposition (PSVD) can be
C   used instead of the classical SVD. The dimension of this subspace
C   may be determined by the rank of C or by an upper bound for those
C   smallest singular values.
C
C   The PTLS algorithm proceeds as follows:
C
C   Step 1: Bidiagonalization phase
C    (a) If m is large enough than n + l, transform C into upper
C        triangular form R by Householder transformations.
C    (b) Transform C (or R) into upper bidiagonal form (p = min(m,n+l)):
C
C                   |q(1) e(1)  0   ...  0   |
C              (0)  | 0   q(2) e(2)      .   |
C             J   = | .                  .   |
C                   | .                e(p-1)|
C                   | 0             ... q(p) |
C
C        if m >= n + l, or lower bidiagonal form:
C
C                   |q(1)  0    0   ...  0     0   |
C              (0)  |e(1) q(2)  0        .     .   |
C             J   = | .                  .     .   |
C                   | .                 q(p)   .   |
C                   | 0             ... e(p-1) q(p)|
C
C        if m < n + l, using Householder transformations. In the second
C        case, transform the matrix to the upper bidiagonal form by
C        applying Givens rotations.
C    (c) Initialize the right singular base matrix with the identity
C        matrix.
C
C   Step 2: Partial diagonalization phase
C   If the upper bound t is not given, then compute t such that
C   precisely p - r singular values (p=min(m,n+l)) of the bidiagonal
C   matrix are less than or equal to t, using a bisection method.
C   Diagonalize the given bidiagonal matrix J partially, using either QL
C   iterations (if the upper left diagonal element of the considered
C   bidiagonal submatrix is smaller than the lower right diagonal
C   element) or QR iterations, such that J is split into unreduced
C   bidiagonal submatrices whose singular values are either all larger
C   than t or are all less than or equal to t. Accumulate the Givens
C   rotations in V.
C
C   Step 3: Back transformation phase
C   Apply the Householder transformations of Step 1(b) onto the base
C   vectors of V associated with the bidiagonal submatrices with all
C   singular values less than or equal to t.
C
C   Step 4: Computation of F and Y
C   Compute the factorization in (1). If F is singular, then reduce the
C   value of r by one and repeat Steps 2, 3 and 4.
C
C   Step 5: Computation of the TLS solution
C   If F is non-singular then the solution X is obtained by solving
C   the equations (2) by forward elimination.
C
C Comments
C   1. The TLS solution is unique if r = n, F is nonsingular and
C   s(n) > s(n+1).
C   If F is singular, however, then the computed solution is infinite
C   and hence does not satisfy the second TLS criterion (see TLS
C   definition). These problems are called nongeneric and the TLS
C   computations have been generalized in order to solve them. The
C   generalization satisfies the TLS criteria for any number l of
C   observation vectors in B provided that, in addition, the solution
C   | X| is constrained to be orthogonal to all vectors of the form |w|
C   |-I|                                                            |0|
C   which belong to the space generated by the columns of the
C   submatrix |Y|.
C             |F|
C
C   2. If r is lowered in Step 4 of the PTLS algorithm, some additional
C   base vectors must be computed in Step 2. The additional computations
C   are kept to a minimum. If r is lowered in Step 4 but the
C   multiplicity of the r-th singular value is larger than 1, then the
C   value of r is further lowered with its multiplicity defined by the
C   parameter tol. This is done at the beginning of Step 2, using a
C   bisection method to estimate t.
C
C   3. The computational efficiency of the PTLS algorithm compared with
C   the classical TLS algorithm is obtained by making use of PSVD
C   instead of performing the entire SVD. Depending on the gap between
C   the r-th and the (r+1)-th singular values of C, the number (n+l-r)
C   of base vectors to be computed with respect to the column dimension
C   (n+l) of C and the desired accuracy reltol, the algorithm used by
C   the PTLS algorithm is approximately twice as fast as the classical
C   TLS algorithm at the expense of extra storage requirements, namely:
C     (n + l) x (n + l - 1)/2,  if m >= n + l or
C     m x (n + l - (m - 1)/2),  if m <  n + l.
C   This is because the Householder transformations performed on the
C   rows of C in the bidiagonalization phase (see Step 1) must be kept
C   until the end (Step 5).
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Sep. 2003.
C
C Revisions:
C   V. Sima, Feb. 2004, Apr. 2009, July 2011, Aug. 2012.
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
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         JOB
      INTEGER           INFO, IWARN, L, LDC, LDWORK, LDX, M, N, RANK
      DOUBLE PRECISION  RELTOL, THETA, TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL, ALLOCATABLE ::          BWORK(:), INUL(:)
      INTEGER, ALLOCATABLE ::          IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: DW(:)
      mwPointer                        C, DWORK, Q, S, X
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           PRINTW
      CHARACTER*120     TEXT
      INTEGER           I, IJOB, IP, ITMP, J, LIWORK, LQ, LW, METH, MNL,
     $                  NL
      DOUBLE PRECISION  RCND(1), TEMP
C
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C
C .. External subroutines ..
      EXTERNAL          DLACPY, MA02KV, MB02MD, MB02ND
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN, MOD
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'TOTALLS requires at least 2 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'TOTALLS requires at least 1 output argument' )
      END IF
C
      M = mxGetM( PRHS(1) )
      N = mxGetN( PRHS(1) )
      L = mxGetN( PRHS(2) )
C
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF ( mxGetM( PRHS(2) ).NE.M ) THEN
         CALL mexErrMsgTxt( 'B must have the same number of rows as A' )
      END IF
C
      IP = 3
      IF ( NRHS.GE.IP ) THEN
C
C   meth
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'METH must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'METH must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         METH = TEMP
         IF ( METH.LT.1 .OR. METH.GT.2 ) THEN
            CALL mexErrMsgTxt
     $      ( 'METH has 1 or 2 the only admissible values')
         END IF
         IP = IP + 1
      ELSE
         METH = 2
      END IF
C
      IF ( METH.EQ.1 ) THEN
         IF ( NRHS.GE.IP ) THEN
C
C   job
C
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'JOB must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            IJOB = TEMP
            IF ( IJOB.LT.0 .OR. IJOB.GT.3 ) THEN
               CALL mexErrMsgTxt
     $           ( 'JOB has 0, 1, 2, or 3 the only admissible values' )
            END IF
            IP = IP + 1
         ELSE
            IJOB = 3
         END IF
C
         IF ( IJOB.EQ.0 ) THEN
            JOB = 'N'
         ELSEIF ( IJOB.EQ.1 ) THEN
            JOB = 'R'
         ELSEIF ( IJOB.EQ.2 ) THEN
            JOB = 'T'
         ELSE
            JOB = 'B'
         END IF
      ELSE
         IJOB = 0
      END IF
C
      IF ( NRHS.GE.IP .AND. MOD( IJOB, 2 ).EQ.0 ) THEN
C
C   r
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'R must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'R must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         RANK = TEMP
         IF ( RANK.GT.MIN( M, N ) ) THEN
            CALL mexErrMsgTxt( 'RANK cannot exceed min( size( A ) )' )
         ELSE IF ( METH.EQ.1 .AND. RANK.LT.0 ) THEN
            CALL mexErrMsgTxt( 'RANK must be a non-negative integer' )
         END IF
         IP = IP + 1
      ELSE
         RANK = -1
      END IF
C
      IF ( NRHS.GE.IP .AND. METH.EQ.2 ) THEN
C
C   t
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'T must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'T must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), THETA, 1 )
         IP = IP + 1
      ELSE IF ( RANK.LT.0 ) THEN
         THETA = DLAMCH( 'Epsilon' )
      ELSE
         THETA = -ONE
      END IF
C
      IF ( NRHS.GE.IP ) THEN
C
C   tol
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TOL, 1 )
         IP = IP + 1
      ELSE
         TOL = ZERO
      END IF
C
      IF ( NRHS.GE.IP .AND. METH.EQ.2 ) THEN
C
C   reltol
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'RELTOL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'RELTOL must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), RELTOL, 1 )
         IP = IP + 1
      ELSE
         RELTOL = ZERO
      END IF
C
C     printw
C
      PRINTW = .FALSE.
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'PRINTW must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'PRINTW must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ITMP = TEMP
         IF ( ITMP.LT.0 .OR. ITMP.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'PRINTW has 0 or 1 the only admissible values' )
         END IF
         PRINTW = ITMP.EQ.1
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms
C in MB02MD and MB02ND.
C
      NL  = N + L
      LDC = MAX( 1, M, NL )
      LDX = MAX( 1, N )
      MNL = MIN( M, NL )
      IF ( METH.EQ.1 ) THEN
         LIWORK = L
         IF ( M.GE.NL ) THEN
            LDWORK = MAX( 2, 3*NL + M, 5*NL )
         ELSE
            LDWORK = MAX( 2, M*NL + MAX( 3*M + NL, 5*M ), 3*L )
         END IF
      ELSE
         LIWORK = NL + L
         LQ = MAX( 1, 2*MNL - 1 )
         IF ( M.GE.NL ) THEN
            LW = ( NL*( NL - 1 ) )/2
         ELSE
            LW = M*NL - ( M*( M - 1 ) )/2
         END IF
         LDWORK = MAX( 2, MAX( M, NL ) + 2*MNL,
     $                    MNL + LW + MAX( 6*NL - 5,
     $                                    L*L + MAX( NL, 3*L ) ) )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( IWORK( LIWORK ) )
      C     = mxCalloc( LDC*NL, 8 )
      DWORK = mxCalloc( LDWORK, 8 )
      X     = mxCalloc(  LDX*L, 8 )
      IF ( METH.EQ.1 ) THEN
         S = mxCalloc( MNL, 8 )
      ELSE
         ALLOCATE ( BWORK( NL ), INUL( NL ) )
         Q = mxCalloc( LQ, 8 )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      ALLOCATE ( DW( M*NL ) )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), DW, M*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), DW(M*N+1), M*L )
      CALL DLACPY( 'Full', M, NL, DW, MAX( 1, M ), %VAL(C), LDC )
      DEALLOCATE ( DW )
C
C Do the actual computations.
C
      IF ( METH.EQ.1 ) THEN
         CALL MB02MD( JOB, M, N, L, RANK, %VAL(C), LDC, %VAL(S),
     $                %VAL(X), LDX, TOL, IWORK, %VAL(DWORK), LDWORK,
     $                IWARN, INFO )
      ELSE
         CALL MB02ND( M, N, L, RANK, THETA, %VAL(C), LDC, %VAL(X), LDX,
     $                %VAL(Q), INUL, TOL, RELTOL, IWORK, %VAL(DWORK),
     $                LDWORK, BWORK, IWARN, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         CALL MA02KV( 1, 2, 1, %VAL( DWORK ), RCND )
         IP   = 1
         PLHS(IP) = mxCreateDoubleMatrix( N, L, 0 )
         CALL mxCopyReal8ToPtr( %VAL(X), mxGetPr( PLHS(IP) ), N*L )
         IP = IP + 1
C
         IF ( NLHS.GE.IP ) THEN
            IF ( LDC.GT.NL )
     $         CALL DLACPY( 'Full', NL, NL, %VAL(C), LDC, %VAL(C),
     $                      MAX( 1, NL ) )
            PLHS(IP) = mxCreateDoubleMatrix( NL, NL, 0 )
            CALL mxCopyReal8ToPtr( %VAL(C), mxGetPr( PLHS(IP) ),
     $                             NL*NL )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            IF ( METH.EQ.1 ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( MNL, MIN( 1, MNL ),
     $                                          0 )
               CALL mxCopyReal8ToPtr( %VAL(S), mxGetPr( PLHS(IP) ),
     $                                MNL*MIN( 1, MNL ) )
               IP = IP + 1
            ELSE
               PLHS(IP) = mxCreateDoubleMatrix( MAX( 1, 2*MNL - 1 ),
     $                                          1, 0 )
               CALL mxCopyReal8ToPtr( %VAL(Q), mxGetPr( PLHS(IP) ),
     $                                MAX( 1, 2*MNL - 1 ) )
               IP = IP + 1
C
               IF ( NLHS.GE.IP ) THEN
C
                  ALLOCATE ( DW( NL ) )
                  DO 10 J = 1, NL
                     IF ( INUL(J) ) THEN
                        DW(J) = ONE
                     ELSE
                        DW(J) = ZERO
                     END IF
   10             CONTINUE
C
                  PLHS(IP) = mxCreateDoubleMatrix( NL, MIN( 1, NL ),
     $                                             0 )
                  CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(IP) ),
     $                                   NL*MIN( 1, NL ) )
                  DEALLOCATE ( DW )
                  IP = IP + 1
               END IF
            END IF
         END IF
C
         IF ( NLHS.GE.IP ) THEN
C
            TEMP = RANK
            PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS(IP) ), 1 )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. METH.EQ.2 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( THETA, mxGetPr( PLHS(IP) ), 1 )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( RCND, mxGetPr( PLHS(IP) ), 1 )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE ( IWORK )
      CALL mxFree( C     )
      CALL mxFree( DWORK )
      CALL mxFree( X     )
      IF ( METH.EQ.1 ) THEN
         CALL mxFree( S  )
      ELSE
         DEALLOCATE ( BWORK, INUL )
         CALL mxFree( Q )
      END IF
C
C Error and warning handling.
C
      IF ( IWARN.NE.0 .AND. PRINTW ) THEN
         IF ( METH.EQ.1 ) THEN
            WRITE( TEXT, '('' IWARN = '',I4,'' ON EXIT FROM MB02MD'')'
     $           ) IWARN
         ELSE
            WRITE( TEXT, '('' IWARN = '',I4,'' ON EXIT FROM MB02ND'')'
     $           ) IWARN
         END IF
         CALL mexPrintf( TEXT )
      END IF
C
      IF ( INFO.NE.0 ) THEN
         IF ( METH.EQ.1 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB02MD'')'
     $           ) INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB02ND'')'
     $           ) INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of TOTALLS ***
      END
