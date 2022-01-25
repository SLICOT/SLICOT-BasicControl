% TOTALLS.F - MEX-function function for solving the Total Least Squares
%             (TLS) problem using a singular value decomposition (SVD)
%             approach or a Partial SVD (PSVD) approach, based on
%             SLICOT routines MB02MD and MB02ND.
%
%   [X(,V(,S)(,Q,NL),ro(,to),rcnd)] = TotalLS(A,B(,meth(,job),r(,t),
%                                                 tol(,reltol),printw))
%   [X(,V,S,ro,rcnd)]               = TotalLS(A,B,1(,job,r,tol,printw))
%   [X(,V,Q,NL,ro,to,rcnd)]         = TotalLS(A,B(,2,r,t,tol,reltol,
%                                                 printw))
%
%   TotalLS solves the Total Least Squares (TLS) problem using a
%   Singular Value Decomposition (SVD) approach, or a Partial SVD (PSVD)
%   approach. The TLS problem assumes an overdetermined set of linear
%   equations AX = B, where both the data matrix A as well as the
%   observation matrix B are inaccurate. The gateway also solves
%   determined and underdetermined sets of equations by computing the
%   minimum norm solution. It is assumed that all preprocessing measures
%   (scaling, coordinate transformations, whitening, ... ) of the data
%   have been performed in advance.
%
%   Description of input parameters:
%   A      - the m-by-n data matrix A.
%   B      - the m-by-l observation matrix B (right-hand sides).
%            Define the matrix C = [A|B] by concatenating A and B.
%   meth   - (optional) scalar indicating the method to be used for
%            solving the TLS problem, as follows:
%            = 1 :  SVD  method;
%            = 2 :  PSVD method.
%            Default: meth = 2.
%   job    - (optional) if meth = 1, scalar determining whether the
%            values of r and tol are given as input parameters or are
%            computed by the gateway, as follows:
%            = 0 :  Compute neither r nor tol;
%            = 1 :  Compute r only;
%            = 2 :  Compute tol only;
%            = 3 :  Compute both r and tol.
%            Default: job = 3.
%   r      - (optional) if meth = 1, and job = 0 or job = 2, or if
%            meth = 2 and r >= 0, scalar specifying the rank of the
%            TLS approximation [A+DA|B+DB];  r <= min(m,n) and r >= 0
%            for meth = 1.
%            Otherwise, the rank is computed by the gateway.
%            If meth = 1, and job = 1 or job = 3, then r must not be
%            specified as an input parameter.
%            Default for meth = 2: r = -1.
%   t      - (optional) if meth = 2, and r < 0, then the rank of the
%            TLS approximation [A+DA|B+DB] is computed using the real
%            number t as (min(m,n+l) - d), where d is the number of
%            singular values of [A|B] less than or equal to t;  t >= 0.
%            If meth = 2 and r >= 0, t is an initial estimate for
%            computing a lower bound on the r largest singular values
%            of [A|B]. If t < 0 on entry however, then t is computed
%            by the gateway.
%            Default: t = EPS, if r <  0, where EPS is the machine
%                              precision (see LAPACK Library routine
%                              DLAMCH);
%                     t = -1,  if r >= 0.
%            If meth = 1, t must not be specified as an input parameter.
%   tol    - (optional) if meth = 1, real scalar defining the tolerance
%            used to determine the rank of the TLS approximation
%            [A+DA|B+DB] and to check the multiplicity of the singular
%            values of matrix C. Specifically, S(i) and S(j) (i < j) are
%            considered to be equal if SQRT(S(i)**2 - S(j)**2) <= tol,
%            and the TLS approximation [A+DA|B+DB] has rank r if
%            S(i) > tol*S(1) (or S(i) > tol, if tol specifies sdev
%            (see below)), for i = 1,2,...,r. The value tol is also used
%            to check the singularity of the upper triangular matrix F
%            (as defined in Method).
%            If job = 0 or job = 1, then tol must specify the desired
%            tolerance. If tol <= 0, the tolerance is taken as EPS.
%            Otherwise, the tolerance is computed by the gateway and the
%            input value of tol defines a non-negative value sdev, i.e.,
%            the estimated standard deviation of the error on each
%            element of the matrix C.
%            If meth = 2, real scalar defining the multiplicity of
%            singular values by considering all singular values within
%            an interval of length tol as coinciding. The value tol is
%            used in checking how many singular values are less than or
%            equal to t. Also in computing an appropriate upper bound t
%            by a bisection method, tol is used as a stopping criterion
%            defining the minimum (absolute) subinterval width. The
%            value tol is also taken as an absolute tolerance for
%            negligible elements in the QR/QL iterations. If tol <= 0,
%            then the tolerance is taken as specified in SLICOT Library
%            routine MB04YD.
%            Default: tol = 0.
%   reltol - (optional) if meth = 2, real scalar specifying the minimum
%            relative width of an interval. When an interval is narrower
%            than tol, or than reltol times the larger (in magnitude)
%            endpoint, then it is considered to be sufficiently small
%            and bisection has converged. If reltol <= BASE * EPS, where
%            BASE is machine radix, then the tolerance is taken
%            as BASE * EPS.
%            Default: reltol = 0.
%            If meth = 1, reltol must not be specified as an input
%            parameter.
%   printw - (optional) switch for printing the warning messages.
%            = 1:  print warning messages;
%            = 0:  do not print warning messages.
%            Default: printw = 0.
%
%   Description of output parameters:
%   X      - the n-by-l solution matrix to the TLS problem specified
%            by A and B.
%   V      - (optional) if meth = 1, the (n+l)-by-(n+l) matrix of
%            (transformed) right singular vectors, including null space
%            vectors, if any, of C. Specifically, the first ro columns
%            (see the description of ro below) always contain the first
%            ro right singular vectors, corresponding to the largest
%            singular values of C. If l = 0 or ro = 0, the remaining
%            n+l-ro columns contain the remaining right singular
%            vectors. Otherwise, these n+l-ro columns contain the matrix
%            V2 transformed as described in Step 3 of the TLS algorithm
%            (see Method).
%            If meth = 2, the columns of V whose index i corresponds
%            with NL(i) = 1, are the possibly transformed n+l-ro base
%            vectors of the right singular subspace corresponding to the
%            singular values of C = [A|B] which are less than or equal
%            to to (see below the description of the parameters ro and
%            to). Specifically, if l = 0 or ro = 0, these vectors are
%            indeed the base vectors above. Otherwise, these vectors
%            form the matrix V2, transformed as described in Step 4 of
%            the PTLS algorithm (see Method). The TLS solution is
%            computed from these vectors. The other columns contain no
%            useful information.
%   S      - (optional) if meth = 1, the min(m,n+l) singular values of
%            the matrix C, ordered such that S(1) >= S(2) >= ... >=
%            >= S(p-1) >= S(p) >= 0, where p = min(m,n+l).
%            If meth = 2, then S must not be specified as an output
%            parameter.
%   Q      - (optional) if meth = 2, the max(1,2*p-1) vector containing
%            the partially diagonalized bidiagonal matrix J computed
%            from C, at the moment that the desired singular subspace
%            has been found. Specifically, the leading p entries of Q
%            contain the diagonal elements q(1),q(2),...,q(p) and the
%            entries Q(p+1), ..., Q(2*p-1) contain the superdiagonal
%            elements e(1),...,e(p-1) of J.
%            If meth = 1, then Q must not be specified as an output
%            parameter.
%   NL     - (optional) if meth = 2, an (n+l) vector; the indices of
%            the elements of this vector with value 1 indicate the
%            columns in V containing the (transformed) base vectors of
%            the right singular subspace of C from which the TLS
%            solution has been computed.
%            If meth = 1, then NL must not be specified as an output
%            parameter.
%   ro     - (optional) if meth = 1, and job = 1 or job = 3, or if
%            meth = 2, and r < 0, then ro contains the computed
%            (effective) rank of the TLS approximation [A+DA|B+DB].
%            If meth = 1, and job = 0 or job = 2, or if meth = 2, and
%            r >= 0, the input value r may be changed by the gateway if
%            the r-th and the (r+1)-th singular values of C = [A|B] are
%            considered to be equal, or if the upper triangular matrix F
%            (as defined in Method) is (numerically) singular.
%   to     - (optional) if meth = 2 and r >= 0 on entry, then to
%            contains the computed bound such that precisely ro singular
%            values of C = [A|B] are greater than to + tol.
%            If meth = 2, and r < 0, then to = t.
%            If meth = 1, then to must not be specified as an output
%            parameter.
%   rcnd   - (optional) real scalar containing the reciprocal of the
%            condition number of the matrix F.
%
%   Method:
%   The methods used are extensions of the classical TLS algorithm.
%
%   Let [A|B] denote the matrix formed by adjoining the columns of B to
%   the columns of A on the right.
%
%   Total Least Squares (TLS) definition:
%
%     Given matrices A and B, find a matrix X satisfying
%
%          (A + DA) X = B + DB,
%
%     where A and DA are m-by-n matrices, B and DB are m-by-l matrices
%     and X is an n-by-l matrix.
%     The solution X must be such that the Frobenius norm of [DA|DB] is
%     a minimum and each column of B + DB is in the range of A + DA.
%     Whenever the solution is not unique, the minimum norm solution X
%     is singled out.
%
%   Define matrix C = [A|B] and s(i) as its i-th singular value for
%   i = 1,...,min(m,n+l). If m < n+l, then s(j) = 0 for j = m+1,...,n+l.
%
%   The Classical TLS algorithm proceeds as follows:
%
%   Step 1: Compute part of the singular value decomposition (SVD)
%           USV' of C = [A|B], namely compute S and V'. (An initial
%           QR factorization of C is used when m is larger enough
%           than n+l.)
%
%   Step 2: If the rank r is not given, compute the rank r0 of the data
%           [A|B] based on tol as follows: if job = 0 or job = 1,
%
%              s(1) >= ... >= s(r0) > tol*s(1) >= ... >= s(n+l).
%
%           Otherwise, tol can be computed from the standard deviation
%           sdev of the errors on [A|B]:
%
%              tol = sqrt(2 * max(m,n+l)) * sdev,
%
%           and the rank r0 is determined (if job = 1 or 3) using
%
%              s(1) >= ... >= s(r0) > tol >= ... >= s(n+l).
%
%           The rank r of the approximation [A+DA|B+DB] is then equal
%           to the minimum of n and r0.
%
%   Step 3: Let V2 be the matrix of the columns of V corresponding to
%           the (n+l-r) smallest singular values of C, i.e., the last
%           (n+l-r) columns of V. Using Householder transformations,
%           compute the orthogonal matrix Q such that:
%
%                     |VH   Y|
%            V2 x Q = |      |,                                      (1)
%                     |0    F|
%
%           where VH is an n-by-(n-r) matrix, Y is an n-by-l matrix and
%           F is an l-by-l upper triangular matrix. If F is singular,
%           then lower the rank r with the multiplicity of s(r) and
%           repeat this step.
%
%   Step 4: If F is nonsingular then the solution X is obtained by
%           solving the following equations by forward elimination:
%
%              X F = -Y.                                             (2)
%
%   Consider now the Partial Total Least Squares (PTLS) approach. Since
%   the TLS solution can be computed from any orthogonal basis of the
%   right singular subspace corresponding to the smallest singular
%   values of C, the Partial Singular Value Decomposition (PSVD) can be
%   used instead of the classical SVD. The dimension of this subspace
%   may be determined by the rank of C or by an upper bound for those
%   smallest singular values.
%
%   The PTLS algorithm proceeds as follows:
%
%   Step 1: Bidiagonalization phase
%    (a) If m is large enough than n + l, transform C into upper
%        triangular form R by Householder transformations.
%    (b) Transform C (or R) into upper bidiagonal form (p = min(m,n+l)):
%
%                   |q(1) e(1)  0   ...  0   |
%              (0)  | 0   q(2) e(2)      .   |
%             J   = | .                  .   |
%                   | .                e(p-1)|
%                   | 0             ... q(p) |
%
%        if m >= n + l, or lower bidiagonal form:
%
%                   |q(1)  0    0   ...  0     0   |
%              (0)  |e(1) q(2)  0        .     .   |
%             J   = | .                  .     .   |
%                   | .                 q(p)   .   |
%                   | 0             ... e(p-1) q(p)|
%
%        if m < n + l, using Householder transformations. In the second
%        case, transform the matrix to the upper bidiagonal form by
%        applying Givens rotations.
%    (c) Initialize the right singular base matrix with the identity
%        matrix.
%
%   Step 2: Partial diagonalization phase
%   If the upper bound t is not given, then compute t such that
%   precisely p - r singular values (p=min(m,n+l)) of the bidiagonal
%   matrix are less than or equal to t, using a bisection method.
%   Diagonalize the given bidiagonal matrix J partially, using either QL
%   iterations (if the upper left diagonal element of the considered
%   bidiagonal submatrix is smaller than the lower right diagonal
%   element) or QR iterations, such that J is split into unreduced
%   bidiagonal submatrices whose singular values are either all larger
%   than t or are all less than or equal to t. Accumulate the Givens
%   rotations in V.
%
%   Step 3: Back transformation phase
%   Apply the Householder transformations of Step 1(b) onto the base
%   vectors of V associated with the bidiagonal submatrices with all
%   singular values less than or equal to t.
%
%   Step 4: Computation of F and Y
%   Compute the factorization in (1). If F is singular, then reduce the
%   value of r by one and repeat Steps 2, 3 and 4.
%
%   Step 5: Computation of the TLS solution
%   If F is non-singular then the solution X is obtained by solving
%   the equations (2) by forward elimination.
%
%   Comments
%   1. The TLS solution is unique if r = n, F is nonsingular and
%   s(n) > s(n+1).
%   If F is singular, however, then the computed solution is infinite
%   and hence does not satisfy the second TLS criterion (see TLS
%   definition). These problems are called nongeneric and the TLS
%   computations have been generalized in order to solve them. The
%   generalization satisfies the TLS criteria for any number l of
%   observation vectors in B provided that, in addition, the solution
%   | X| is constrained to be orthogonal to all vectors of the form |w|
%   |-I|                                                            |0|
%   which belong to the space generated by the columns of the
%   submatrix |Y|.
%             |F|
%
%   2. If r is lowered in Step 4 of the PTLS algorithm, some additional
%   base vectors must be computed in Step 2. The additional computations
%   are kept to a minimum. If r is lowered in Step 4 but the
%   multiplicity of the r-th singular value is larger than 1, then the
%   value of r is further lowered with its multiplicity defined by the
%   parameter tol. This is done at the beginning of Step 2, using a
%   bisection method to estimate t.
%
%   3. The computational efficiency of the PTLS algorithm compared with
%   the classical TLS algorithm is obtained by making use of PSVD
%   instead of performing the entire SVD. Depending on the gap between
%   the r-th and the (r+1)-th singular values of C, the number (n+l-r)
%   of base vectors to be computed with respect to the column dimension
%   (n+l) of C and the desired accuracy reltol, the algorithm used by
%   the PTLS algorithm is approximately twice as fast as the classical
%   TLS algorithm at the expense of extra storage requirements, namely:
%     (n + l) x (n + l - 1)/2,  if m >= n + l or
%     m x (n + l - (m - 1)/2),  if m <  n + l.
%   This is because the Householder transformations performed on the
%   rows of C in the bidiagonalization phase (see Step 1) must be kept
%   until the end (Step 5).
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Sep. 2003.
%
% Revisions:
%   -
%
