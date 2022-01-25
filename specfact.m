% SPECFACT.F - MEX-function for computing the spectral factorization
%              of a real polynomial, arising from optimality problems,
%              using SLICOT routines SB08MD and SB08ND.
%
%   [E(,res,B)] = specfact(task,A(,form))
%
%   1. SPECFACT computes a real polynomial E(s) such that
%
%        (a)  E(-s) * E(s) = A(-s) * A(s) and
%        (b)  E(s) is stable - that is, all the zeros of E(s) have
%             non-positive real parts,
%
%   which corresponds to computing the spectral factorization of the
%   real polynomial A(s) arising from continuous optimality problems.
%
%   The input polynomial may be supplied either in the form
%
%   A(s) = a(0) + a(1) * s + ... + a(DA) * s**DA
%
%   or as
%
%   B(s) = A(-s) * A(s)
%        = b(0) + b(1) * s**2  + ... + b(DA) * s**(2*DA)             (1)
%
%   2. SPECFACT computes a real polynomial E(z) such that
%
%        (a)  E(1/z) * E(z) = A(1/z) * A(z) and
%        (b)  E(z) is stable - that is, E(z) has no zeros with modulus
%             greater than 1,
%
%   which corresponds to computing the spectral factorization of the
%   real polynomial A(z) arising from discrete optimality problems.
%
%   The input polynomial may be supplied either in the form
%
%   A(z) = a(0) + a(1) * z + ... + a(DA) * z**DA
%
%   or as
%
%   B(z) = A(1/z) * A(z)
%        = b(0) + b(1) * (z + 1/z) + ... + b(DA) * (z**DA + 1/z**DA) (2)
%
%   Description of input parameters:
%   task   - integer specifying the computations to be performed.
%            = 1 :  compute the spectral factorization for the
%                   continuous optimality problems;
%            = 2 :  compute the spectral factorization for the
%                   discrete optimality problems.
%   A      - a (DA+1)-vector containing either the coefficients of the
%            polynomial A in increasing powers of the indeterminate s
%            or z, if form = 1, or the coefficients b(0), ..., b(DA) of
%            the polynomial B in the formulas (1) or (2), if form = 2.
%   form   - (optional) indicates whether the coefficients of A or B are
%            supplied, as follows:
%            = 1 :  the coefficients of A are supplied;
%            = 2 :  the coefficients of B are supplied.
%            Default: form = 1.
%
%   Description of output parameters:
%   E      - a (DA+1)-vector containing the coefficients of the spectral
%            factor E in increasing powers of s or z.
%   res    - (optional) an estimate of the accuracy with which the
%            coefficients of the polynomial E have been computed (see
%            also Method and Comments).
%   B      - (optional) the (DA+1)-vector containing the coefficients
%            b(0), ..., b(DA) of the polynomial B in the formulas (1)
%            or (2).
%
%   Method
%       _                                               _
%   Let A(s) be the conjugate polynomial of A(s), i.e., A(s) = A(-s), if
%                     _
%   task = 1, and let A(z) be the conjugate polynomial of A(z), i.e., 
%   _
%   A(z) = A(1/z), if task = 2. The method used is based on applying the
%   Newton-Raphson iteration to the function
%             _       _
%      F(e) = A * A - e * e,
%
%   which leads to the iteration formulae:
%
%      _(i)   (i)  _(i)   (i)     _      )
%      q   * x   + x   * q    = 2 A * A  )
%                                        )   for i = 0, 1, 2,...
%       (i+1)    (i)   (i)               )
%      q     = (q   + x   )/2            )
%
%                                (0)         DA
%   For task = 1, starting from q   = (1 + s)   (which has no zeros in
%                                                           (1)   (2)
%   the closed right half-plane), the sequence of iterates q   , q   ,
%   ..., converges to a solution of F(e) = 0 which has no zeros in the
%   open right half-plane. Similarly, for task = 2, the iteration starts
%   from
%
%         (0)                                        DA
%        q   (z) = (b(0) + b(1) * z + ... + b(DA) * z  ) / SQRT( b(0))
%
%   which is a Hurwitz polynomial that has no zeros in the closed unit
%                     (i)
%   circle. Then lim q   = e, the convergence is uniform and e is a
%   Hurwitz polynomial.
%
%   The iterates satisfy the following conditions, if task = 1:
%
%            (i)
%      (a)  q   is a stable polynomial (no zeros in the closed right
%           half-plane) and
%
%            (i)        (i-1)
%      (b)  q   (1) <= q     (1),
%
%   or, if task = 2:
%            (i)
%      (a)  q    has no zeros in the closed unit circle,
%            (i)        (i-1)
%      (b)  q   (0) <= q     (0) and
%
%            DA   (i) 2    DA     2
%      (c)  SUM (q   )  - SUM (A )  >= 0.
%           k=0   k       k=0   k
%                                   (i)
%   The iterative process stops if q    violates (a) or (b) (if
%   task = 1), or (a), (b) or (c) (if task = 2), or if the condition
%                     _(i) (i)  _
%      (d)  RES  = ||(q   q   - A A)|| < tol,
%
%   is satisfied, where || . || denotes the largest coefficient of
%                   _(i) (i)  _
%   the polynomial (q   q   - A A) and tol is an estimate of the
%                                                  _(i)  (i)
%   rounding error in the computed coefficients of q    q   .  If
%                                          (i-1)
%   condition (a) or (b) is violated then q      is taken; if
%                                   (i)
%   condition (c) is violated then q    is used. If there is no 
%   convergence after 30 iterations, then an error indicator is set, and
%   the value of res may indicate whether or not the last computed
%   iterate is close to the solution.
%
%   If task = 1 and form = 2, then it is possible that the equation
%   e(-s) * e(s) = B(s) has no real solution, which will be the case
%   if A(1) < 0 or if ( -1)**DA * A(DA+1) < 0.
%                                                       (0)
%   If task = 2 and form = 2, then it is possible that q    is not a
%   Hurwitz polynomial, in which case the equation e(1/z) * e(z) = B(z)
%   has no real solution.
%
%   Comments:
%   1. In order for the problem e(-s) * e(s) = B(s) (for task = 1) to
%   have a real solution e(s), it is necessary and sufficient that
%   B(j*omega) >= 0 for any purely imaginary argument j*omega.
%   2. The conditioning of the problem if task = 1 depends upon the
%   distance of the zeros of A(s) from the imaginary axis and on their
%   multiplicity. For a well-conditioned problem the accuracy of the
%   computed coefficients of E(s) is of the order of res. However, for
%   problems with zeros near the imaginary axis or with multiple zeros,
%   the value of res may be an overestimate of the true accuracy.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2003.
%
% Revisions:
%   -
