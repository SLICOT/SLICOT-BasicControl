% HNORM.F  - MEX-function function for computing various system norms
%            (Hankel norm, H2 norm) and complex stability radius of a
%            continuous-time or discrete-time system in standard form,
%            using SLICOT routines AB13AD, AB13BD, AB13ED and AB13FD.
%
%   [Hnorm(,ns,Ao,Bo,Co(,Do)(,hsv))]
%                              = Hnorm(job,A,B,C(,D)(,dico,equil
%                                      (,alpha)(,tol)))
%   [low(,high))]              = Hnorm(job,A(,tol))
%
%   [Hknorm(,ns,Ao,Bo,Co,hsv)] = Hnorm(1,A,B,C(,dico,equil,alpha))
%   [Hnorm(,nq,Ao,Bo,Co,Do)]   = Hnorm(2,A,B,C,D(,dico,jobn,tol))
%   [low(,high)]               = Hnorm(3,A(,tol))
%   [beta(,omega)]             = Hnorm(4,A(,tol))
%
%   1) HNORM computes the Hankel norm, H2 or L2 norm of the continuous-time
%   or discrete-time system in standard form
%                                     -1
%        G(lambda) = C*( lambda*I - A ) *B + D .
%
%   For computing the Hankel norm, or for computing the H2/L2 norm for
%   a continuous-time system, D is assumed 0. For computing the H2 or
%   L2 norm, the transfer-function matrix G must not have poles on the
%   imaginary axis, for a continuous-time system, or on the unit circle,
%   for a discrete-time system. If the H2-norm is computed, the system
%   must be stable.
%
%   2) HNORM estimates beta(A), the 2-norm distance from a real matrix A to
%   the nearest complex matrix with an eigenvalue on the imaginary axis.
%   The number beta(A) is the minimum of the smallest singular value of
%   the matrix (A - jwI), where I is the identity matrix and j**2 = -1,
%   and the minimum is taken over all real w. The estimate is given as
%
%            low <= beta(A) <= high,
%
%   where either
%
%            (1 + tol) * low >= high,
%   or
%            low = 0   and   high = delta,
%
%   and delta is a small number approximately equal to the square root
%   of machine precision times the Frobenius norm of A. If all
%   eigenvalues of A lie in the open left half complex plane, then
%   beta(A) is the distance to the nearest unstable complex matrix,
%   i.e., the complex stability radius.
%
%   Description of input parameters:
%   job    - option parameter indicating the task to be performed:
%            = 1 :  Hankel-norm of the alpha-stable projection of the
%                   transfer-function matrix G (D assumed 0);
%            = 2 :  H2 or L2 norm of a system;
%            = 3 :  complex stability radius, using bisection;
%            = 4 :  complex stability radius, using bisection and SVD.
%   A      - the n-by-n system state matrix A.
%   B      - the n-by-m system input matrix B.
%   C      - the p-by-n system output matrix C.
%   D      - the p-by-m system input/output matrix D.
%   dico   - (optional) specifies the type of the system:
%            = 1 :  continuous-time system;
%            = 2 :  discrete-time system.
%            Default: dico = 1.
%   equil  - (optional) specifies whether the system (A,B,C) should be
%            preliminarily equilibrated:
%            = 1 :  do not perform equilibration;
%            = 2 :  perform equilibration (scaling).
%            Default: equil = 1.
%   alpha  - (optional) constant specifying the alpha-stability boundary
%            for the eigenvalues of the state dynamics matrix A. For a
%            continuous-time system (dico = 1), alpha <= 0 is the
%            boundary value for the real parts of eigenvalues, while for
%            a discrete-time system (dico = 2), 0 <= alpha <= 1 is the
%            boundary value for the moduli of eigenvalues. The alpha-
%            stability domain is defined either as the open half complex
%            plane left to alpha, if dico = 1, or the interior of the
%            alpha-radius circle centered in the origin, if dico = 2.
%            Default: alpha =  -sqrt(epsilon_machine), if dico = 1;
%                     alpha = 1-sqrt(epsilon_machine), if dico = 2,
%            where epsilon_machine is the relative machine precision.
%   jobn   - (optional) constant specifying the norm to be computed:
%            = 1 :  the H2-norm;
%            = 2 :  the L2-norm.
%            Default: jobn = 1.
%   tol    - (optional) if job = 2, absolute tolerance level below
%            which the elements of B are considered zero (used for
%            controllability tests). If tol <= 0, default value is used.
%            Default: n*norm(B,1)*epsilon_machine.
%            If job = 3, tol specifies the accuracy with which low and
%            high approximate beta(A). If tol is less than
%            sqrt(epsilon_machine), then that value is used instead.
%            The recommended value is tol = 9, which gives an estimate
%            of beta(A) correct to within an order of magnitude.
%            Default: sqrt(epsilon_machine).
%            If job = 4, the accuracy with which beta(A) is to be
%            calculated. (See Comment 3 below.) If tol is set less than
%            epsilon_machine, then that value is used instead.
%            Default: epsilon_machine.
%
%   Description of output parameters:
%   Hknorm - the Hankel-norm of the alpha-stable projection of G.
%   Hnorm  - the H2-norm of G, if jobn = 1, or the L2-norm of G,
%            if jobn = 2.
%   ns     - the order of the alpha-stable subsystem.
%   nq     - the order of the resulting numerator Q of the right coprime
%            factorization with inner denominator of G (see Comment 1).
%            Generally, nq = n - nus, where nus is the number of
%            uncontrollable unstable eigenvalues.
%   Ao     - if job = 1, the n-by-n resulted state dynamics matrix in
%            a block diagonal real Schur form with its eigenvalues
%            reordered and separated. Ao has two diagonal blocks: the
%            leading ns-by-ns part has eigenvalues in the alpha-
%            stability domain and the trailing (n-ns)-by-(n-ns) part
%            has eigenvalues outside the alpha-stability domain.
%            If job = 2, the nq-by-nq state dynamics matrix (in a real
%            Schur form) of the numerator factor Q of the right coprime
%            factorization with inner denominator of G (see Comment 1).
%   Bo     - if job = 1, the n-by-m input matrix of the transformed
%            system.
%            If job = 2, the nq-by-m input matrix of the numerator
%            factor Q of the right coprime factorization with inner
%            denominator of G.
%   Co     - if job = 1, the p-by-n output matrix of the transformed
%            system.
%            If job = 2, the p-by-nq output matrix of the numerator
%            factor Q of the right coprime factorization with inner
%            denominator of G.
%   Do     - if job = 2, the p-by-m input/output matrix of the numerator
%            factor Q of the right coprime factorization with inner
%            denominator of G.
%   hsv    - (optional) the ns Hankel singular values (ordered
%            decreasingly) of the alpha-stable part of the system.
%            hsv(1) is the Hankel norm of the alpha-stable subsystem.
%   low    - a lower bound for beta(A).
%   high   - (optional) an upper bound for beta(A).
%   beta   - the computed value of beta(A), which actually is an upper
%            bound.
%   omega  - (optional) the value of w such that the smallest singular
%            value of (A - jwI) equals beta(A).
%
%   Comments
%   1) For job = 2, if the transfer-function matrix G is unstable,
%      then a right coprime factorization with inner denominator of G
%      is first computed
%               -1
%        G = Q*R  ,
%
%      where Q and R are stable transfer-function matrices and R is
%      inner. If G is stable, then Q = G and R = I.
%      Let (AQ,BQ,CQ,DQ) be the state-space representation of Q.
%
%      If dico = 1, then the L2-norm of G is computed as
%
%         norm2(G) = norm2(Q) = sqrt(trace(BQ'*X*BQ)),
%
%      where X satisfies the continuous-time Lyapunov equation
%
%         AQ'*X + X*AQ + CQ'*CQ = 0.
%
%      If dico = 2, then the l2-norm of G is computed as
%
%         norm2(G) = norm2(Q) = sqrt(trace(BQ'*X*BQ+DQ'*DQ)),
%
%      where X satisfies the discrete-time Lyapunov equation
%
%         AQ'*X*AQ - X + CQ'*CQ = 0.
%
%   2) For job = 3, the algorithm computes a lower bound low and an
%      upper bound high for beta(A) by a bisection method in the
%      following way. Given a non-negative real number sigma, the
%      Hamiltonian matrix H(sigma) is constructed:
%
%                     |   A      -sigma*I |     | A   G  |
%         H(sigma) =  |                   | :=  |        | .
%                     | sigma*I    -A'    |     | F  -A' |
%
%      It can be shown that H(sigma) has an eigenvalue whose real
%      part is zero if and only if sigma >= beta. Any lower and upper
%      bounds on beta(A) can be improved by choosing a number between
%      them and checking to see if H(sigma) has an eigenvalue with zero
%      real part.  This decision is made by computing the eigenvalues of
%      H(sigma) using the square reduced algorithm of Van Loan.
%
%   3) For job = 4, a reliable, quadratically convergent algorithm is
%      used, which finds (using a bisection strategy) an interval which
%      contains beta(A), and then switches to a modified bisection
%      strategy which converges quadratically to a minimizer. Note that
%      the efficiency of the strategy degrades if there are several
%      local minima that are near or equal the global minimum.
%      The computed function value beta satisfies
%
%         beta(A) <= beta + epsilon,
%
%         beta/(1+tol) - delta <= max(beta(A), sqrt(2*n*eps)*norm(A)),
%
%      where norm(A) is the Frobenius norm of A,
%
%         epsilon = p(n) * eps * norm(A),
%      and
%         delta   = p(n) * sqrt(eps) * norm(A),
%
%      eps = epsilon_machine, and p(n) is a low degree polynomial. It
%      is recommended to choose tol greater than sqrt(eps). Although
%      rounding errors can cause failure for smaller values of tol,
%      nevertheless, the calculation usually succeeds. Regardless of
%      success or failure, the first inequality holds.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2002.
%
% Revisions:
%   -
%
