% POLEZEROZ.F - MEX-function for computing the poles and zeros
%               of a standard or descriptor system, and the Kronecker
%               structure of the system pencil, using SLICOT routines
%               AB08MZ, AB08NZ, and AG08BZ. Complex systems.
%
%   [(poles(,betap))(,zeros,betaz)(,rnk) ...
%                            (,infz,Kronl,Kronr(,infe,niz,Af,Ef))] ...
%                           = polezeroz(task,A(,B,C,D)(,E)(,tol,bal))
%   [poles(,betap)]         = polezeroz(0,A(,E))
%   [rnk]                   = polezeroz(1,A,B,C,D(,tol,bal))
%   [rnk,infz,Kronl,Kronr(,infe,niz,Af,Ef)] = ...
%                           = polezeroz(2,A,B,C,D(,E)(,tol,bal))
%   [zeros,betaz,rnk,infz,Kronl,Kronr(,infe,niz,Af,Ef)] ...
%                           = polezeroz(3,A,B,C,D(,E)(,tol,bal))
%   [poles(,betap),zeros,betaz,rnk,infz,Kronl, ...
%                             Kronr(,infe,niz,Af,Ef)] ...
%                           = polezeroz(4,A,B,C,D(,E)(,tol,bal))
%
%   POLEZERO computes the normal rank, poles, zeros, and the Kronecker
%   structure of the system pencil for a standard or descriptor system.
%
%   Description of input parameters:
%   task   - integer specifying the computations to be performed.
%            = 0 :  compute the poles of the system;
%            = 1 :  compute the normal rank of the transfer-function
%                   matrix of a standard system;
%            = 2 :  compute the normal rank and Kronecker structure of
%                   the system;
%            = 3 :  compute the normal rank, Kronecker structure, and
%                   zeros of the system;
%            = 4 :  compute the normal rank, Kronecker structure, poles,
%                   and zeros of the system.
%   A      - the l-by-n state dynamics matrix A.
%            If task = 0, task = 1, or task = 4, it is assumed that
%            l = n.
%   B      - the l-by-m input/state matrix B.
%            If task = 0, B must not be specified as an input parameter.
%   C      - the p-by-n state/output matrix C.
%            If task = 0, C must not be specified as an input parameter.
%   D      - the p-by-m input/output matrix D.
%            If task = 0, D must not be specified as an input parameter.
%   E      - (optional) the l-by-n descriptor matrix E.
%            For a standard system, one must have l = n, and matrix E
%            may not be specified as an input parameter, or it may be
%            set to either an empty matrix or an identity matrix of
%            order n.
%            If task > 1, l = n = 1, and the number of input parameters
%            is 6 or larger, then the 6-th parameter is taken as E,
%            not tol.
%   tol    - (optional) if task > 0, real scalar containing the
%            tolerance to be used in rank decisions to determine the
%            effective rank, which is defined as the order of the
%            largest leading (or trailing) triangular submatrix in the
%            QR (or RQ) factorization with column (or row) pivoting
%            whose estimated condition number is less than 1/tol.
%            If tol <= 0, then default tolerances are used instead,
%            defined in terms of the size of the system matrix and eps,
%            where eps is the machine precision.
%            Default: tol = 0.
%            If task = 0, tol must not be specified as an input
%            parameter.
%   bal    - (optional) if task > 0, integer indicating whether the
%            system should be balanced (scaled).
%            = 0 :  use balancing;
%            = 1 :  do not use balancing.
%            Default: bal = 0.
%            If task = 0, bal must not be specified as an input
%            parameter.
%
%   Description of output parameters:
%   poles  - (optional) if task = 0 or task = 4, an n-vector containing
%            the (alpha part, if E is general, of the) system poles.
%            Otherwise, poles must not be specified as an output
%            parameter.
%   betap  - (optional) if E is general, and task = 0 or task = 4, an
%            n-vector containing the beta part of the system poles.
%            Otherwise, betap must not be specified as an output
%            parameter.
%   zeros  - (optional) if task >= 3, a vector containing the alpha part
%            of the system zeros.
%            Otherwise, zeros must not be specified as an output
%            parameter.
%   betaz  - (optional) if task >= 3, a vector containing the beta part
%            of the system zeros.
%            Otherwise, betaz must not be specified as an output
%            parameter.
%   rnk    - (optional) if task > 0, the normal rank of the system
%            pencil.
%            Otherwise, rnk must not be specified as an output
%            parameter.
%   infz   - (optional) if task >= 2, integer vector of length dinfz
%            containing information on the infinite elementary divisors
%            as follows: the system has infz(i) infinite elementary
%            divisors of degree i (in the Smith form), where
%            i = 1,2,...,dinfz.
%            Otherwise, infz must not be specified as an output
%            parameter.
%   Kronl  - (optional) if task >= 2, integer vector containing the left
%            Kronecker (row) indices.
%            Otherwise, Kronl must not be specified as an output
%            parameter.
%   Kronr  - (optional) if task >= 2, integer vector containing the
%            right Kronecker (column) indices.
%            Otherwise, Kronr must not be specified as an output
%            parameter.
%   infe   - (optional) if task >= 2, integer vector containing the
%            multiplicities of infinite eigenvalues.
%            Otherwise, infe must not be specified as an output
%            parameter.
%   niz    - (optional) if task >= 2, the number of infinite zeros.
%            Otherwise, niz must not be specified as an output
%            parameter.
%   Af     - (optional) if task >= 2, the nfz-by-nfz matrix Af of the
%            reduced pencil (see Method). If task >= 3, Af is in a
%            Schur form.
%            Otherwise, Af must not be specified as an output
%            parameter.
%   Ef     - (optional) if task >= 2, the nfz-by-nfz matrix Ef of the
%            reduced pencil (see Method). If task >= 3, Ef is in an
%            upper triangular form.
%            Otherwise, Ef must not be specified as an output
%            parameter.
%
%   Method:
%   If task = 0, the (generalized) eigenvalues of the system (pencil)
%   are computed.
%
%   If task = 1, the (n+p)-by-(m+n) compound matrix (B  A) is reduced
%                                                   (D  C)
%
%   to one with the same invariant zeros and with D of full row rank.
%   The normal rank of the transfer-function matrix is the rank of D.
%
%   If task >= 2 and identity E, the gateway extracts from the system
%   matrix of a state-space system (A,B,C,D) a regular pencil
%   Af - lambda*Ef which has the invariant zeros of the system as
%   generalized eigenvalues, as follows:
%
%        (a) construct the (n+p)-by-(m+n) compound matrix (B  A);
%                                                         (D  C)
%
%        (b) reduce the above system to one with the same invariant
%            zeros and with D of full row rank;
%
%        (c) pertranspose the system;
%
%        (d) reduce the system to one with the same invariant zeros and
%            with D square invertible;
%
%        (e) perform a unitary transformation on the columns of
%            (A - lambda*I  B) in order to reduce it to
%            (      C       D)
%
%            (Af - lambda*Ef  X)
%            (                 ), with Y and Ef square invertible;
%            (      0         Y)
%
%        (f) compute the right and left Kronecker indices of the system
%            (A,B,C,D), which together with the orders of the infinite
%            zeros (determined by steps (a) - (e)) constitute the
%            complete set of structural invariants under strict
%            equivalence transformations of a linear system.
%
%   If task >= 2 and general E, the gateway extracts from the system
%   matrix of a descriptor system (A-lambda*E,B,C,D) a regular pencil
%   Af-lambda*Ef which has the finite zeros of the system as generalized
%   eigenvalues. The procedure has the following computational steps:
%
%        (a) construct the (l+p)-by-(m+n) system pencil
%
%             S(lambda) = ( B  A )-lambda*( 0  E );
%                         ( D  C )        ( 0  0 )
%
%        (b) reduce S(lambda) to S1(lambda) with the same finite
%            zeros and right Kronecker structure but with E upper
%            triangular and nonsingular;
%
%        (c) reduce S1(lambda) to S2(lambda) with the same finite zeros
%            and right Kronecker structure but with D of full row rank;
%
%        (d) reduce S2(lambda) to S3(lambda) with the same finite zeros
%            and with D square invertible;
%
%        (e) perform a unitary transformation on the columns of
%
%            S3(lambda) = (A-lambda*E   B) in order to reduce it to
%                         (     C       D)
%
%            (Af-lambda*Ef   X), with Y and Ef square invertible;
%            (     0         Y)
%
%        (f) compute the right and left Kronecker indices of the system
%            matrix, which together with the multiplicities of the
%            finite and infinite eigenvalues constitute the
%            complete set of structural invariants under strict
%            equivalence transformations of a linear system.
%
%   Comments
%   1. If only the normal rank of the system pencil of a standard system
%      is desired, using task = 1 will ensure the maximum efficiency.
%   2. If E is general, the poles and zeros are returned as generalized
%      eigenvalues, with complex numerators alpha and real denominators
%      beta, but the ratios are not computed.
%   3. If Af and Ef are not needed, they should not be specified as
%      output arguments for task >= 3, since then the zeros are computed
%      slightly faster.

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2008.
%
% Revisions:
%   V. Sima, May 2009.
%
