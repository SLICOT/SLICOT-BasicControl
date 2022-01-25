% SKEWHAMILDEFL.f - Gateway function for computing the eigenvalues of a
%                   skew-Hamiltonian/Hamiltonian pencil and the right
%                   deflating subspace corresponding to the eigenvalues
%                   with strictly negative real part, using SLICOT
%                   routine MB03LD.
%
%   [ALPHAR,ALPHAI,BETA(,Q,neig)] = ...
%                  skewHamildefl(A,DE,B,FG(,compq(,orthm)))
%
%   SKEWHAMILDEFL computes the eigenvalues of an n-by-n skew-Hamiltonian/Hamiltonian
%   pencil aS - bH, with
%
%         (  A  D  )         (  B  F  )
%     S = (      T ) and H = (      T ),                           (1)
%         (  E  A  )         (  G -B  )
%
%   using the structured Schur form of an embedded pencil.
%   Optionally, the right deflating subspace of aS - bH corresponding to
%   the eigenvalues with strictly negative real part is computed.
%
% Description of input parameters:
%   A      - the m-by-m matrix A, with m = n/2.
%   DE     - an  m-by-(m+1) matrix containing the strict triangles of
%            the skew-symmetric matrices D and E, as follows:
%            the leading m-by-m strictly lower triangular part contains
%            the strictly lower triangle of the matrix E, and the
%            m-by-m strictly upper triangular part of the submatrix in
%            the columns 2 to m+1 contains the strictly upper triangle
%            of the matrix D of S in (1).
%            So, if i > j, then E(i,j) = -E(j,i) is stored in DE(i,j)
%            and D(j,i) = -D(i,j) is stored in DE(j,i+1).
%            The entries on the diagonal and the first superdiagonal of
%            DE need not be set, but are assumed to be zero.
%            DE is an empty matrix if m = 0.
%   B      - the m-by-m matrix B.
%   FG     - an  m-by-(m+1) matrix containing the triangles of the
%            symmetric matrices F and G, as follows:
%            the leading m-by-m lower triangular part contains the lower
%            triangle of the matrix G, and the m-by-m upper triangular
%            part of the submatrix in the columns 2 to m+1 contains the
%            upper triangle of the matrix F of H in (1).
%            So, if i >= j, then G(i,j) = G(j,i) is stored in FG(i,j)
%            and F(j,i) = F(i,j) is stored in FG(j,i+1).
%            FG is an empty matrix if m = 0.
%   compq  - (optional) scalar indicating whether the orthogonal
%            transformation matrix Q is returned, or if Q is not
%            required, as follows:
%            = 0 :  Q is not required (default);
%            = 1 :  on exit, Q contains a matrix Q with orthogonal
%                   columns.
%   orthm  - (optional) if compq = 1, scalar indicating the technique
%            for computing the orthogonal basis of the deflating
%            subspace, as follows:
%            = 0 :  QR factorization, the fastest technique (default);
%            = 1 :  QR factorization with column pivoting;
%            = 2 :  singular value decomposition.
%            If compq = 0, the orthm value is not used.
%            Usually, orthm = 0 gives acceptable results, but badly
%            scaled or ill-conditioned problems might need to set
%            orthm = 1 or even orthm = 2.
%
% Description of output parameters:
%   ALPHAR,- the m-vectors of real parts and imaginary parts,
%   ALPHAI   respectively, of each scalar alpha defining an eigenvalue
%            of the pencil aS - bH.
%            If ALPHAI(j) is zero, then the j-th eigenvalue is real.
%   BETA     the m-vector of the scalars beta that define the eigenvalues
%            of the pencil aS - bH.
%            Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
%            beta = BETA(j) represent the j-th eigenvalue of the pencil
%            aS - bH, in the form lambda = alpha/beta. Since lambda may
%            overflow, the ratios should not, in general, be computed.
%            Due to the skew-Hamiltonian/Hamiltonian structure of the
%            pencil, for every eigenvalue lambda, -lambda is also an
%            eigenvalue. Only eigenvalues with imaginary parts greater
%            than or equal to zero are stored; their conjugate
%            eigenvalues are not stored. If imaginary parts are zero
%            (i.e., for real eigenvalues), only positive eigenvalues
%            are stored.
%   Q      - if compq = 1, an n-by-neig matrix containing the computed
%            basis of the right deflating subspace.
%   neig   - if compq = 1, the number of eigenvalues with negative real
%            parts.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, July 2012.
%
% Revisions:
%   -
%
