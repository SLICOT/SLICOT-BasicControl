% SKEWHAMILDEFLF.F - Gateway function for computing the eigenvalues of a
%                    skew-Hamiltonian/Hamiltonian pencil, the right
%                    deflating subspace and the companion subspace
%                    corresponding to the eigenvalues with strictly
%                    negative real part, using SLICOT routine MB03LF.
%
%   [ALPHAR,ALPHAI,BETA(,Q,U,neig)] = skewHamildeflf(Z,B,FG(,...
%                                              compq,compu(,orthm)))
%
%   SKEWHAMILDEFLF computes the relevant eigenvalues of a real n-by-n skew-
%   Hamiltonian/Hamiltonian pencil aS - bH, with
%
%                                 (  B  F  )      (  0  I  )
%     S = T Z = J Z' J' Z and H = (        ), J = (        ),        (1)
%                                 (  G -B' )      ( -I  0  )
%
%   where the notation M' denotes the transpose of the matrix M.
%   Optionally, if compq = 1, an orthogonal basis of the right deflating
%   subspace of aS - bH corresponding to the eigenvalues with strictly
%   negative real part is computed. Optionally, if compu = 1, an
%   orthogonal basis of the companion subspace, range(P_U), which
%   corresponds to the eigenvalues with strictly negative real part, is
%   computed.
%
% Description of input parameters:
%   Z      - the n-by-n matrix Z.
%   B      - the m-by-m matrix B, with m = n/2.
%   FG     - an  m-by-(m+1) matrix containing the triangles of the
%            symmetric matrices F and G, as follows:
%            the leading m-by-m lower triangular part contains the lower
%            triangle of the matrix G, and the m-by-m upper triangular
%            part of the submatrix in the columns 2 to m+1 contains the
%            upper triangle of the matrix F of H in (1).
%            So, if i >= j, then G(i,j) = G(j,i) is stored in FG(i,j)
%            and F(j,i) = F(i,j) is stored in FG(j,i+1).
%            FG is an empty matrix if m = 0.
%   compq  - (optional) scalar indicating whether the orthogonal basis Q
%            of the right deflating subspace corresponding to the
%            eigenvalues with strictly negative real part is returned,
%            or if Q is not required, as follows:
%            = 0 :  Q is not required (default);
%            = 1 :  on exit, Q contains a matrix Q with orthogonal
%                   columns.
%   compu  - (optional) scalar indicating whether the orthogonal basis U
%            of the companion subspace corresponding to the eigenvalues
%            with strictly negative real part is returned, or if U is
%            not required, as follows:
%            = 0 :  U is not required (default);
%            = 1 :  on exit, U contains a matrix U with orthogonal
%                   columns.
%   orthm  - (optional) if compq = 1 or compu = 1, scalar indicating the
%            technique for computing the orthogonal basis of the
%            deflating subspace, or the companion subspace, as follows:
%            = 0 :  QR factorization, the fastest technique (default);
%            = 1 :  QR factorization with column pivoting;
%            = 2 :  singular value decomposition.
%            If compq = 0 and compu = 0, the orthm value is not used.
%            Usually, orthm = 0 gives acceptable results, but badly
%            scaled or ill-conditioned problems might need to set
%            orthm = 1 or even orthm = 2.
%
% Description of output parameters:
%   ALPHAR,- the m-vectors of real parts and imaginary parts,
%   ALPHAI   respectively, of each scalar alpha defining an eigenvalue
%            of the pencil aS - bH.
%            If ALPHAI(j) is zero, then the j-th eigenvalue is real.
%   BETA     the m-vector of the scalars beta that define the
%            eigenvalues of the pencil aS - bH.
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
%   U      - if compu = 1, an n-by-neig matrix containing the computed
%            basis of the companion subspace.
%   neig   - if compq = 1 or compu = 1, the number of eigenvalues with
%            negative real parts.
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
