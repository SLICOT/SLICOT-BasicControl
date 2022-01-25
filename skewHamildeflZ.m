% SKEWHAMILDEFLZ.F - Gateway function for computing the eigenvalues of a
%                    complex skew-Hamiltonian/Hamiltonian pencil and the
%                    right deflating subspace corresponding to the
%                    eigenvalues with strictly negative real part, using
%                    SLICOT routine MB03LZ.f.
%                    The gateway accepts real or complex input matrices.
%
%   [ALPHA,BETA(,Q,neig)] = skewHamildeflZ(A,DE,B,FG(,compq(,orthm)))
%
%   SKEWHAMILDEFLZ computes the eigenvalues of a complex n-by-n skew-Hamiltonian/
%   Hamiltonian pencil aS - bH, with n = 2m,
%
%         (  A  D  )         (  B  F  )
%     S = (      H ),    H = (      H ).                             (1)
%         (  E  A  )         (  G -B  )
%
%   The structured Schur form of the embedded real skew-Hamiltonian/
%   skew-Hamiltonian pencil aB_S - bB_T, defined as
%
%
%           (  Re(A)  -Im(A)  |  Re(D)  -Im(D)  )
%           (                 |                 )
%           (  Im(A)   Re(A)  |  Im(D)   Re(D)  )
%           (                 |                 )
%     B_S = (-----------------+-----------------) , and
%           (                 |      T       T  )
%           (  Re(E)  -Im(E)  |  Re(A )  Im(A ) )
%           (                 |      T       T  )
%           (  Im(E)   Re(E)  | -Im(A )  Re(A ) )
%                                                                    (2)
%           ( -Im(B)  -Re(B)  | -Im(F)  -Re(F)  )
%           (                 |                 )
%           (  Re(B)  -Im(B)  |  Re(F)  -Im(F)  )
%           (                 |                 )
%     B_T = (-----------------+-----------------) ,  T = i*H,
%           (                 |      T       T  )
%           ( -Im(G)  -Re(G)  | -Im(B )  Re(B ) )
%           (                 |      T       T  )
%           (  Re(G)  -Im(G)  | -Re(B ) -Im(B ) )
%
%   is determined and used to compute the eigenvalues. Optionally,
%   if compq = 1, an orthonormal basis of the right deflating
%   subspace of the pencil aS - bH, corresponding to the eigenvalues
%   with strictly negative real part, is computed. Namely, after
%   transforming aB_S - bB_H by unitary matrices, we have
%
%              ( Ao  Do  )              ( Bo  Fo  )
%     B_Sout = (       H ) and B_Hout = (       H ),                 (3)
%              (  0  Ao  )              (  0 -Bo  )
%
%   and the eigenvalues with strictly negative real part of the
%   complex pencil aB_Sout - bB_Hout are moved to the top.
%
% Description of input parameters:
%   A      - the m-by-m matrix A, with m = n/2.
%   DE     - an  m-by-(m+1) matrix containing the strict triangles of
%            the skew-Hermitian matrices D and E, as follows:
%            the leading m-by-m lower triangular part contains the lower
%            triangle of the matrix E, and the m-by-m upper triangular
%            part of the submatrix in the columns 2 to m+1 contains the
%            upper triangle of the matrix D of S in (1).
%            So, if i >= j, then E(i,j) = -E(j,i) is stored in DE(i,j)
%            and D(j,i) = -D(i,j) is stored in DE(j,i+1).
%            DE is an empty matrix if m = 0.
%   B      - the m-by-m matrix B.
%   FG     - an  m-by-(m+1) matrix containing the triangles of the
%            Hermitian matrices F and G, as follows:
%            the leading m-by-m lower triangular part contains the lower
%            triangle of the matrix G, and the m-by-m upper triangular
%            part of the submatrix in the columns 2 to m+1 contains the
%            upper triangle of the matrix F of H in (1).
%            So, if i >= j, then G(i,j) = -G(j,i) is stored in FG(i,j)
%            and F(j,i) = -F(i,j) is stored in FG(j,i+1).
%            FG is an empty matrix if m = 0.
%   compq  - (optional) scalar indicating whether an orthonormal basis Q
%            of the right deflating subspace is returned, or if Q is not
%            required, as follows:
%            = 0 :  Q is not required (default);
%            = 1 :  on exit, Q contains a matrix with orthonormal
%                   columns.
%   orthm  - (optional) if compq = 1, scalar indicating the technique
%            for computing the unitary basis of the deflating subspace,
%            as follows:
%            = 0 :  QR factorization, the fastest technique (default);
%            = 1 :  QR factorization with column pivoting;
%            = 2 :  singular value decomposition.
%            If compq = 0, the orthm value is not used.
%            Usually, orthm = 0 gives acceptable results, but badly
%            scaled or ill-conditioned problems might need to set
%            orthm = 1 or even orthm = 2.
%
% Description of output parameters:
%   ALPHA  - the n-vector of the numerators alpha defining the
%            eigenvalues of the pencil aS - bT.
%   BETA     the n-vector of the denominators beta defining the
%            eigenvalues of the pencil aS - bT.
%            Together, the quantities alpha = ALPHA(j) and
%            beta = BETA(j) represent the j-th eigenvalue of the pencil
%            aS - bT, in the form lambda = alpha/beta. Since lambda may
%            overflow, the ratios should not, in general, be computed.
%   Q      - if compq = 1, an n-by-neig matrix containing an orthonormal
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
