% SKEWHAMILDEFLFZ.F - Gateway function for computing the eigenvalues of
%                     a complex skew-Hamiltonian/Hamiltonian pencil and
%                     the right deflating subspace and/or the companion
%                     subspace corresponding to the eigenvalues with
%                     strictly negative real part, using SLICOT routine
%                     MB03FZ.f.
%                     The gateway accepts real or complex input
%                     matrices.
%
%   [ALPHA,BETA(,Q,U,neig)] = ...
%                skewHamildeflfZ(Z,B,FG(,compq,compu(,orthm)))
%
%   SKEWHAMILDEFLFZ computes the eigenvalues of a complex n-by-n skew-Hamiltonian/
%   Hamiltonian pencil aS - bH, with n = 2m,
%
%              H  T           (  B  F  )       (  Z11  Z12  )
%       S = J Z  J  Z and H = (      H ), Z =: (            ).       (1)
%                             (  G -B  )       (  Z21  Z22  )
%
%   The structured Schur form of the embedded real skew-Hamiltonian/
%                                                         H  T
%   skew-Hamiltonian pencil, aB_S - bB_T, with B_S = J B_Z  J  B_Z,
%
%             (  Re(Z11)  -Im(Z11)  |  Re(Z12)  -Im(Z12)  )
%             (                     |                     )
%             (  Im(Z11)   Re(Z11)  |  Im(Z12)   Re(Z12)  )
%             (                     |                     )
%       B_Z = (---------------------+---------------------) ,
%             (                     |                     )
%             (  Re(Z21)  -Im(Z21)  |  Re(Z22)  -Im(Z22)  )
%             (                     |                     )
%             (  Im(Z21)   Re(Z21)  |  Im(Z22)   Re(Z22)  )
%                                                                    (2)
%             ( -Im(B)  -Re(B)  | -Im(F)  -Re(F)  )
%             (                 |                 )
%             (  Re(B)  -Im(B)  |  Re(F)  -Im(F)  )
%             (                 |                 )
%       B_T = (-----------------+-----------------) ,  T = i*H,
%             (                 |      T       T  )
%             ( -Im(G)  -Re(G)  | -Im(B )  Re(B ) )
%             (                 |      T       T  )
%             (  Re(G)  -Im(G)  | -Re(B ) -Im(B ) )
%
%   is determined and used to compute the eigenvalues. Optionally, if
%   compq = 1, an orthonormal basis of the right deflating subspace,
%   Def_-(S, H), of the pencil aS - bH in (1), corresponding to the
%   eigenvalues with strictly negative real part, is computed. Namely,
%   after transforming aB_S - bB_H, in the factored form, by unitary
%                                      H  T
%   matrices, we have B_Sout = J B_Zout  J  B_Zout,
%
%                ( Ao  Do  )              ( Bo  Fo  )
%       B_Zout = (         ) and B_Hout = (       H ),               (3)
%                (  0  Co  )              (  0 -Bo  )
%
%   and the eigenvalues with strictly negative real part of the
%   complex pencil aB_Sout - bB_Hout are moved to the top. Optionally,
%   if compu = 1, an orthonormal basis of the companion subspace,
%   range(P_U), which corresponds to the eigenvalues with negative
%   real part, is computed.
%
% Description of input parameters:
%   Z      - the n-by-n matrix Z.
%   B      - the m-by-m matrix B, with m = n/2.
%   FG     - an  m-by-(m+1) matrix containing the triangles of the
%            Hermitian matrices F and G, as follows:
%            the leading m-by-m lower triangular part contains the lower
%            triangle of the matrix G, and the m-by-m upper triangular
%            part of the submatrix in the columns 2 to m+1 contains the
%            upper triangle of the matrix F of H in (1).
%            So, if i >= j, then G(i,j) = -G(j,i) is stored in FG(i,j)
%            and F(j,i) = -F(i,j) is stored in FG(j,i+1).
%            FG is an empty matrix if m = 0.
%   compq  - (optional) scalar indicating whether to compute a basis Q
%            of the deflating subspace corresponding to the eigenvalues
%            of the pencil aS - bH with strictly negative real part.
%            = 0 :  Q is not required (default);
%            = 1 :  on exit, Q contains a matrix with orthonormal
%                   columns.
%   compu  - (optional) scalar indicating whether to compute a basis U
%            of the companion subspace corresponding to the eigenvalues
%            of the pencil aS - bH with strictly negative real part.
%            = 0 :  U is not required (default);
%            = 1 :  on exit, U contains a matrix with orthonormal
%                   columns.
%   orthm  - (optional) if compq = 1 or compu = 1, scalar indicating the
%            technique for computing the unitary basis of the deflating
%            subspace, or the companion subspace, as follows:
%            = 0 :  QR factorization, the fastest technique (default);
%            = 1 :  QR factorization with column pivoting;
%            = 2 :  singular value decomposition.
%            If compq = 0 and compu = 0, the orthm value is not used.
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
%   U      - if compu = 1, an n-by-neig matrix containing an orthonormal
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
