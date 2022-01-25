% SKEWHAMILEIGZ.F - Gateway function for computing the eigenvalues of a
%                   complex skew-Hamiltonian/Hamiltonian pencil, using
%                   SLICOT routine MB04BZ.
%                   The gateway accepts real or complex input matrices.
%
%   [ALPHA,BETA(,Ao,Do,Bo,Fo(,Q))] =
%                              skewHamileigZ(A,DE,B,FG(,job,compq))
%
%   [ALPHA,BETA]             = skewHamileigZ(A,DE,B,FG)
%   [ALPHA,BETA,Ao,Do,Bo,Fo] = skewHamileigZ(A,DE,B,FG,1)
%
%   SKEWHAMILEIGZ computes the eigenvalues of a complex N-by-N skew-Hamiltonian/
%   Hamiltonian pencil aS - bH, with
%
%         (  A  D  )         (  B  F  )
%     S = (      H ) and H = (      H ).                           (1)
%         (  E  A  )         (  G -B  )
%
%   The structured Schur form of the embedded real skew-Hamiltonian/
%   skew-Hamiltonian pencil aB_S - bB_T, defined as
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
%                                                                  (2)
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
%   if JOB = 'T', the pencil aB_S - bB_H (B_H = -i*B_T) is transformed
%   by a unitary matrix Q to the structured Schur form
%
%              ( Ao  Do  )              ( Bo  Fo  )
%     B_Sout = (       H ) and B_Hout = (       H ),               (3)
%              (  0  Ao  )              (  0 -Bo  )
%
%   where Ao and Bo are upper triangular, Do is skew-Hermitian, and
%   Fo is Hermitian. The embedding doubles the multiplicities of the
%   eigenvalues of the pencil aS - bH. Optionally, if COMPQ = 'C', the
%   unitary matrix Q is computed.
%
% Description of input parameters:
%   A      - the m-by-m matrix A, with m = n/2.
%   DE     - an  m-by-(m+1) matrix containing the triangles of the
%            skew-Hermitian matrices D and E, as follows:
%            the leading m-by-m lower triangular part contains the
%            lower triangle of the matrix E, and the m-by-m upper
%            triangular part of the submatrix in the columns 2 to m+1
%            contains the upper triangle of the matrix D of S in (1).
%            So, if i >= j, then E(i,j) = -conj(E(j,i)) is stored in
%            DE(i,j) and D(j,i) = -conj(D(i,j)) is stored in DE(j,i+1).
%            DE is an empty matrix if m = 0.
%   B      - the m-by-m matrix B.
%   FG     - an  m-by-(m+1) matrix containing the triangles of the
%            Hermitian matrices F and G, as follows:
%            the leading m-by-m lower triangular part contains the lower
%            triangle of the matrix G, and the m-by-m upper triangular
%            part of the submatrix in the columns 2 to m+1 contains the
%            upper triangle of the matrix F of H in (1).
%            So, if i >= j, then G(i,j) = conj(G(j,i)) is stored in
%            FG(i,j) and F(j,i) = conj(F(i,j)) is stored in FG(j,i+1).
%            FG is an empty matrix if m = 0.
%   job    - (optional) scalar indicating the computation to be
%            performed, as follows:
%            = 0 :  compute the eigenvalues only (default);
%            = 1 :  compute the eigenvalues and the matrices of the
%                   transformed pencil in (3).
%   compq  - (optional) scalar indicating whether the unitary
%            transformation matrix Q is required, as follows:
%            = 0 :  Q is not required (default);
%            = 1 :  on exit, Q contains the unitary transformation
%                   matrix.
%
% Description of output parameters:
%   ALPHA  - the n-vector of the numerators alpha defining the
%            eigenvalues of the pencil aS - bH.
%   BETA     the n-vector of the denominators beta defining the
%            eigenvalues of the pencil aS - bH.
%            Together, the quantities alpha = ALPHA(j) and
%            beta = BETA(j) represent the j-th eigenvalue of the pencil
%            aS - bH, in the form lambda = alpha/beta. Since lambda may
%            overflow, the ratios should not, in general, be computed.
%   Ao     - if job = 1, the computed n-by-n submatrix Ao in (3).
%            The strictly lower triangular part is not zeroed.
%   Do     - if job = 1, the computed n-by-n upper triangular part of
%            the skew-Hermitian matrix Do in (3). The strictly lower
%            triangular part is not set.
%   Bo     - if job = 1, the computed n-by-n submatrix Bo in (3).
%            The strictly lower triangular part is not zeroed.
%   Fo     - if job = 1, the computed n-by-n upper triangular part of
%            the Hermitian matrix Fo in (3). The strictly lower
%            triangular part is not set.
%   Q      - if compq = 1, an 2*n-by-2*n matrix containing the computed
%            unitary matrix.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, July 2012.
%
% Revisions:
%   V. Sima, Nov. 2012.
%
