% SKEWHAMIL2EIG.F - Gateway function for MB04FD.
%
%   [ALPHAR,ALPHAI,BETA(,Ao,Do,Bo,Fo(,Q))] =
%                      skewHamil2eig(A,DE,B,FG(,job,compq(,Q)))
%
%   [ALPHAR,ALPHAI,BETA]             = skewHamil2eig(A,DE,B,FG)
%   [ALPHAR,ALPHAI,BETA,Ao,Do,Bo,Fo] = skewHamil2eig(A,DE,B,FG,1)
%
%   SKEWHAMIL2EIG computes the eigenvalues of a real n-by-n 
%   skew-Hamiltonian/skew-Hamiltonian pencil aS - bT, with n = 2m,
%
%         (  A  D  )         (  B  F  )
%     S = (        ) and T = (        ).                             (1)
%         (  E  A' )         (  G  B' )
%
%   Optionally, if job = 1, decompositions of S and T will be
%   computed via an orthogonal transformation Q as follows:
%
%                   (  Ao  Do  )
%     J Q' J' S Q = (          ),
%                   (   0  Ao' )
%
%                   (  Bo  Fo  )            (  0  I  )
%     J Q' J' T Q = (          ), where J = (        ),              (2)
%                   (   0  Bo' )            ( -I  0  )
%
%   and Ao is upper triangular, Bo is upper quasi-triangular, and Do and
%   Fo are skew-symmetric. The notation M' denotes the transpose of the
%   matrix M.
%   Optionally, if compq = 1, the orthogonal transformation matrix Q
%   will be computed.
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
%            skew-symmetric matrices F and G, as follows:
%            the leading m-by-m strictly lower triangular part contains
%            the strictly lower triangle of the matrix G, and the
%            m-by-m strictly upper triangular part of the submatrix in
%            the columns 2 to m+1 contains the strictly upper triangle
%            of the matrix F of T in (1).
%            So, if i > j, then G(i,j) = -G(j,i) is stored in FG(i,j)
%            and F(j,i) = -F(i,j) is stored in FG(j,i+1).
%            The entries on the diagonal and the first superdiagonal of
%            FG need not be set, but are assumed to be zero.
%            FG is an empty matrix if m = 0.
%   job    - (optional) scalar indicating the computation to be
%            performed, as follows:
%            = 0 :  compute the eigenvalues only (default);
%            = 1 :  compute the eigenvalues and the matrices of the
%                   transformed pencil in (2).
%   compq  - (optional) scalar indicating whether the orthogonal
%            transformation matrix Q is returned or accumulated
%            into an orthogonal matrix, or if Q is not required,
%            as follows:
%            = 0 :  Q is not required (default);
%            = 1 :  on entry, Q need not be set;
%                   on exit, Q contains the orthogonal matrix Q;
%            = 2 :  the orthogonal transformations are accumulated
%                   into Q;
%                   on input, Q must contain an orthogonal matrix Q0;
%                   on exit, Q contains Q0*Q.
%   Q      - if compq = 2, the n-by-n orthogonal matrix Q0.
%
% Description of output parameters:
%   ALPHAR,- the m-vectors of real parts and imaginary parts,
%   ALPHAI   respectively, of each scalar alpha defining an eigenvalue
%            of the pencil aS - bT.
%            If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
%            positive, then the j-th and (j+1)-st eigenvalues are a
%            complex conjugate pair.
%   BETA     the m-vector of the scalars beta that define the
%            eigenvalues of the pencil aS - bT.
%            Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
%            beta = BETA(j) represent the j-th eigenvalue of the pencil
%            aS - bT, in the form lambda = alpha/beta. Since lambda may
%            overflow, the ratios should not, in general, be computed.
%            Due to the skew-Hamiltonian/skew-Hamiltonian structure of
%            the pencil, every eigenvalue occurs twice and thus it has
%            only to be saved once in ALPHAR, ALPHAI and BETA.
%   Ao     - if job = 1, the computed m-by-m submatrix Ao in (2).
%   Do     - if job = 1, the computed m-by-m strictly upper triangular
%            part of the skew-symmetric matrix Do in (2). The lower
%            triangle is not set.
%   Bo     - if job = 1, the computed m-by-m submatrix Bo in (2).
%   Fo     - if job = 1, the computed m-by-m strictly upper triangular
%            part of the skew-symmetric matrix Fo in (2). The lower
%            triangle is not set.
%   Q      - if compq > 0, an n-by-n matrix containing the computed
%            orthogonal matrix Q (if compq = 1) or Q0*Q (if compq = 2).
%            If job = 0, Do and Fo contain the corresponding matrices
%            just before the application of the QZ algorithm, but Ao and
%            Bo contain meaningless elements.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2012.
%
% Revisions:
%   -
%
