% SKEWHAMIL2FEIG.F - Gateway function for MB04ED.
%
%   [ALPHAR,ALPHAI,BETA(,Zo,Bo,Fo(,Q,U1,U2))] =
%           skewHamil2feig(Z,B,FG(,job,compq,compu(,U01,U02)))
%
%   [ALPHAR,ALPHAI,BETA]          = skewHamil2feig(Z,B,FG)
%   [ALPHAR,ALPHAI,BETA,Zo,Bo,Fo] = skewHamil2feig(Z,B,FG,1)
%
%   SKEWHAMIL2FEIG computes the eigenvalues of a real n-by-n 
%   skew-Hamiltonian/skew-Hamiltonian pencil aS - bT, with n = 2m, and
%
%                           (  B  F  )            (  0  I  )
%     S = J Z' J' Z and T = (        ), where J = (        ).        (1)
%                           (  G  B' )            ( -I  0  )
%
%   Optionally, if job = 1, aS - bT will be transformed to
%   structured Schur form: an orthogonal transformation matrix Q and
%   an orthogonal symplectic transformation matrix U are computed,
%   such that
%
%              (  Z11  Z12  )
%     U' Z Q = (            ) = Zo, and
%              (   0   Z22' )
%                                                                    (2)
%                   (  Bo  Fo  )
%     J Q' J' T Q = (          ),
%                   (   0  Bo' )
%
%   where Z11 and Z22 are upper triangular and Bo is upper quasi-
%   triangular. The notation M' denotes the transpose of the matrix M.
%   Optionally, if compq <> 0, the orthogonal transformation matrix Q
%   will be computed.
%   Optionally, if compu > 0, the orthogonal symplectic transformation
%   matrix
%
%         (  U1  U2  )
%     U = (          )
%         ( -U2  U1  )
%
%   will be computed.
%
% Description of input parameters:
%   Z      - the n-by-n matrix Z.
%   B      - the m-by-m matrix B.
%   FG     - an  m-by-(m+1) matrix containing the strict triangles of
%            the skew-symmetric matrices F and G, as follows:
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
%            transformation matrix Q is returned, as follows:
%            = 0 :  Q is not required (default);
%            = 1 :  on exit, Q contains the orthogonal matrix Q;
%            =-1 :  on exit, for job = 0, Q contains the orthogonal
%                   matrix Q1 which reduced Z to 
%
%                            (  Z11  Z12  )
%                     Z*Q1 = (            ),
%                            (   0   Z22  )
%
%                   where Z11 and Z22' are upper triangular (the first
%                   step of the algorithm).
%   compu  - (optional) scalar indicating whether the orthogonal
%            symplectic transformation matrix U is returned, as follows:
%            = 0 :  U is not required (default);
%            = 1 :  on exit, U1 and U2 contain the submatrices 
%                   of the orthogonal symplectic matrix U;
%            = 2 :  the orthogonal transformations are accumulated
%                   into U;
%                   on input, U01 and U02 must contain the corresponding
%                   submatrices of an orthogonal symplectic matrix U0;
%                   on exit, U1 and U2 contain the updated submatrices
%                   U1 and U2 of the matrix product U0*U, with U
%                   returned for compu = 1.
%   U01,   - if compu = 2, the m-by-m submatrices U01 and U02 of U0,
%   U02      respectively.
%
% Description of output parameters:
%   ALPHAR,- the m-vectors of real parts and imaginary parts,
%   ALPHAI   respectively, of each scalar alpha defining an eigenvalue
%            of the pencil aS - bT.
%            If ALPHAI(j) is zero, then the j-th eigenvalue is real.
%   BETA     the m-vector of the scalars beta that define the eigenvalues
%            of the pencil aS - bT.
%            Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
%            beta = BETA(j) represent the j-th eigenvalue of the pencil
%            aS - bT, in the form lambda = alpha/beta. Since lambda may
%            overflow, the ratios should not, in general, be computed.
%            Due to the skew-Hamiltonian/skew-Hamiltonian structure of
%            the pencil, every eigenvalue occurs twice and thus it has
%            only to be saved once in ALPHAR, ALPHAI and BETA.
%   Zo     - if job = 1, the computed n-by-n matrix Zo in (2).
%   Bo     - if job = 1, the computed m-by-m strictly upper triangular
%            part of the skew-symmetric matrix Bo in (2).
%   Fo     - if job = 1, the computed m-by-m strictly upper triangular
%            part of the skew-symmetric matrix Fo in (2).
%   Q      - if compq <> 0, an n-by-n matrix containing the computed
%            orthogonal matrix Q or Q1.
%   U1     - if compu > 0, an m-by-m matrix containing the computed
%            matrix U1.
%   U2     - if compu > 0, an m-by-m matrix containing the computed
%            matrix U2.
%            If job = 0, Zo, Bo, Fo, Q, U1, and U2 contain the
%            corresponding matrices just before the application of the
%            periodic QZ algorithm.
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
