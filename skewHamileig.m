% SKEWHAMILEIG.F - Gateway function for computing the eigenvalues of a
%                  skew-Hamiltonian/Hamiltonian pencil, using SLICOT
%                  routine MB04BD.
%
%   [ALPHAR,ALPHAI,BETA(,Ao,Do,Bo,Fo,C1o,Vo,C2o(,Q1,Q2))] =
%                      skewHamileig(A,DE,C,VW(,job,compq1(,Q),compq2))
%
%   [ALPHAR,ALPHAI,BETA]                 = skewHamileig(A,DE,C,VW)
%   [ALPHAR,ALPHAI,BETA,Ao,Do,Bo,Fo,C1o,Vo,C2o] ...
%                                        = skewHamileig(A,DE,C,VW,1)
%
%   SKEWHAMILEIG computes the eigenvalues of a real n-by-n skew-Hamiltonian/
%   Hamiltonian pencil aS - bH with n = 2m,
%
%         (  A  D  )         (  C  V  )
%     S = (        ) and H = (        ).                             (1)
%         (  E  A' )         (  W -C' )
%
%   Optionally, if job = 1, decompositions of S and H will be
%   computed via orthogonal transformations Q1 and Q2 as follows:
%
%                     (  Ao  Do  )
%     Q1' S J Q1 J' = (          ),
%                     (   0  Ao' )
%
%                     (  Bo  Fo  )
%     J' Q2' J S Q2 = (          ) =: T,                             (2)
%                     (   0  Bo' )
%
%                (  C1o  Vo   )            (  0  I  )
%     Q1' H Q2 = (            ), where J = (        ),
%                (  0    C2o' )            ( -I  0  )
%
%   and Ao, Bo, C1o are upper triangular, C2o is upper quasi-triangular
%   and Do and Fo are skew-symmetric. The notation M' denotes the
%   transpose of the matrix M.
%   Optionally, if compq1 = 1, the orthogonal transformation matrix
%   Q1 will be computed.
%   Optionally, if compq2 = 1, the orthogonal transformation matrix
%   Q2 will be computed.
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
%   C      - the m-by-m matrix C.
%   VW     - an  m-by-(m+1) matrix containing the triangles of the
%            symmetric matrices V and W, as follows:
%            the leading m-by-m lower triangular part contains the lower
%            triangle of the matrix W, and the m-by-m upper triangular
%            part of the submatrix in the columns 2 to m+1 contains the
%            upper triangle of the matrix V of H in (1).
%            So, if i >= j, then W(i,j) = W(j,i) is stored in VW(i,j)
%            and V(j,i) = V(i,j) is stored in VW(j,i+1).
%            VW is an empty matrix if m = 0.
%   job    - (optional) scalar indicating the computation to be
%            performed, as follows:
%            = 0 :  compute the eigenvalues only (default);
%            = 1 :  compute the eigenvalues and the matrices of the
%                   transformed pencil in (2).
%   compq1 - (optional) scalar indicating whether the orthogonal
%            transformation matrix Q1 is returned or accumulated
%            into an orthogonal matrix, or if Q1 is not required,
%            as follows:
%            = 0 :  Q1 is not required (default);
%            = 1 :  on entry, Q1 need not be set;
%                   on exit, Q1 contains the orthogonal matrix Q1;
%            = 2 :  the orthogonal transformations are accumulated
%                   into Q1;
%                   on input, Q1 must contain an orthogonal matrix Q;
%                   on exit, Q1 contains Q*Q1.
%   Q      - if compq1 = 2, the n-by-n orthogonal matrix Q.
%   compq2 - (optional) scalar indicating whether the orthogonal
%            transformation matrix Q2 is returned or accumulated
%            into an orthogonal matrix, or if Q2 is not required,
%            as follows:
%            = 0 :  Q2 is not required (default);
%            = 1 :  on exit, Q2 contains the orthogonal matrix Q2;
%            = 2 :  on exit, Q2 contains the matrix product J*Q*J'*Q2,
%                   where Q2 is the product of the orthogonal
%                   transformations that are applied to the pencil
%                   aS - bH to reduce S and H to the forms in (2),
%                   for compq2 = 1.
%                   If nonzero, compq2 should coincide with compq1.
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
%            eigenvalue, and thus it has only to be saved once in
%            ALPHAR, ALPHAI and BETA.
%            Specifically, only eigenvalues with imaginary parts greater
%            than or equal to zero are stored; their conjugate 
%            eigenvalues are not stored. If imaginary parts are zero
%            (i.e., for real eigenvalues), only positive eigenvalues
%            are stored.
%   Ao     - if job = 1, the computed m-by-m submatrix Ao in (2).
%   Do     - if job = 1, the computed m-by-m strictly upper triangular
%            part of the skew-symmetric matrix Do in (2).
%   Bo     - if job = 1, the computed m-by-m submatrix Bo in (2).
%   Fo     - if job = 1, the computed m-by-m strictly upper triangular
%            part of the skew-symmetric matrix Fo in (2).
%   C1o    - if job = 1, the computed m-by-m submatrix C1o in (2).
%   Vo     - if job = 1, the computed m-by-m matrix Vo in (2).
%   C2o    - if job = 1, the computed m-by-m submatrix C2o in (2).
%   Q1     - if compq1 > 0, an n-by-n matrix containing the computed
%            orthogonal matrix Q1 (if compq1 = 1) or Q*Q1
%            (if compq1 = 2).
%   Q2     - if compq2 > 0, an n-by-n matrix containing the computed
%            orthogonal matrix Q2 (if compq2 = compq1 = 1) or J*Q*J'*Q2
%            (if compq2 = compq1 = 2).
%            If job = 0, Ao, Do, Bo, Fo, C1o, Vo, and C2o contain the
%            corresponding matrices just before the application of the
%            periodic QZ algorithm.
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
