% SYMPLURV.F - Gateway function for computing the eigenvalues of a real
%              skew-Hamiltonian/Hamiltonian pencil in factored form,
%              using SLICOT routine MB04AD.
%
%   [ALPHAR,ALPHAI,BETA(,To,Zo,Ho(,Q1,Q2,U11,U12,U21,U22))] =
%           symplURV(Z,H(,job,compq1(,Q01),compq2(,Q02),
%                             compu1(,U011,U012),compu2(,U021,U022)))
%
%   SYMPLURV computes the eigenvalues of a real skew-Hamiltonian/Hamiltonian
%   pencil aS - bH with
%
%                                    (  0  I  )
%     S = T Z = J Z' J' Z, where J = (        ),                   (1)
%                                    ( -I  0  )
%
%   via generalized symplectic URV decomposition. That is, orthogonal
%   matrices Q1 and Q2 and orthogonal symplectic matrices U1 and U2
%   are computed such that
%
%                                 (  T11  T12 )
%     Q1' T U1 = Q1' J Z' J' U1 = (           ) = To,
%                                 (   0   T22 )
%
%                (  Z11  Z12 )
%     U2' Z Q2 = (           ) = Zo,                             (2)
%                (   0   Z22 )
%
%                ( H11  H12 )
%     Q1' H Q2 = (          ) = Ho,
%                (  0   H22 )
%
%   where T11, T22', Z11, Z22', H11 are upper triangular and H22' is
%   upper quasi-triangular.
%   Optionally, if compq1 = 1 (compq2 = 1), the orthogonal matrix Q1
%   (Q2) that fulfills (2) is computed.
%   Optionally, if compu1 = 1 (compu2 = 1), the orthogonal symplectic
%   matrix 
%   
%          (  U11  U12  )          (  U21  U22  )
%     U1 = (            )   ( U2 = (            ) )
%          ( -U12  U11  )          ( -U22  U21  )
%
%   that fulfills (2) is computed. Only U11 and U12 (U21 and U22 ) are
%   returned.
%
% Description of input parameters:
%   Z      - the n-by-n matrix Z, n even.
%   H      - the n-by-n Hamiltonian matrix H.
%   job    - (optional) scalar indicating the computation to be
%            performed, as follows:
%            = 0 :  compute the eigenvalues only (default);
%            = 1 :  compute the eigenvalues and the matrices of the
%                   transformed pencil in (2).
%   compq1 - (optional) scalar indicating whether the orthogonal
%            transformation matrix Q1 is returned, or if Q1 is not
%            required, as follows:
%            = 0 :  Q1 is not required (default);
%            = 1 :  on exit, Q1 contains the orthogonal matrix Q1;
%            = 2 :  the orthogonal transformations are accumulated
%                   into Q1;
%                   on input, Q01 must contain an orthogonal matrix Q01;
%                   on exit, Q1 contains Q01*Q1, with Q1 returned for
%                   compq1 = 1.
%   Q01    - if compq1 = 2, the n-by-n orthogonal matrix Q01.
%   compq2 - (optional) scalar indicating whether the orthogonal
%            transformation matrix Q2 is returned, or if Q2 is not
%            required, as follows:
%            = 0 :  Q2 is not required (default);
%            = 1 :  on exit, Q2 contains the orthogonal matrix Q2;
%            = 2 :  the orthogonal transformations are accumulated
%                   into Q2;
%                   on input, Q02 must contain an orthogonal matrix Q02;
%                   on exit, Q2 contains Q02*Q2, with Q2 returned for
%                   compq2 = 1.
%   Q02    - if compq2 = 2, the n-by-n orthogonal matrix Q02.
%   compu1 - (optional) scalar indicating whether the orthogonal
%            symplectic transformation matrix U1 is returned, or if U1
%            is not required, as follows:
%            = 0 :  U1 is not required (default);
%            = 1 :  on exit, U11 and U12 contain the corresponding
%                   submatrices of the orthogonal symplectic matrix U1;
%            = 2 :  the orthogonal transformations are accumulated
%                   into U1;
%                   on input, U011, U012 must contain the corresponding
%                   submatrices of an orthogonal symplectic matrix U01;
%                   on exit, U11 and U12 contain the updated submatrices
%                   U11 and U12 of the matrix product U01*U1, with U1
%                   returned for compu1 = 1.
%   U011,  - if compu1 = 2, the m-by-m submatrices U011 and U012 of U01,
%   U012,    respectively.
%   compu2 - (optional) scalar indicating whether the orthogonal
%            symplectic transformation matrix U2 is returned, or if U2
%            is not required, as follows:
%            = 0 :  U2 is not required (default);
%            = 1 :  on exit, U21 and U22 contain the corresponding
%                   submatrices of the orthogonal symplectic matrix U2;
%            = 2 :  the orthogonal transformations are accumulated
%                   into U2;
%                   on input, U021, U022 must contain the corresponding
%                   submatrices of an orthogonal symplectic matrix U02;
%                   on exit, U21 and U22 contain the updated submatrices
%                   U21 and U22 of the matrix product U02*U2, with U2
%                   returned for compu2 = 1.
%   U021,  - if compu2 = 2, the m-by-m submatrices U021 and U022 of U02,
%   U022,    respectively.
%
% Description of output parameters:
%   ALPHAR,- the m-vectors of real parts and imaginary parts,
%   ALPHAI   respectively, of each scalar alpha defining an eigenvalue
%            of the pencil aS - bH, with m = n/2.
%            If ALPHAI(j) is zero, then the j-th eigenvalue is real.
%   BETA     the m-vector (or (n+1)-vector) of the scalars beta (and
%            possibly the scaling factors) that define the eigenvalues
%            of the pencil aS - bH. The length of beta is decided by the
%            MEX-file and it should be checked. Usually, the quantities
%            alpha = (ALPHAR(j),ALPHAI(j)), and beta = BETA(j) represent
%            together the j-th eigenvalue of the pencil aS - bH, in the
%            form lambda = alpha/beta. Since lambda may overflow, the
%            ratios should not, in general, be computed. Due to the
%            skew-Hamiltonian/Hamiltonian structure of the pencil, only
%            half of the spectrum is saved in ALPHAR, ALPHAI and BETA.
%            Specifically, only eigenvalues with imaginary parts greater
%            than or equal to zero are stored; their conjugate
%            eigenvalues are not stored. If imaginary parts are zero
%            (i.e., for real eigenvalues), only positive eigenvalues
%            are stored.
%            If one or more BETA(j) is not representable, BETA(m+j),
%            ..., BETA(n) return the scaling parameters for the
%            eigenvalues of the pencil aS - bH. Specifically, the j-th
%            eigenvalue is represented by the quantities alpha,
%            beta, and gamma = BETA(m+j) in the form
%            lambda = (alpha/beta) * b**gamma, where b is the
%            machine base, stored in BETA(n+1).
%   To     - the computed n-by-n matrix To in (2).
%   Zo     - the computed n-by-n matrix Zo in (2). Z(m+1:n,1:m) is not
%            set to 0.
%   Ho     - the computed n-by-n matrix Ho in (2).
%   Q1     - if compq1 > 0, an n-by-n matrix containing the computed
%            or updated orthogonal matrix Q1.
%   Q2     - if compq2 > 0, an n-by-n matrix containing the computed
%            or updated orthogonal matrix Q2.
%   U11,   - if compu1 > 0, the m-by-m matrices containing the computed
%   U12      or updated submatrices U11 and U12 of the orthogonal
%            symplectic matrix U1.
%   U21,   - if compu2 > 0, the m-by-m matrices containing the computed
%   U22      or updated submatrices U21 and U22 of the orthogonal
%            symplectic matrix U2.
%            If job = 0, To, Zo, Ho, Q1, Q2, U11, U12, U21, and U22
%            contain the corresponding matrices just before the
%            application of the periodic QZ algorithm.
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
