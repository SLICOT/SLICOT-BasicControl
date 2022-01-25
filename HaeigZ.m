% HAEIGZ.F - Gateway function for computing a balancing transformation
%            or the eigenvalues of a complex Hamiltonian matrix, using
%            SLICOT Library routines MB04DZ and MB03XZ.
%            The gateway accepts real or complex input matrices.
%
% Matlab call:
%   [W(,Ao,Go(,U1,U2(,l,scal)))]   = HaeigZ(A,QG(,job,jobu,balanc))
%   [W(,Se,De(,Ue1,Ue2(,l,scal)))] = HaeigZ(A,QG(,job,jobu,balanc))
%   [l,scal(,Ab,QGb)]              = HaeigZ(A,QG,job,balanc)
%
%   [l,scal]                 = HaeigZ(A,QG,-1,balanc)
%   [l,scal,Ab,QGb]          = HaeigZ(A,QG,-1,balanc)
%   [W]                      = HaeigZ(A,QG)
%   [W,Se,De]                = HaeigZ(A,QG)
%   [W,Se,De,Ue1,Ue2]        = HaeigZ(A,QG,0,1)
%   [W,Se,De,Ue1,Ue2,l,scal] = HaeigZ(A,QG,0,1,balanc)
%   [W,Ao]                   = HaeigZ(A,QG,1)
%   [W,Ao,U1,U2]             = HaeigZ(A,QG,1,1)
%   [W,Ao,Go]                = HaeigZ(A,QG,2)
%   [W,Ao,Go,U1,U2]          = HaeigZ(A,QG,2,1)
%   [W,Ao,Go,U1,U2,l,scal]   = HaeigZ(A,QG,2,1,balanc)
%
%   HAEIGZ computes the eigenvalues of a complex n-by-n Hamiltonian matrix H,
%   with
%
%                 [  A   G  ]         H        H
%           H  =  [       H ],   G = G ,  Q = Q ,                 (1)
%                 [  Q  -A  ]
%
%   where A, G and Q are complex n-by-n matrices.
%
%   Due to the structure of H, if lambda is an eigenvalue, then
%   -conjugate(lambda) is also an eigenvalue. This does not mean that
%   purely imaginary eigenvalues are necessarily multiple. The function
%   computes the eigenvalues of H using an embedding to a skew-
%   Hamiltonian matrix He,
%
%                  [  Ae   Ge  ]            T            T
%           He  =  [         T ],   Ge = -Ge ,   Qe = -Qe ,        (2)
%                  [  Qe   Ae  ]
%
%   where Ae, Ge, and Qe are real 2*n-by-2*n matrices. Then, an
%   orthogonal symplectic matrix Ue is used to reduce He to the
%   structured real Schur form,
%
%         T          [  Se   De ]            T
%        Ue He Ue =  [        T ],   De = -De ,                    (3)
%                    [  0    Se ]
%
%   where Ue is a 4n-by-4n real symplectic matrix, and Se is upper
%   quasi-triangular (real Schur form).
%
%   Optionally, if job > 0, the matrix i*He is further transformed to
%   the structured complex Schur form
%
%         H            [  Ao  Go ]          H
%        U (i*He) U =  [       H ],  Go = Go ,                     (4)
%                      [  0  -Ao ]
%
%   where U is a 4n-by-4n unitary symplectic matrix, and Ao is upper
%   triangular (Schur form). Optionally, if jobu = 1, the unitary
%   symplectic transformation matrix
%
%         (  U1  U2  )
%     U = (          )
%         ( -U2  U1  )
%
%   is computed.
%
%   If job = -1, an accurate similarity transformation T such that
%   Hb = T\H*T has, as nearly as possible, approximately equal row and
%   column norms. T is a permutation of a diagonal matrix and
%   symplectic. T is stored in an n-vector scal as described in MB04DZ.
%
% Description of input parameters:
%   A      - the n-by-n matrix A.
%   QG     - an  n-by-(n+1) matrix containing the triangles of the
%            Hermitian matrices G and Q, as follows:
%            the leading n-by-n lower triangular part contains the
%            lower triangle of the matrix Q, and the n-by-n upper
%            triangular part of the submatrix in the columns 2 to n+1
%            contains the upper triangle of the matrix G of H in (1).
%            So, if i >= j, then Q(i,j) = conj(Q(j,i)) is stored in
%            QG(i,j) and G(j,i) = conj(G(i,j)) is stored in QG(j,i+1).
%            QG is an empty matrix if n = 0.
%   job    - (optional) scalar indicating the computation to be
%            performed, as follows:
%            = -1 :  compute a balancing transformation only;
%            =  0 :  compute the eigenvalues only (default);
%            =  1 :  compute the eigenvalues and the matrix Ao in (4);
%            =  2 :  compute the eigenvalues and the matrices Ao and Go
%                    in (4).
%   jobu   - (optional) if job > 0, scalar indicating whether the
%            unitary transformation matrix U is returned, as follows:
%            = 0 :  U is not required (default);
%            = 1 :  on exit, U contains the unitary transformation
%                   matrix.
%   balanc - determines whether H should be permuted (balanc = 1),
%            scaled (balanc = 2), or permuted and scaled (balanc = 3)
%            prior to eigenvalue computations. Otherwise balanc = 0
%            (default). This parameter is optional if job >= 0, but
%            compulsory if job < 0.
%
% Description of output parameters:
%   l      - if job = -1 or balanc > 0, an integer determined when H was
%            balanced. The balanced A(I,J) = 0 if I > J and J = 1:l-1.
%            The balanced Q(I,J) = 0 if J = 1:l-1 or I = 1:l-1.
%   scal   - if job = -1 or balanc > 0, an n-vector containing details
%            of the permutation and/or scaling factors applied when
%            balancing. See MB04DZ, for details.
%   Ab     - if job = -1, the matrix A of the balanced Hamiltonian Hb.
%            The lower triangular part of the first l-1 columns of A is
%            zero.
%   QGb    - if job = -1, the matrices Q and G of the balanced
%            Hamiltonian Hb, stored compactly. The lower triangular and
%            diagonal part of the first l-1 columns of QGb is zero.
%   W      - the 2n-vector of the eigenvalues of the matrix H.
%   Se     - if job = 0, the computed 2n-by-2n upper real Schur
%            submatrix Se in (3).
%   De     - if job = 0, the computed 2n-by-2n skew-symmetric submatrix
%            De in (3).
%   Ue1    - if job = 0, the computed 2n-by-2n (1,1) block of the matrix
%            Ue in (3).
%   Ue2    - if job = 0, the computed 2n-by-2n (2,1) block of the matrix
%            Ue in (3).
%   Ao     - if job > 0, the computed 2n-by-2n upper triangular
%            submatrix Ao in (4).
%   Go     - if job = 2, the computed 2n-by-2n Hermitian matrix Go
%            in (4).
%   U1     - if jobu = 1, a 2n-by-2n matrix containing the computed
%            matrix U1.
%   U2     - if jobu = 1, an 2n-by-2n matrix containing the computed
%            matrix U2.
%
% Contributors:
%   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2011.
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
