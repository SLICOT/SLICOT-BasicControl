% HAPACK_HAEIG.F - Gateway function for Hamiltonian eigenvalue computations.
%
%                    [...] = hapack_haeig(task,...)
%
%   PVL decomposition
%   task = 1 :    [Ar,QGr] = hapack_haeig(1,A,QG)
%           [U1,U2,Ar,QGr] = hapack_haeig(1,A,QG)
%
%   URV decomposition
%   task = 2 :     [S,T,G] = hapack_haeig(2,A,QG)
%            [U1,U2,S,T,G] = hapack_haeig(2,A,QG)
%      [U1,U2,V1,V2,S,T,G] = hapack_haeig(2,A,QG)
%                     [Br] = hapack_haeig(2,B)
%               [U1,U2,Br] = hapack_haeig(2,B)
%         [U1,U2,V1,V2,Br] = hapack_haeig(2,B)
%
%   URV/PS decomposition
%   task = 3 : [S,T,wr,wi] = hapack_haeig(3,A,QG)
%            [S,T,G,wr,wi] = hapack_haeig(3,A,QG)
%        [U1,U2,S,T,wr,wi] = hapack_haeig(3,A,QG)
%      [U1,U2,S,T,G,wr,wi] = hapack_haeig(3,A,QG)
%        [U1,U2,V1,V2,...
%               S,T,wr,wi] = hapack_haeig(3,A,QG)
%       [U1,U2,V1,V2,...
%             S,T,G,wr,wi] = hapack_haeig(3,A,QG)
%
%   Balancing
%   task = 4 : [scale,ilo] = hapack_haeig(4,A,QG,balanc)
%       [Ar,QGr,scale,ilo] = hapack_haeig(4,A,QG,balanc)
%
%   Simply eigenvalues:
%   task = 5:      [wr,wi] = hapack_haeig(5,A,QG,balanc)
%
%   Reordering Hamiltonian Schur form:
%   task = 6 :    [Ar,QGr] = hapack_haeig(6,A,QG,select,lower)
%           [U1,U2,Ar,QGr] = hapack_haeig(6,A,QG,select,lower)
%           [U1,U2,Ar,QGr] = hapack_haeig(6,U1,U2,A,QG,select,lower)
%
%   Selected stable/unstable invariant subspaces
%   task = 7 :        [US] = hapack_haeig(7,S,T,U1,U2,V1,V2,select)
%                  [US,UU] = hapack_haeig(7,S,T,U1,U2,V1,V2,select)
%
%   task = 8 :        [UU] = hapack_haeig(8,S,T,U1,U2,V1,V2,select)
%
%   Complete stable/unstable invariant subspaces:
%   task = 9 :        [US] = hapack_haeig(9,S,T,U1,U2,V1,V2)
%                     [US] = hapack_haeig(9,S,T,G,U1,U2,V1,V2)
%                  [US,UU] = hapack_haeig(9,S,T,U1,U2,V1,V2)
%                  [US,UU] = hapack_haeig(9,S,T,G,U1,U2,V1,V2)
%
%   task = 10 :       [UU] = hapack_haeig(10,S,T,U1,U2,V1,V2)
%                     [UU] = hapack_haeig(10,S,T,G,U1,U2,V1,V2)
% Purpose:
%   To compute eigenvalues and invariant subspaces of a 2*N-by-2*N
%   Hamiltoninan matrix H = [A, G; Q, -A']. On input as well as on
%   output Hamiltonian matrices are stored in compressed format in
%   an N-by-N array Ar and an N-by-(N+1) array QGr.
%
%   task =  1: Computes the PVL form Hr = [Ar, Gr; Qr, -Ar'] of a
%              Hamiltonian H with diagonal Qr and Ar in upper Hessenberg
%              form. Optionally, N-by-N matrices U1 and U2 are computed
%              so that U = [U1, U2; -U2, U1] satisfies U'*U = I and
%              H = U*Hr*U'.
%
%   task =  2: Computes the symplectic URV form Br = [T, G; 0, S'] of a
%              Hamiltonian H or a general matrix B so that T is upper
%              triangular and S is in upper Hessenberg form. Optionally,
%              N-by-N matrices U1, U2, V1 and V2 are computed so that
%              U = [U1, U2; -U2, U1], V = [V1', V2'; -V2', V1'] (or
%              V = [V1, V2'; -V2', V1] if B is Hamiltonian and stored in
%              compressed format) satisfy U'*U = I, V'*V = I and
%              B = U*Br*V'.
%
%   task =  3: Computes the symplectic URV/periodic Schur form
%              Hr = [T, G; 0, S'] of a Hamiltonian matrix H so that T is
%              upper triangular and S is in real Schur form. Optionally,
%              N-by-N matrices U1, U2, V1 and V2 are computed so that
%              U = [U1, U2; -U2, U1], V = [V1, V2; -V2, V1] satisfy
%              U'*U = I, V'*V = I and H = U*Hr*V'. Also, the nonpositive
%              square roots of the eigenvalues of -S*T, which are equal
%              to the stable eigenvalues of H, are computed and returned
%              in vectors wr and wi containing the real and imaginary
%              parts, respectively.
%
%   task =  4: Finds an accurate similarity transformation T such that
%              Hr = T\H*T has, as nearly as possible, approximately
%              equal row and column norms. T is a permutation of a
%              diagonal matrix and symplectic. T is stored in an
%              n-vector SCALE as described in MB04DD.
%
%   task =  5: Computes the stable eigenvalues of H and returns vectors
%              wr and wi containing the real and imaginary parts,
%              respectively.
%
%   task =  6: Reorders the eigenvalues selected by the arrays SELECT
%              and LOWER to the top left part of a matrix H in
%              Hamiltonian Schur form. Optionally, the corresponding
%              orthogonal symplectic transformation matrix
%              U = [U1, U2; -U2, U1] is computed or a given matrix
%              U = [U1, U2; -U2, U1] is post-multiplied by this matrix.
%
%   task =  7
%   task =  8: Computes the invariant subspace US (UU) belonging to
%              a set of selected stable (unstable) eigenvalues of H.
%              The user must provide the matrices U1, U2, V1, V2, S
%              and T as computed by MB03XD (task 3).
%
%   task =  9
%   task = 10: Computes the invariant subspace US (UU) belonging to
%              all stable (unstable) eigenvalues of H. The user must
%              provide the matrices U1, U2, V1, V2, S, T and
%              possibly G as computed by MB03XD (task 3). Note that
%              if G is provided then US and UU are computed from a
%              larger set of vectors which may imply more accurate
%              bases.
%
% Description of input parameters:
%   task   - integer option to determine the computation to perform as
%            described above.
%   A      - real N-by-N matrix. If task = 6 this matrix must be in real
%            Schur form.
%   B      - real 2*N-by-2*N matrix.
%   G      - real N-by-N matrix.
%   QG     - real N-by-(N+1) matrix containing the symmetric matrices
%            Q and G in compressed format. If task = 6 the part
%            containing the matrix Q must be zero.
%   S      - a real N-by-N matrix in real Schur form.
%   T      - a real, upper triangular N-by-N matrix.
%   U1,U2  - real N-by-N matrices.
%   V1,V2  - real N-by-N matrices.
%   balanc - determines whether H should be permuted (balanc = 'p'),
%            scaled (balanc = 's') or permuted and scaled (balanc = 'b')
%            prior to eigenvalue computations. Otherwise balanc = 'n'.
%   lower  - N-vector controlling whether lambda or -lambda of each
%            selected eigenvalue pair (lambda,-lambda) is to be
%            reordered to the top part of a Hamiltonian Schur form.
%   select - if task = 6: N-vector specifying the eigenvalues to be
%            reordered to the top part of a Hamiltonian Schur form; and
%            if task > 6: N-vector specifying the eigenvalue pairs to
%            which stable/unstable invariant subspaces are to be
%            computed.
%
% Description of output parameters:
%   Ar     - real N-by-N matrix containing the (1,1) block of the
%            returned Hamiltonian matrix Hr.
%   QGr    - real N-by-(N+1) matrix containing the symmetric (2,1) and
%            (1,2) blocks of the returned (skew-)Hamiltonian matrix Hr
%            in compressed format.
%   Br     - real 2*N-by-2*N matrix in symplectic URV form.
%   T,S,G  - submatrices of a symplectic URV(/periodic Schur) form of H.
%   U1,U2  - real N-by-N matrices so that U = [U1, U2; -U2, U1]
%            satisfies U'*U = I.
%   V1,V2  - real N-by-N matrices so that V = [V1', V2'; -V2', V1'] (if
%            task = 2) or V = [V1, V2; -V2, V1] (otherwise) satisfies
%            V'*V = I.
%   US     - real 2*N-by-M matrix containing the orthonormal basis of a
%            stable invariant subspace.
%   UU     - real 2*N-by-M matrix containing the orthonormal basis of an
%            unstable invariant subspace.
%   ilo    - ilo-1 is the number of deflated eigenvalues in the balanced
%            Hamiltonian matrix.
%   scale  - real N-vector containing the scaling factors as returned by
%            MB04DD.
%   wr     - real N-vector containing the real parts of the eigenvalues
%            of H.
%   wi     - real N-vector containing the imaginary parts of the
%            eigenvalues of H.

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   D. Kressner, UCL CESAME, Louvain-la-Neuve, June 2003.
%
% Revisions:
%   V. Sima, Nov. 2012.
%
