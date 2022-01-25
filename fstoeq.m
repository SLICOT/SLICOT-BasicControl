%FSTOEQ  MEX-function for computing orthogonal-triangular decompositions
%        of (block) Toeplitz matrices and/or solving associated linear
%        least-squares systems using SLICOT routines MB02HD, MB02ID,
%        MB02JD, MB02JX, MB02KD.
%
%            [(Q,)R(,X,Y)] = fstoeq(task,TC,TR(,B,C))
%                  [(Q,)R] = fstoeq(task,TC,TR,tol)
%                [R(,X,Y)] = fstoeq(task,m,n,TCB,TRB(,B,C))
%                  [X(,Y)] = fstoeq(task,TC,TR,B(,C))
%                      [Y] = fstoeq(task,TC,TR,C)
%
%   task =  1 :        [R] = fstoeq( 1,TC,TR)
%                    [Q,R] = fstoeq( 1,TC,TR)
%   task =  2 :      [R,X] = fstoeq( 2,TC,TR,B)
%                  [Q,R,X] = fstoeq( 2,TC,TR,B)
%   task =  3 :      [R,Y] = fstoeq( 3,TC,TR,C)
%                  [Q,R,Y] = fstoeq( 3,TC,TR,C)
%   task =  4 :    [R,X,Y] = fstoeq( 4,TC,TR,B,C)
%                [Q,R,X,Y] = fstoeq( 4,TC,TR,B,C)
%   task =  5 :      [R,E] = fstoeq( 5,TC,TR,tol)
%                  [Q,R,E] = fstoeq( 5,TC,TR,tol)
%   task =  6 :        [R] = fstoeq( 6,m,n,TCB,TRB)
%   task =  7 :      [R,X] = fstoeq( 7,m,n,TCB,TRB,B)
%   task =  8 :      [R,Y] = fstoeq( 8,m,n,TCB,TRB,C)
%   task =  9 :    [R,X,Y] = fstoeq( 9,m,n,TCB,TRB,B,C)
%   task = 10 :        [X] = fstoeq(10,TC,TR,B)
%   task = 11 :        [Y] = fstoeq(11,TC,TR,C)
%   task = 12 :      [X,Y] = fstoeq(12,TC,TR,B,C)
%   task = 13 :        [X] = fstoeq(13,TC,TR,B)
%
% Purpose:
%   To compute an orthogonal-triangular/trapezoidal decomposition
%   (with partial column pivoting) of a (banded) block Toeplitz matrix.
%   If task <= 5 or task >= 10, the first block column of T is contained
%   in TC and has the dimension M*K-by-L and the first block row of T is
%   contained in TR and has the dimension K-by-N*L, while for task >= 6
%   and task <= 9 the leading nonzero blocks of the first block column
%   of T are contained in TCB which has the dimension (ML+1)*K-by-L
%   and the leading nonzero blocks of the first block row of T are
%   contained in TRB which has the dimesion K-by-(NU+1)*L.
%   If task = 13, the product of a block Toeplitz matrix T with a block
%   column vector C is computed.
%
%   task =  1: Compute R, the Cholesky factor of T'*T, i.e., R is lower
%              triangular, so that R*R' = T'*T. Optionally, a matrix Q
%              is computed so that Q'*Q = I and T = Q*R'.
%
%   task =  2: Compute R and the least-squares solution of T*X = B for a
%              given right hand side matrix B.
%
%   task =  3: Compute R and the minimum norm solution of T'*Y = C for a
%              given right hand side matrix C.
%
%   task =  4: Compute R, X and Y.
%
%   task =  5: Compute R, a low rank Cholesky factor of T'*T, i.e.,
%              R is lower triangular, so that R*R' = (T*E)'*(T*E)
%              for a permutation matrix E. Optionally, a matrix
%              Q is computed so that Q'*Q = I and T*E = Q*R'.
%
%   task =  6: Compute R, the Cholesky factor of T'*T, where T is a
%              banded block Toeplitz matrix.
%
%   task =  7: Compute R and the least-squares solution of T*X = B for a
%              banded block Toeplitz matrix T and given right hand side
%              matrix B.
%
%   task =  8: Compute R and the minimum norm solution of T'*Y = C for a
%              banded block Toeplitz matrix and given right hand side
%              matrix C.
%
%   task =  9: Compute R, X and Y for a banded block Toeplitz matrix.
%
%   task = 10: Compute the least-squares solution of T*X = B for a block
%              Toeplitz matrix T and a given right hand side matrix B.
%
%   task = 11: Compute the minimum norm solution of T'*Y = C for a given
%              right hand side matrix C.
%
%   task = 12: Compute X and Y.
%
%   task = 13: Compute the matrix-vector products X = T*B.
%
% Input parameters: 
%   task  - integer option to determine the computation to perform as
%           described above.
%   TC    - real M*K-by-L matrix, containing the first block column of
%           the block Toeplitz matrix T with K-by-L blocks.
%   TR    - real K-by-N*L matrix, containing the first block row of the
%           block Toeplitz matrix T with K-by-L blocks. If the first
%           block of TR differs from the first block of TC a warning
%           will be displayed that the first block of TC will be used
%           as block diagonal element of the matrix T.
%   B     - real right hand side M*K-by-NCB (task < 13) / N*L-by-NCB
%           (task = 13) matrix.
%   C     - real right hand side N*L-by-NCC matrix.
%   tol   - real tolerance used to estimate the numerical rank of T.
%           (See Section METHOD of the SLICOT Library routine MB02JX,
%           argument TOL1.)
%   m     - integer containing the number of block rows of T.
%   n     - integer containing the number of block columns of T.
%   TCB   - real (ML+1)*K-by-L matrix, containing the leading nonzero
%           blocks of the first block column of the banded block
%           Toeplitz matrix T with K-by-L blocks.
%   TRB   - real K-by-(NU+1)*L matrix, containing the leading nonzero
%           blocks of the first block row of the banded block Toeplitz
%           matrix T with K-by-L blocks. If the first block of TRB
%           differs from the first block of TCB a warning will be
%           displayed that the first block of TCB will be used as block
%           diagonal element of the matrix T.
%
% Output parameters:
%
%   Q     - the M*K-by-RNK orthogonal factor satisfying T = Q*R' /
%           T*E = Q*R', where RNK is the numerical rank of T. Note that
%           if task<>5 then RNK = MIN(M*K,N*L).
%   R     - the N*L-by-RNK Cholesky factor of T'*T / (T*E)'*(T*E)
%           (task < 6),
%           the MIN(ML+NU+1,N)*L-by-MIN(M*K,N*L) matrix containing the
%           Cholesky factor of T'*T in banded storage format
%           (6 <= task <= 9).
%   X     - the N*L-by-NCB least-squares solution of T*X = B (task < 13),
%           the M*K-by-NCB matrix containing the matrix-vector products
%           T*B (task = 13).
%   Y     - the M*K-by-NCC minimum norm solution of T'*Y = C.
%   E     - MIN(M*K,N*L) vector recording the column pivoting performed.
%           If E(j) = k, then the j-th column of T*P was the k-th column
%           of T, where P is the permutation matrix.
%
% Comments:
%   1.    - If only least-squares/minimum norm solutions are desired
%           then task = 6..9 should be used for efficiency.
%   2.    - Note that for the computation of the least-squares/minimum
%           norm solution, a semi-normal equation approach is used which
%           does not yield numerically backward stable solutions.
%           Iterative refinement should be used if improved accuracy
%           is required.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   D. Kressner, TU Berlin, Aug. 2002.
%
% Revisions:
%   V. Sima, Research Institute for Informatics, Bucharest, Aug. 2002.
%
