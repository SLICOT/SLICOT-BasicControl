%FSTOEP  MEX-function for factoring symmetric positive definite
%        (block) Toeplitz matrices and/or solving associated linear
%        systems using SLICOT routines MB02CD, MB02DD, MB02ED.
%
%                   [R(,X)] = fstoep(task,T(,B))
%              [G(,Li,R,X)] = fstoep(task,T(,B))
%           [G(,Li),R,H,CS] = fstoep(task,T)
%                 [Ru(,Xa)] = fstoep(task,Ta,H,CS,G,R(,Ba))
%   [Gu(,Liu,Ru,Hu,CSu,Xa)] = fstoep(task,Ta,H,CS,G,R(,Li,Ba))
%                       [X] = fstoep(task,T,B)
%
%   task =  1 :                 [R] = fstoep( 1,T)
%                             [R,X] = fstoep( 1,T,B)
%   task =  2 :                 [G] = fstoep( 2,T)
%   task =  3 :               [G,R] = fstoep( 3,T)
%                           [G,R,X] = fstoep( 3,T,B)
%   task =  4 :              [G,Li] = fstoep( 4,T)
%                          [G,Li,X] = fstoep( 4,T,B)
%   task =  5 :            [G,Li,R] = fstoep( 5,T)
%                        [G,Li,R,X] = fstoep( 5,T,B)
%   task =  6 :          [G,R,H,CS] = fstoep( 6,T)
%   task =  7 :       [G,Li,R,H,CS] = fstoep( 7,T)
%   task =  8 :                [Ru] = fstoep( 8,Ta,H,CS,G,R)
%                           [Ru,Xa] = fstoep( 8,Ta,H,CS,G,R,Ba)
%   task =  9 :             [Gu,Ru] = fstoep( 9,Ta,H,CS,G,R)
%                        [Gu,Ru,Xa] = fstoep( 9,Ta,H,CS,G,R,Ba)
%   task = 10 :      [Gu,Ru,Hu,CSu] = fstoep(10,Ta,H,CS,G,R)
%                [Gu,Liu,Ru,Hu,CSu] = fstoep(10,Ta,H,CS,G,R,Li)
%             [Gu,Liu,Ru,Hu,CSu,Xa] = fstoep(10,Ta,H,CS,G,R,Li,Ba)
%   task = 11 :                 [X] = fstoep(11,T,B)
%
%   FSTOEP factors a symmetric positive definite block Toeplitz matrix BT
%   and/or its inverse, and/or computes/updates the generator of its
%   inverse, given the first block row / column T of BT, and/or solves
%   the associated linear systems X*BT = B / BT*X = B.
%   For task <= 5 or task = 11, the first block row / column of BT
%   contains T and has the dimension K-by-N*K / N*K-by-K, while for
%   task >= 6 and task <= 10, the first block row / column of the block
%   Toeplitz matrix, denoted BTA, contains [ T Ta ] / [ T  ] and has the
%   dimension K-by-(N+M)*K / (N+M)*K-by-K.            [ Ta ]
%
%   task =  1: Compute R, the Cholesky factor of BT, i.e., R is
%              upper / lower triangular, so that R'*R = BT / R*R' = BT;
%
%   task =  2: Compute the generator G of inv(BT), of dimension
%              2*K-by-N*K / N*K-by-2*K;
%
%   task =  3: Compute both G and R;
%
%   task =  4: Compute both G and Li, where Li is the Cholesky factor
%              of inv(BT), i.e., Li is lower / upper triangular, so that
%              Li'*Li = inv(BT) / Li*Li' = inv(BT);
%
%   task =  5: Compute both G, Li, and R;
%
%   task =  6: Compute G and R, and deliver the information needed to
%              update a Cholesky factorization.
%
%   task =  7: Compute G, Li, and R, and deliver the information needed
%              to update a Cholesky factorization.
%
%   task =  8: Compute updated R, Ru, given additional blocks of
%              data, Ta, and previous factorization results.
%
%   task =  9: Compute updated G and R, Gu and Ru, given additional
%              blocks of data, Ta, and previous factorization results.
%
%   task = 10: Compute updated G, (Li,) and R, as well as details of the
%              transformations used (needed for subsequent updating),
%              given additional blocks of data, Ta, and previous 
%              factorization results.
%
%   task = 11: Solve BT*X = B or X*BT = B.
%
%   Note:      The linear systems BT*X = B or X*BT = B can also be
%              solved when task <= 5, task <> 2, by specifying an
%              additional input and output parameter.
%              Moreover, when task = 8, 9, or 10, one can solve the 
%              linear systems BTA*Xa = Ba or Xa*BTA = Ba, by specifying
%              an additional input and output parameter.
%
%   Description of input parameters:
%   task  - integer option to determine the computation to perform:
%           =  1 : compute the Cholesky factor R of the matrix BT;
%           =  2 : compute the generator G of inv(BT);
%           =  3 : compute both G and R;
%           =  4 : compute both G and Li, where Li is the Cholesky
%                  factor of inv(BT);
%           =  5 : compute both G, Li and R;
%           =  6 : compute G and R, and deliver the information needed
%                  to update the factorization;
%           =  7 : compute G, Li, and R, and deliver the information
%                  needed to update the factorization.
%           =  8 : compute updated R, given additional blocks of data;
%           =  9 : compute updated G and R, given additional data;
%           = 10 : compute updated G, (Li,) and R, as well as details of
%                  the transformations used, given additional data;
%           = 11 : solve X*BT = B or BT*X = B.
%   T     - real K-by-N*K / N*K-by-K matrix, containing (part) of the
%           first block row / column of the symmetric positive definite
%           block Toeplitz matrix BT (BTA), with K-by-K blocks.
%   B     - real right hand side NRHS-by-N*K / N*K-by-NRHS matrix.
%   Ta    - real K-by-M*K / M*K-by-K matrix containing additional data
%           of the first block row / column of the symmetric positive
%           definite block Toeplitz matrix BTA.
%   H     - real K-by-N*K / N*K-by-K matrix, containing part of the
%           information about the Householder transformations used.
%   CS    - real 3*(N-1)*K vector containing further details on the
%           Householder transformations and hyperbolic rotations used.
%   G     - real 2*K-by-N*K / N*K-by-2*K generator G of inv(BT).
%   R     - real N*K-by-N*K Cholesky factor of BT.
%   Li    - real N*K-by-N*K Cholesky factor of inv(BT).
%   Ba    - real right hand side NRHS-by-(N+M)*K / (N+M)*K-by-NRHS 
%           matrix.
%
%   Description of output parameters:
%   R     - the N*K-by-N*K Cholesky factor of BT.
%   G     - the 2*K-by-N*K / N*K-by-2*K generator G of inv(BT).
%   Li    - the N*K-by-N*K Cholesky factor of inv(BT).
%   X     - the NRHS-by-N*K / N*K-by-NRHS solution matrix of the system
%           X*BT = B / BT*X = B.
%   H     - the K-by-N*K / N*K-by-K matrix, containing part of the
%           information about the Householder transformations used.
%   CS    - the 3*(N-1)*K vector containing further details on the
%           Householder transformations and hyperbolic rotations used.
%   Ru    - the (N+M)*K-by-(N+M)*K Cholesky factor of BTA.
%   Gu    - the 2*K-by-(N+M)*K / (N+M)*K-by-2*K generator of inv(BTA).
%   Liu   - the (N+M)*K-by-(N+M)*K Cholesky factor of inv(BTA).
%   Hu    - the K-by-(N+M)*K / (N+M)*K-by-K matrix containing part of
%           the information about the Householder transformations used.
%   CSu   - the 3*(N+M-1)*K vector containing further details on the
%           Householder transformations and hyperbolic rotations used.
%   Xa    - the NRHS-by-(N+M)*K / (N+M)*K-by-NRHS solution matrix of the
%           system Xa*BTA = Ba / BTA*Xa = Ba.
%
% Comments:
%   1. If the second dimension of T (Ta) is larger than the first one,
%      but it is not a multiple of its first dimension, then the number
%      of columns considered is the largest possible multiple of K.
%      If the first dimension of T (Ta) is larger than the second one,
%      but it is not a multiple of its second dimension, then the number
%      of rows considered is the largest possible multiple of K.
%   2. If one dimension of T is zero, then both K and N are set to 0.
%   3. If task >= 6 and task <= 10, T and Ta must have the same shape.
%   4. If task = 6, 7, or 10, on output, G or Gu contain information 
%      needed to perform (further) updates. To get the true generator, 
%      the submatrix (K+1:2*K,1:K) / (1:K,K+1:2*K) of G or Gu must be
%      set to zeros(K).
%   5. If N = 1, then the returned results satisfy R'*R = BT
%      (or Ru'*Ru = BTA), Li'*Li = inv(BT) (or Liu'*Liu = inv(BTA)),
%      and G has the size 2*K-by-K (or Gu has the size 2*K-by-(1+M)*K).
%   6. If only the solution of the system X*BT = B / BT*X = B is
%      desired, then task = 11 should be used, for efficiency.
%

%   RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2000.
%
%   Revisions:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000.
%
