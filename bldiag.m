% BLDIAG.F   - MEX-function for reducing a matrix to a block-diagonal
%              form using well-conditioned non-orthogonal similarity
%              transformations, based on SLICOT routine MB03RD.
%
%   [Ao,blsize(,Wr,Wi,Xo)] = BLDIAG(A(,flaga,sorta,bound,jobx(,X)
%                                    (,tol)))
%
%   BLDIAG reduces a matrix A (possibly in a real Schur form) to a
%   block-diagonal form using well-conditioned non-orthogonal
%   similarity transformations. A general matrix A is first reduced to
%   a real Schur form, using orthogonal transformations. The condition
%   numbers of the transformations used for further reduction to a
%   block-diagonal form are roughly bounded by bound*bound, where
%   bound is a given value, or is taken as 100, by default. The
%   transformations are optionally postmultiplied in a given matrix X.
%   The real Schur form is optionally ordered, so that clustered
%   eigenvalues are grouped in the same block.
%
%   Description of input parameters:
%   A      - the n-by-n matrix A. The matrix is assumed to be in a real
%            Schur form if flaga = 1.
%   flaga  - (optional) scalar indicating whether or not the given
%            matrix A is in a real Schur form, as follows:
%            = 0 :  the matrix is general (default);
%            = 1 :  the matrix is in a real Schur form.
%   sorta  - (optional) scalar indicating whether or not the diagonal
%            blocks of the real Schur form are reordered, as follows:
%            = 0 :  the diagonal blocks are not reordered (default);
%            = 1 :  the diagonal blocks are reordered before each
%                   step of reduction, so that clustered eigenvalues
%                   appear in the same block;
%            = 2 :  the diagonal blocks are not reordered, but the
%                   "closest-neighbour" strategy is used instead of the
%                   standard "closest to the mean" strategy;
%            = 3 :  the diagonal blocks are reordered before each
%                   step of reduction, and the "closest-neighbour"
%                   strategy is used (see METHOD).
%   bound  - (optional) real scalar defining the upper bound for the
%            infinity norm of elementary submatrices of the individual
%            transformations used for reduction (see below). bound >= 1.
%            A large value allows more ill-conditioned transformations.
%            Default:  bound = 100.
%   jobx   - (optional) scalar indicating whether or not the
%            transformations performed are accumulated, as follows:
%            = 0 :  the transformations are not accumulated (default);
%            = 1 :  the transformations are accumulated in X (the
%                   given matrix X, if flaga = 1, is updated).
%   X      - (optional) if flaga = 1 and jobx = 1, an n-by-n matrix
%            containing the given initial transformation matrix. If
%            flaga = 0 or jobx = 0, X must not be specified on input.
%            If flaga = 0 and jobx = 1, X is taken as an identity
%            matrix.
%   tol    - (optional) real scalar indicating the tolerance to be used
%            in the ordering of the diagonal blocks of the real Schur
%            form matrix. If tol > 0, then the given value of tol is
%            used as an absolute tolerance: a block i and a temporarily
%            fixed block 1 (the first block of the current trailing
%            submatrix to be reduced) are considered to belong to the
%            same cluster if their eigenvalues satisfy
%
%               | lambda_1 - lambda_i | <= tol.
%
%            If tol < 0, then the given value of tol is used for a
%            relative tolerance: a block i and a temporarily fixed
%            block 1 are considered to belong to the same cluster
%            if their eigenvalues satisfy, for j = 1, ..., N,
%
%               | lambda_1 - lambda_i | <= | tol | * max | lambda_j |.
%
%            If tol = 0 or tol is missing, then an implicitly computed,
%            default tolerance, defined by tol = sqrt( sqrt( EPS ) ) is
%            used instead, as a relative tolerance, where EPS is the
%            machine precision (see LAPACK Library routine DLAMCH).
%            If sorta = 0 or sorta = 2, this parameter is not used.
%
%   Description of output parameters:
%   Ao     - the n-by-n computed block-diagonal matrix, in real Schur
%            canonical form. The non-diagonal blocks are set to zero.
%   blsize - integer vector of length at most n containing the orders
%            of the resulting diagonal blocks of the matrix Ao.
%   Wr,Wi  - (optional) n-vectors of real and imaginary parts,
%            respectively, of the eigenvalues of the matrix A.
%   Xo     - (optional) if flaga = 1 and jobx = 1, the computed n-by-n
%            matrix containing the product of the given matrix X and the
%            transformation matrix that reduced A to block-diagonal
%            form. The transformation matrix is itself a product of
%            non-orthogonal similarity transformations (having elements
%            with magnitude less than or equal to bound). If flaga = 0
%            and jobx = 1, Xo is the computed transformation matrix that
%            reduced A to block-diagonal form (including the initial
%            reduction to a real Schur form).
%
%   Method
%   Consider first that sorta = 0 and assume that A is in a real Schur
%   form,
%
%          ( A    A   )
%          (  11   12 )
%      A = (          ),
%          ( 0    A   )
%          (       22 )
%
%   where initially A   is the first diagonal block of dimension 1-by-1
%                    11
%   or 2-by-2. An attempt is made to compute a transformation matrix X
%   of the form
%
%          ( I   P )
%      X = (       )                                               (1)
%          ( 0   I )
%
%   (partitioned as A), so that
%
%               ( A     0  )
%       -1      (  11      )
%      X  A X = (          ),
%               ( 0    A   )
%               (       22 )
%
%   and the elements of P do not exceed the value bound in magnitude.
%   An adaptation of the standard Bartels-Stewart method for solving
%   Sylvester equations, which controls the magnitude of the individual
%   elements of the computed solution, is used to obtain matrix P.
%   When this attempt failed, an 1-by-1 (or 2-by-2) diagonal block of
%   A  , whose eigenvalue(s) is (are) the closest to the mean of those
%    22
%   of A   is selected, and moved by orthogonal similarity
%       11
%   transformations in the leading position of A  ; the moved diagonal
%                                               22
%   block is then added to the block A  , increasing its order by 1
%                                     11
%   (or 2). Another attempt is made to compute a suitable transformation
%   matrix X with the new definitions of the blocks A   and A  . After a
%                                                    11      22
%   successful transformation matrix X has been obtained, it
%   postmultiplies the current transformation matrix (if jobx = 1), and
%   the whole procedure is repeated for the matrix A  .
%                                                   22
%   When sorta = 1, the diagonal blocks of the real Schur form are
%   reordered before each step of the reduction, so that each cluster
%   of eigenvalues, defined as specified in the definition of tol,
%   appears in adjacent blocks. The blocks for each cluster are merged
%   together, and the procedure described above is applied to the
%   larger blocks. Using the option sorta = 1 will usually provide
%   better efficiency than the standard option (sorta = 0), because
%   there could be no or few unsuccessful attempts to compute individual
%   transformation matrices X of the form (1). However, the resulting
%   dimensions of the blocks are usually larger; this could make
%   subsequent calculations less efficient.
%
%   When sorta = 2 or 3, the procedure is similar to that for sorta = 0
%   or 1, respectively, but the block of A   whose eigenvalue(s) is
%                                         22
%   (are) the closest to those of A   (not to their mean) is selected
%                                  11
%   and moved to the leading position of A  . This is called the
%                                         22
%   "closest-neighbour" strategy.
%
%   Comments
%   The individual non-orthogonal transformation matrices used in the
%   reduction of A in a real Schur form to a block-diagonal form could
%   have condition numbers of the order bound*bound. This does not
%   guarantee that their product is well-conditioned enough.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003.
%
% Revisions:
%   -
%
