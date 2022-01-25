function [Ao,blsize,W,Xo] = bdiag(A,sorta,bound,flaga,jobx,X,tol)
% BDIAG  reduces a matrix A to a block-diagonal form using
%        well-conditioned non-orthogonal similarity transformations.
%        A general matrix A is first reduced to a real Schur form,
%        using orthogonal transformations. The condition numbers of
%        the transformations used for further reduction are bounded.
%
%        [Ao,BLSIZE] = BDIAG(A)  computes the block-diagonal matrix Ao,
%        in real Schur canonical form; the integer vector BLSIZE
%        contains the orders of the resulting diagonal blocks.
%
%        [Ao,BLSIZE,W] = BDIAG(A,SORTA,BOUND,FLAGA)  returns in W the
%        eigenvalues of the matrix A. The additional input arguments
%        specify various options.
%
%        SORTA indicates whether or not the diagonal blocks of the
%        real Schur form are reordered, as follows:
%        SORTA = 0 :  the diagonal blocks are not reordered;
%        SORTA = 1 :  the diagonal blocks are reordered before each
%                     step of reduction, so that clustered eigenvalues
%                     appear in the same block;
%        SORTA = 2 :  the diagonal blocks are not reordered, but the
%                     "closest-neighbour" strategy is used instead of
%                     the standard "closest to the mean" strategy;
%        SORTA = 3 :  the diagonal blocks are reordered before each
%                     step of reduction, and the "closest-neighbour"
%                     strategy is used.
%        Default: SORTA = 0.
%
%        BOUND defines the upper bound for the infinity norm of
%        elementary submatrices of the individual transformations used
%        for reduction. BOUND >= 1.
%        A large value allows more ill-conditioned transformations.
%        Default:  BOUND = 100.
%
%        FLAGA indicates whether or not the given matrix A is in a
%        real Schur form, as follows:
%        FLAGA = 0 :  the matrix is general;
%        FLAGA = 1 :  the matrix is in a real Schur form.
%        Default:  FLAGA = 0.
%
%        [Ao,BLSIZE,W,Xo] = BDIAG(A,SORTA,BOUND,FLAGA,JOBX,X,TOL)  
%        could return the transformation matrix Xo. Specifically:
%
%        JOBX indicates whether or not the transformations performed
%        are accumulated, as follows:
%        JOBX = 0 :  the transformations are not accumulated;
%        JOBX = 1 :  the transformations are accumulated in X (the
%                    given matrix X, if FLAGA = 1, is updated).
%        Default:  JOBX = 0.
%
%        If FLAGA = 1 and JOBX = 1, X contains the given initial
%        transformation matrix. If FLAGA = 0 or JOBX = 0, X must
%        not be specified on input. If FLAGA = 0 and JOBX = 1,
%        X is taken as an identity matrix.
%
%        TOL is the tolerance to be used in the ordering of the
%        diagonal blocks of the real Schur form matrix. If TOL > 0,
%        then the given value of TOL is used as an absolute tolerance:
%        a block i and a temporarily fixed block 1 (the first block of
%        the current trailing submatrix to be reduced) are considered
%        to belong to the same cluster if their eigenvalues satisfy
%
%               | lambda_1 - lambda_i | <= TOL.
%
%        If TOL < 0, then the given value of TOL is used for a
%        relative tolerance: a block i and a temporarily fixed
%        block 1 are considered to belong to the same cluster
%        if their eigenvalues satisfy, for all j,
%
%               | lambda_1 - lambda_i | <= | TOL | * max | lambda_j |.
%
%        If SORTA = 0 or SORTA = 2, this parameter is not used.
%        Default (also used if TOL = 0):  TOL = sqrt( sqrt( eps ) ).
%
%        If FLAGA = 1 and JOBX = 1, Xo contains the product of the
%        given matrix X and the transformation matrix that reduced A
%        to block-diagonal form. The transformation matrix is itself
%        a product of non-orthogonal similarity transformations
%        (having elements with magnitude less than or equal to BOUND).
%        If FLAGA = 0 and JOBX = 1, Xo is the computed transformation
%        matrix that reduced A to block-diagonal form (including the
%        initial reduction to a real Schur form).
%
%        See also BLDIAG
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, April 2003.
%
%        Revisions:
%        V. Sima, March 2009.
%

ni   = nargin;
nout = nargout;
if ni < 1 || nout < 2,
    error(['Usage: [Ao,BLSIZE]      = BDIAG(A)',  sprintf('\n'),...
           '       [Ao,BLSIZE,W]    = BDIAG(A,SORTA,BOUND,FLAGA)',  sprintf('\n'),...
           '       [Ao,BLSIZE,W,Xo] = BDIAG(A,SORTA,BOUND,FLAGA,JOBX,X,TOL)'])
end
% 
if ni == 1,
    sorta = 0;
    bound = 100;
    flaga = 0;
    jobx  = 0;
    tol   = 0;
elseif ni == 2,
    bound = 100;
    flaga = 0;
    jobx  = 0;
    tol   = 0;
elseif ni == 3,
    flaga = 0;
    jobx  = 0;
    tol   = 0;
elseif ni == 4,
    jobx  = 0;
    tol   = 0;
elseif flaga*jobx == 1 && ni == 6,
    tol   = 0;
elseif ni == 5,
    tol   = 0;
end
%
if nout == 2,
    [ Ao, blsize ] = bldiag( A, flaga, sorta, bound );
elseif nout == 3,
    [ Ao, blsize, Wr, Wi ] = bldiag( A, flaga, sorta, bound );
    W = Wr + 1i*Wi;
elseif nout == 4,
    if flaga*jobx == 1 && ni >= 6,
        [ Ao, blsize, Wr, Wi, Xo ] = bldiag( A, flaga, sorta, bound, jobx, X, tol );
    else
        [ Ao, blsize, Wr, Wi, Xo ] = bldiag( A, flaga, sorta, bound, jobx );
    end
    W = Wr + 1i*Wi;
end
%
% end bdiag
