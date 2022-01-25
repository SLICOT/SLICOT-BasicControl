#include "fintrf.h"
C BLDIAG.F   - Gateway function to reduce a matrix to a block-diagonal
C              form using well-conditioned non-orthogonal similarity
C              transformations, based on SLICOT routine MB03RD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [Ao,blsize(,Wr,Wi,Xo)] = BLDIAG(A(,flaga,sorta,bound,jobx(,X)
C                                    (,tol)))
C
C Purpose:
C   To reduce a matrix A (possibly in a real Schur form) to a
C   block-diagonal form using well-conditioned non-orthogonal
C   similarity transformations. A general matrix A is first reduced to
C   a real Schur form, using orthogonal transformations. The condition
C   numbers of the transformations used for further reduction to a
C   block-diagonal form are roughly bounded by bound*bound, where
C   bound is a given value, or is taken as 100, by default. The
C   transformations are optionally postmultiplied in a given matrix X.
C   The real Schur form is optionally ordered, so that clustered
C   eigenvalues are grouped in the same block.
C
C Input parameters:
C   A      - the n-by-n matrix A. The matrix is assumed to be in a real
C            Schur form if flaga = 1.
C   flaga  - (optional) scalar indicating whether or not the given
C            matrix A is in a real Schur form, as follows:
C            = 0 :  the matrix is general (default);
C            = 1 :  the matrix is in a real Schur form.
C   sorta  - (optional) scalar indicating whether or not the diagonal
C            blocks of the real Schur form are reordered, as follows:
C            = 0 :  the diagonal blocks are not reordered (default);
C            = 1 :  the diagonal blocks are reordered before each
C                   step of reduction, so that clustered eigenvalues
C                   appear in the same block;
C            = 2 :  the diagonal blocks are not reordered, but the
C                   "closest-neighbour" strategy is used instead of the
C                   standard "closest to the mean" strategy;
C            = 3 :  the diagonal blocks are reordered before each
C                   step of reduction, and the "closest-neighbour"
C                   strategy is used (see METHOD).
C   bound  - (optional) real scalar defining the upper bound for the
C            infinity norm of elementary submatrices of the individual
C            transformations used for reduction (see below). bound >= 1.
C            A large value allows more ill-conditioned transformations.
C            Default:  bound = 100.
C   jobx   - (optional) scalar indicating whether or not the
C            transformations performed are accumulated, as follows:
C            = 0 :  the transformations are not accumulated (default);
C            = 1 :  the transformations are accumulated in X (the
C                   given matrix X, if flaga = 1, is updated).
C   X      - (optional) if flaga = 1 and jobx = 1, an n-by-n matrix
C            containing the given initial transformation matrix. If
C            flaga = 0 or jobx = 0, X must not be specified on input.
C            If flaga = 0 and jobx = 1, X is taken as an identity
C            matrix.
C   tol    - (optional) real scalar indicating the tolerance to be used
C            in the ordering of the diagonal blocks of the real Schur
C            form matrix. If tol > 0, then the given value of tol is
C            used as an absolute tolerance: a block i and a temporarily
C            fixed block 1 (the first block of the current trailing
C            submatrix to be reduced) are considered to belong to the
C            same cluster if their eigenvalues satisfy
C
C               | lambda_1 - lambda_i | <= tol.
C
C            If tol < 0, then the given value of tol is used for a
C            relative tolerance: a block i and a temporarily fixed
C            block 1 are considered to belong to the same cluster
C            if their eigenvalues satisfy, for j = 1, ..., N,
C
C               | lambda_1 - lambda_i | <= | tol | * max | lambda_j |.
C
C            If tol = 0 or tol is missing, then an implicitly computed,
C            default tolerance, defined by tol = sqrt( sqrt( EPS ) ) is
C            used instead, as a relative tolerance, where EPS is the
C            machine precision (see LAPACK Library routine DLAMCH).
C            If sorta = 0 or sorta = 2, this parameter is not used.
C
C Output parameters:
C   Ao     - the n-by-n computed block-diagonal matrix, in real Schur
C            canonical form. The non-diagonal blocks are set to zero.
C   blsize - integer vector of length at most n containing the orders
C            of the resulting diagonal blocks of the matrix Ao.
C   Wr,Wi  - (optional) n-vectors of real and imaginary parts,
C            respectively, of the eigenvalues of the matrix A.
C   Xo     - (optional) if flaga = 1 and jobx = 1, the computed n-by-n
C            matrix containing the product of the given matrix X and the
C            transformation matrix that reduced A to block-diagonal
C            form. The transformation matrix is itself a product of
C            non-orthogonal similarity transformations (having elements
C            with magnitude less than or equal to bound). If flaga = 0
C            and jobx = 1, Xo is the computed transformation matrix that
C            reduced A to block-diagonal form (including the initial
C            reduction to a real Schur form).
C
C Method
C   Consider first that sorta = 0 and assume that A is in a real Schur
C   form,
C
C          ( A    A   )
C          (  11   12 )
C      A = (          ),
C          ( 0    A   )
C          (       22 )
C
C   where initially A   is the first diagonal block of dimension 1-by-1
C                    11
C   or 2-by-2. An attempt is made to compute a transformation matrix X
C   of the form
C
C          ( I   P )
C      X = (       )                                               (1)
C          ( 0   I )
C
C   (partitioned as A), so that
C
C               ( A     0  )
C       -1      (  11      )
C      X  A X = (          ),
C               ( 0    A   )
C               (       22 )
C
C   and the elements of P do not exceed the value bound in magnitude.
C   An adaptation of the standard Bartels-Stewart method for solving
C   Sylvester equations, which controls the magnitude of the individual
C   elements of the computed solution, is used to obtain matrix P.
C   When this attempt failed, an 1-by-1 (or 2-by-2) diagonal block of
C   A  , whose eigenvalue(s) is (are) the closest to the mean of those
C    22
C   of A   is selected, and moved by orthogonal similarity
C       11
C   transformations in the leading position of A  ; the moved diagonal
C                                               22
C   block is then added to the block A  , increasing its order by 1
C                                     11
C   (or 2). Another attempt is made to compute a suitable transformation
C   matrix X with the new definitions of the blocks A   and A  . After a
C                                                    11      22
C   successful transformation matrix X has been obtained, it
C   postmultiplies the current transformation matrix (if jobx = 1), and
C   the whole procedure is repeated for the matrix A  .
C                                                   22
C   When sorta = 1, the diagonal blocks of the real Schur form are
C   reordered before each step of the reduction, so that each cluster
C   of eigenvalues, defined as specified in the definition of tol,
C   appears in adjacent blocks. The blocks for each cluster are merged
C   together, and the procedure described above is applied to the
C   larger blocks. Using the option sorta = 1 will usually provide
C   better efficiency than the standard option (sorta = 0), because
C   there could be no or few unsuccessful attempts to compute individual
C   transformation matrices X of the form (1). However, the resulting
C   dimensions of the blocks are usually larger; this could make
C   subsequent calculations less efficient.
C
C   When sorta = 2 or 3, the procedure is similar to that for sorta = 0
C   or 1, respectively, but the block of A   whose eigenvalue(s) is
C                                         22
C   (are) the closest to those of A   (not to their mean) is selected
C                                  11
C   and moved to the leading position of A  . This is called the
C                                         22
C   "closest-neighbour" strategy.
C
C Comments
C   The individual non-orthogonal transformation matrices used in the
C   reduction of A in a real Schur form to a block-diagonal form could
C   have condition numbers of the order bound*bound. This does not
C   guarantee that their product is well-conditioned enough.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009,
C   Dec. 2010.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HUNDRD
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, HUNDRD = 1.0D2 )
C
C .. Mex-file interface parameters ..
      mwPointer         PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      mwPointer         mxCalloc, mxCreateDoubleMatrix,
     $                  mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         JOBX, SORT
      INTEGER           INFO, LDA, LDX, N, NBLCKS
      DOUBLE PRECISION  PMAX, TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: BLSIZE(:)
      DOUBLE PRECISION, ALLOCATABLE :: DW( : )
      mwPointer         A, DWORK, WI, WR, X
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           BWORK(1)
      CHARACTER         JOBVS
      CHARACTER*120     TEXT
      INTEGER           I, IFLAG, IP, ISORT, JOB, LDWORK, SDIM
      DOUBLE PRECISION  TEMP
C
C .. External Functions ..
      LOGICAL           SELECT
      EXTERNAL          SELECT
C
C .. External subroutines ..
      EXTERNAL          DGEES, MB03RD
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN, MOD
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'BLDIAG requires at least 1 input argument' )
      ELSE IF ( NLHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'BLDIAG requires at least 2 output arguments' )
      END IF
C
C   A(nxn)(, flaga, sorta, bound, jobx(, X(nxn)) (,tol))
C
      N =  mxGetM( PRHS(1) )
      IF ( mxGetN( PRHS(1) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
C   flaga
C
      IP = 2
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'FLAGA must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'FLAGA must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IFLAG = TEMP
         IF ( IFLAG.LT.0 .OR. IFLAG.GT.1 ) THEN
            CALL mexErrMsgTxt
     $        ( 'FLAGA has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
      ELSE
         IFLAG = 0
      END IF
C
C   sorta
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'SORTA must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'SORTA must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ISORT = TEMP
         IF ( ISORT.LT.0 .OR. ISORT.GT.3 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SORTA has 0, 1, 2, or 3 the only admissible values')
         END IF
         IP = IP + 1
      ELSE
         ISORT = 0
      END IF
C
      IF ( ISORT.EQ.0 ) THEN
         SORT = 'N'
      ELSEIF ( ISORT.EQ.1 ) THEN
         SORT = 'S'
      ELSEIF ( ISORT.EQ.2 ) THEN
         SORT = 'C'
      ELSE
         SORT = 'B'
      END IF
C
C   bound
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'BOUND must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'BOUND must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), PMAX, 1 )
         IF ( PMAX.LT.ONE ) THEN
            CALL mexErrMsgTxt
     $           ( 'BOUND must be larger than or equal to 1' )
         END IF
         IP = IP + 1
      ELSE
         PMAX = HUNDRD
      END IF
C
C   jobx
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'JOBX must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'JOBX must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         JOB = TEMP
         IF ( JOB.LT.0 .OR. JOB.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'JOBX has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
      ELSE
         JOB = 0
      END IF
C
      IF ( JOB.EQ.0 ) THEN
         JOBX  = 'N'
         JOBVS = 'N'
      ELSE
         JOBX  = 'U'
         JOBVS = 'V'
      END IF
C
C   X
C
      IF ( NRHS.GE.IP .AND. IFLAG*JOB.GT.0 ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.N ) THEN
            WRITE( TEXT, '(''X must have '',I5,'' rows'')' ) N
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxGetN( PRHS(IP) ).NE.N ) THEN
            WRITE( TEXT, '(''X must have '',I5,'' columns'')' ) N
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'X must be a real matrix' )
         END IF
         IP = IP + 1
      END IF
C
C   tol
C
      IF ( MOD( ISORT, 2 ).EQ.0 .AND. NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TOL, 1 )
         IP = IP + 1
      ELSE
         TOL = ZERO
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA = MAX( 1, N )
      IF ( JOB.GT.0 ) THEN
         LDX = LDA
      ELSE
         LDX = 1
      END IF
C
C   ldwork
C
      IF ( IFLAG.EQ.0 ) THEN
         LDWORK = MAX( 1, 3*N )
      ELSE
         LDWORK = N
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( BLSIZE( N ), DW( N ) )
      A     = mxCalloc(  LDA*N,   8 )
      DWORK = mxCalloc( LDWORK,   8 )
      WI    = mxCalloc(      N,   8 )
      WR    = mxCalloc(      N,   8 )
      X     = mxCalloc(  LDX*LDX, 8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), %VAL( A ), N*N )
      IF ( IFLAG*JOB.GT.0 .AND. IP.GE.6 )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), %VAL( X ), N*N )
C
C Do the actual computations.
C
      INFO = 0
      IF ( IFLAG.EQ.0 )
     $   CALL DGEES( JOBVS, 'No order', SELECT, N, %VAL( A ), LDA, SDIM,
     $               %VAL( WR ), %VAL( WI ), %VAL( X ), LDX,
     $               %VAL( DWORK ), LDWORK, BWORK, INFO )
      IF ( INFO.EQ.0 )
     $   CALL MB03RD( JOBX, SORT, N, PMAX, %VAL( A ), LDA, %VAL( X ),
     $                LDX, NBLCKS, BLSIZE, %VAL( WR ), %VAL( WI ), TOL,
     $                %VAL( DWORK ), INFO )
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         IP = 2
         PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
         CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS(1) ), N*N )
         IF ( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( NBLCKS, MIN( N, 1 ), 0 )
            DO 10 I = 1, NBLCKS
               DW(I) = BLSIZE(I)
   10       CONTINUE
            CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(2) ),
     $                             NBLCKS*MIN( N, 1 ) )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, MIN( N, 1 ), 0 )
            CALL mxCopyReal8ToPtr( %VAL( WR ), mxGetPr( PLHS(IP) ),
     $                             N*MIN( N, 1 ) )
            IP = IP + 1
            PLHS(IP) = mxCreateDoubleMatrix( N, MIN( N, 1 ), 0 )
            CALL mxCopyReal8ToPtr( %VAL( WI ), mxGetPr( PLHS(IP) ),
     $                             N*MIN( N, 1 ) )
            IP = IP + 1
         END IF
C
         IF ( JOB.GT.0 .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( X ), mxGetPr( PLHS(IP) ),
     $                             N*N )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( BLSIZE, DW )
      CALL mxFree( A     )
      CALL mxFree( DWORK )
      CALL mxFree( WI    )
      CALL mxFree( WR    )
      CALL mxFree( X     )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DGEES'')' ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of BLDIAG ***
      END
