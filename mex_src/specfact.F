#include "fintrf.h"
C SPECFACT.F - Gateway function for computing the spectral factorization
C              of a real polynomial, arising from optimality problems,
C              using SLICOT routines SB08MD and SB08ND.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [E(,res,B)] = specfact(task,A(,form))
C
C Purpose:
C   1. To compute a real polynomial E(s) such that
C
C        (a)  E(-s) * E(s) = A(-s) * A(s) and
C        (b)  E(s) is stable - that is, all the zeros of E(s) have
C             non-positive real parts,
C
C   which corresponds to computing the spectral factorization of the
C   real polynomial A(s) arising from continuous optimality problems.
C
C   The input polynomial may be supplied either in the form
C
C   A(s) = a(0) + a(1) * s + ... + a(DA) * s**DA
C
C   or as
C
C   B(s) = A(-s) * A(s)
C        = b(0) + b(1) * s**2  + ... + b(DA) * s**(2*DA)             (1)
C
C   2. To compute a real polynomial E(z) such that
C
C        (a)  E(1/z) * E(z) = A(1/z) * A(z) and
C        (b)  E(z) is stable - that is, E(z) has no zeros with modulus
C             greater than 1,
C
C   which corresponds to computing the spectral factorization of the
C   real polynomial A(z) arising from discrete optimality problems.
C
C   The input polynomial may be supplied either in the form
C
C   A(z) = a(0) + a(1) * z + ... + a(DA) * z**DA
C
C   or as
C
C   B(z) = A(1/z) * A(z)
C        = b(0) + b(1) * (z + 1/z) + ... + b(DA) * (z**DA + 1/z**DA) (2)
C
C Input parameters:
C   task   - integer specifying the computations to be performed.
C            = 1 :  compute the spectral factorization for the
C                   continuous optimality problems;
C            = 2 :  compute the spectral factorization for the
C                   discrete optimality problems.
C   A      - a (DA+1)-vector containing either the coefficients of the
C            polynomial A in increasing powers of the indeterminate s
C            or z, if form = 1, or the coefficients b(0), ..., b(DA) of
C            the polynomial B in the formulas (1) or (2), if form = 2.
C   form   - (optional) indicates whether the coefficients of A or B are
C            supplied, as follows:
C            = 1 :  the coefficients of A are supplied;
C            = 2 :  the coefficients of B are supplied.
C            Default: form = 1.
C
C Output parameters:
C   E      - a (DA+1)-vector containing the coefficients of the spectral
C            factor E in increasing powers of s or z.
C   res    - (optional) an estimate of the accuracy with which the
C            coefficients of the polynomial E have been computed (see
C            also Method and Comments).
C   B      - (optional) the (DA+1)-vector containing the coefficients
C            b(0), ..., b(DA) of the polynomial B in the formulas (1)
C            or (2).
C
C Method:
C       _                                               _
C   Let A(s) be the conjugate polynomial of A(s), i.e., A(s) = A(-s), if
C                     _
C   task = 1, and let A(z) be the conjugate polynomial of A(z), i.e.,
C   _
C   A(z) = A(1/z), if task = 2. The method used is based on applying the
C   Newton-Raphson iteration to the function
C             _       _
C      F(e) = A * A - e * e,
C
C   which leads to the iteration formulae:
C
C      _(i)   (i)  _(i)   (i)     _      )
C      q   * x   + x   * q    = 2 A * A  )
C                                        )   for i = 0, 1, 2,...
C       (i+1)    (i)   (i)               )
C      q     = (q   + x   )/2            )
C
C                                (0)         DA
C   For task = 1, starting from q   = (1 + s)   (which has no zeros in
C                                                           (1)   (2)
C   the closed right half-plane), the sequence of iterates q   , q   ,
C   ..., converges to a solution of F(e) = 0 which has no zeros in the
C   open right half-plane. Similarly, for task = 2, the iteration starts
C   from
C
C         (0)                                        DA
C        q   (z) = (b(0) + b(1) * z + ... + b(DA) * z  ) / SQRT( b(0))
C
C   which is a Hurwitz polynomial that has no zeros in the closed unit
C                     (i)
C   circle. Then lim q   = e, the convergence is uniform and e is a
C   Hurwitz polynomial.
C
C   The iterates satisfy the following conditions, if task = 1:
C
C            (i)
C      (a)  q   is a stable polynomial (no zeros in the closed right
C           half-plane) and
C
C            (i)        (i-1)
C      (b)  q   (1) <= q     (1),
C
C   or, if task = 2:
C            (i)
C      (a)  q    has no zeros in the closed unit circle,
C            (i)        (i-1)
C      (b)  q   (0) <= q     (0) and
C
C            DA   (i) 2    DA     2
C      (c)  SUM (q   )  - SUM (A )  >= 0.
C           k=0   k       k=0   k
C                                   (i)
C   The iterative process stops if q    violates (a) or (b) (if
C   task = 1), or (a), (b) or (c) (if task = 2), or if the condition
C                     _(i) (i)  _
C      (d)  RES  = ||(q   q   - A A)|| < tol,
C
C   is satisfied, where || . || denotes the largest coefficient of
C                   _(i) (i)  _
C   the polynomial (q   q   - A A) and tol is an estimate of the
C                                                  _(i)  (i)
C   rounding error in the computed coefficients of q    q   .  If
C                                          (i-1)
C   condition (a) or (b) is violated then q      is taken; if
C                                   (i)
C   condition (c) is violated then q    is used. If there is no
C   convergence after 30 iterations, then an error indicator is set, and
C   the value of res may indicate whether or not the last computed
C   iterate is close to the solution.
C
C   If task = 1 and form = 2, then it is possible that the equation
C   e(-s) * e(s) = B(s) has no real solution, which will be the case
C   if A(1) < 0 or if ( -1)**DA * A(DA+1) < 0.
C                                                       (0)
C   If task = 2 and form = 2, then it is possible that q    is not a
C   Hurwitz polynomial, in which case the equation e(1/z) * e(z) = B(z)
C   has no real solution.
C
C Comments:
C   1. In order for the problem e(-s) * e(s) = B(s) (for task = 1) to
C   have a real solution e(s), it is necessary and sufficient that
C   B(j*omega) >= 0 for any purely imaginary argument j*omega.
C   2. The conditioning of the problem if task = 1 depends upon the
C   distance of the zeros of A(s) from the imaginary axis and on their
C   multiplicity. For a well-conditioned problem the accuracy of the
C   computed coefficients of E(s) is of the order of res. However, for
C   problems with zeros near the imaginary axis or with multiple zeros,
C   the value of res may be an overestimate of the true accuracy.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009,
C   Dec. 2010, July 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
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
      CHARACTER         ACONA
      INTEGER           DA, INFO, LDWORK
      DOUBLE PRECISION  RES
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      mwPointer         A,  DWORK, E
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           FORM, ITASK, LA
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          SB08MD, SB08ND
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SPECFACT requires at least 2 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SPECFACT requires at least 1 output argument' )
      END IF
C
C   task, A(DA), form.
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR.
     $     mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'TASK must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'TASK must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      ITASK = TEMP
      IF ( ITASK.LT.1 .OR. ITASK.GT.2 ) THEN
         CALL mexErrMsgTxt
     $    ( 'ITASK has 1 or 2 the only admissible values')
      END IF
C
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 .OR.
     $     MIN( mxGetM( PRHS(2) ), mxGetN( PRHS(2) ) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real polynomial' )
      END IF
      LA = MAX( mxGetM( PRHS(2) ), mxGetN( PRHS(2) ) )
      DA = LA - 1
      IF ( DA.LT.0 ) THEN
         CALL mexErrMsgTxt( 'DA must be at least 0' )
      END IF
C
      IF ( NRHS.GE.3 ) THEN
         IF ( mxGetM( PRHS(3) ).NE.1 .OR.
     $        mxGetN( PRHS(3) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'FORM must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(3) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'FORM must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TEMP, 1 )
         FORM = TEMP
         IF ( FORM.LT.1 .OR. FORM.GT.2 ) THEN
            CALL mexErrMsgTxt
     $       ( 'FORM has 1 or 2 the only admissible values')
         END IF
      ELSE
         FORM = 1
      END IF
C
      IF ( FORM.EQ.1 ) THEN
         ACONA = 'A'
      ELSE
         ACONA = 'B'
      END IF
C
C Determine the lenghts of working arrays.
C
      LDWORK = 5*DA + 5
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      A     = mxCalloc(     LA, 8 )
      DWORK = mxCalloc( LDWORK, 8 )
      E     = mxCalloc(     LA, 8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( A ), LA )
C
C Do the actual computations.
C
      IF ( ITASK.EQ.1 ) THEN
         CALL SB08MD( ACONA, DA, %VAL( A ), RES, %VAL( E ),
     $                %VAL( DWORK ), LDWORK, INFO )
      ELSE
         CALL SB08ND( ACONA, DA, %VAL( A ), RES, %VAL( E ),
     $                %VAL( DWORK ), LDWORK, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( LA, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( E ), mxGetPr( PLHS(1) ), LA )
         IF ( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( RES, mxGetPr( PLHS(2) ), 1 )
         END IF
         IF ( NLHS.GE.3 ) THEN
            PLHS(3) = mxCreateDoubleMatrix( LA, 1, 0 )
            CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS(3) ),
     $                             LA )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      CALL mxFree( A     )
      CALL mxFree( DWORK )
      CALL mxFree( E     )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( ITASK.EQ.1 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB08MD'')'
     $           ) INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM SB08ND'')'
     $           ) INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SPECFACT ***
      END
