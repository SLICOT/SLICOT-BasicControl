#include "fintrf.h"
C PERSCHUR.F - Gateway function to compute the periodic Schur
C              decomposition using SLICOT routines MB03VD, MB03VY,
C              and MB03WD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [Ao(,WR,WI,Z)] = perschur(A(,job,compz,index))
C
C   [Ao(,Z)]       = perschur(A, 0 (,compz,index))
C   [WR,WI(,Z)]    = perschur(A(,1  ,compz,index))
C   [Ao(,WR,WI,Z)] = perschur(A, 2 (,compz,index))
C
C Purpose:
C   To compute either the upper Hessenberg form or the Schur
C   decomposition and the eigenvalues of a product of n-by-n matrices,
C   A = A_1*A_2*...*A_p, p >= 1, using orthogonal transformations, and
C   without evaluating the product.
C
C Input parameters:
C   A     - the n-by-n-by-p matrices A_1, ..., A_p.
C   job   - (optional) scalar indicating how A should be transformed,
C           as follows:
C           = 0 :  compute the factors H_1, ..., H_p of the upper
C                  Hessenberg form H of A, i.e., H = H_1*H_2*...*H_p,
C                  where H_1 is upper Hessenberg and H_i are upper
C                  triangular, i = 2, ..., p; specifically, the
C                  transformations are
C                     Z_i'*A_i*Z_(i+1) = H_i, i = 1, ..., p-1, and
C                     Z_p'*A_p*Z_1     = H_p;
C           = 1 :  compute the eigenvalues only (default);
C           = 2 :  compute the factors T_1, ..., T_p of the full Schur
C                  form, T = T_1*T_2*...*T_p, where T_1 is in real Schur
C                  form and T_i are upper triangular, i = 2, ..., p;
C                  specifically, the transformations are
C                     Z_i'*A_i*Z_(i+1) = T_i, i = 1, ..., p-1, and
C                     Z_p'*A_p*Z_1     = T_p.
C   compz - (optional) scalar indicating whether or not to compute
C           the orthogonal matrices Z_1, ..., Z_p, as follows:
C           = 0 :  the matrices Z_1,..., Z_p are not computed (default);
C           = 1 :  the matrices Z_1,..., Z_p are computed, i = 1,..., p.
C   index - (optional) vector of length at most 4 containing
C           [ILO; IHI; ILOZ; IHIZ]. If index is specified, it is assumed
C           that all matrices A_j, j = 2, ..., p, are already upper
C           triangular in rows and columns 1:ILO-1 and IHI+1:n, and A_1
C           is upper Hessenberg in rows and columns 1:ILO-1 and IHI+1:n,
C           with A_1(ILO,ILO-1) = 0 (unless ILO = 1), and
C           A_1(IHI+1,IHI) = 0 (unless IHI = n).
C           Otherwise, ILO = 1, IHI = n.
C           1 <= ILO <= max(1,n); min(ILO,n) <= IHI <= n.
C           If compz = 1, ILOZ and IHIZ specify the rows of Z to which
C           the transformations must be applied.
C           1 <= ILOZ <= ILO; IHI <= IHIZ <= n.
C           Default: index = [1; n; 1; n].
C
C Output parameters:
C   Ao    - if job = 0, the computed n-by-n-by-p matrix H, containing
C           the factors H_1, ..., H_p of the periodic Hessenberg form;
C           if job = 2, the computed n-by-n-by-p matrix T, containing
C           the factors T_1, ..., T_p of the periodic Schur form;
C           if job = 1, Ao is not returned, since no useful information
C           is available.
C   WR    - if job >= 1, the n-vector of real parts of the computed
C           eigenvalues.
C   WI    - if job >= 1, the n-vector of imaginary parts of the computed
C           eigenvalues. Complex conjugate pairs appear consecutively.
C           If index is specified, and ILO > 1 and/or IHI < n, then only
C           the locations ILO to IHI are set in WR and WI.
C   Z     - (optional) if compz = 1, the n-by-n-by-p matrix Z,
C           containing the factors Z_i, i = 1 : p, which reduce the
C           matrix product to the upper Hessenberg form, or periodic
C           Schur form, according to job.
C           The transformations are applied only to the submatrices
C           Z_i(ILOZ:IHIZ,ILO:IHI), i = 1, ..., p.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, June 2002.
C
C Revisions:
C   V. Sima, July 2002, Apr. 2009, July 2012, Aug. 2012, Dec. 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
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
      CHARACTER         COMPZ, JOB
      INTEGER           IHI, IHIZ, ILO, ILOZ, INFO, LDA1, LDA2, LDTAU,
     $                  LDWORK, LDZ1, LDZ2, N, P
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      mwPointer         A, DWORK, TAU, WI, WR, Z
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           I, ICOMP, IJOB, IP, ITMP, LI, NCI, NRI
      DOUBLE PRECISION  TEMP, TMP(4)
C
C .. External subroutines ..
      EXTERNAL          DLACPY, MA02LD, MB03VD, MB03VY, MB03WD
C
C .. Intrinsic functions ..
      INTRINSIC         INT, MAX, MIN, MOD
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'PERSCHUR requires at least 1 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'PERSCHUR requires at least 1 output argument' )
      END IF
C
C   A(nxnxp) (, job, compz, index).
C
      N = mxGetM( PRHS(1) )
      IF ( N.EQ.0 ) THEN
         P = 1
      ELSE
         P = mxGetN( PRHS(1) )
         IF ( MOD( P, N ).NE.0 ) THEN
            WRITE( TEXT, '(''The length of A must be a multiple of '',
     $                     I7)' ) N*N
            CALL mexErrMsgTxt( TEXT )
         ELSE
            P = P / N
         END IF
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
C   job
C
      IP   = 2
      IJOB = 1
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'JOB must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ITMP = TEMP
         IF ( ITMP.LT.0 .OR. ITMP.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'JOB has 0, 1, or 2 the only admissible values' )
         END IF
         IJOB = ITMP
         IP   = IP + 1
      END IF
C
      IF ( IJOB.EQ.1 ) THEN
         JOB = 'E'
         IF ( NLHS.LT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'PERSCHUR requires at least 2 output arguments' )
         END IF
      ELSE IF ( IJOB.EQ.2 ) THEN
         JOB = 'S'
         IF ( NLHS.EQ.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'PERSCHUR requires at least 3 output arguments' )
         END IF
      END IF
C
C   compz
C
      ICOMP = 0
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'COMPZ must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'COMPZ must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ITMP = TEMP
         IF ( ITMP.LT.0 .OR. ITMP.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPZ has 0 or 1 the only admissible values' )
         END IF
         ICOMP = ITMP
         IP    = IP + 1
      END IF
C
      IF ( ICOMP.EQ.0 ) THEN
         COMPZ = 'N'
      ELSE
         COMPZ = 'V'
      END IF
C
C   index
C
      ILO  = 1
      IHI  = N
      ILOZ = 1
      IHIZ = N
      IF ( NRHS.GE.IP ) THEN
         NRI = mxGetM( PRHS(IP) )
         NCI = mxGetN( PRHS(IP) )
         LI  = NRI*NCI
         IF ( MIN( NRI, NCI ).NE.1 .AND. LI.GT.4 ) THEN
            WRITE( TEXT, '(''INDEX must be a vector with at most 4 '',
     $                     ''integer elements'')' )
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'INDEX must be an integer vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TMP, LI )
         IF ( LI.GE.1 ) ILO  = MAX( 1, MIN( N,   INT( TMP(1) ) ) )
         IF ( LI.GE.2 ) IHI  = MIN( N, MAX( ILO, INT( TMP(2) ) ) )
         IF ( LI.GE.3 ) ILOZ = MAX( 1, MIN( ILO, INT( TMP(3) ) ) )
         IF ( LI.GE.4 ) IHIZ = MIN( N, MAX( IHI, INT( TMP(4) ) ) )
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA1 = MAX( 1, N )
      LDA2 = LDA1
      IF ( ICOMP.EQ.0 ) THEN
         LDZ1 = 1
         LDZ2 = 1
      ELSE
         LDZ1 = LDA1
         LDZ2 = LDA1
      END IF
      LDTAU = MAX( 1, N - 1 )
C
C   ldwork
C
      IF ( ICOMP.GT.0 ) THEN
         CALL MB03VY( N, P, ILO, IHI, TMP, LDZ1, LDZ2, TMP, LDTAU,
     $                TMP(1), -1, INFO )
         LDWORK = INT( TMP(1) )
      ELSE
         LDWORK = N
      END IF
      IF ( IJOB.GT.0 )
     $   LDWORK = MAX( LDWORK, 1, IHI - ILO + P - 1 )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      A     = mxCalloc( LDA1*LDA2*P, 8 )
      DWORK = mxCalloc(      LDWORK, 8 )
      TAU   = mxCalloc(     LDTAU*P, 8 )
      WI    = mxCalloc(           N, 8 )
      WR    = mxCalloc(           N, 8 )
      Z     = mxCalloc( LDZ1*LDZ2*P, 8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), %VAL( A ),
     $                       N*N*P )
C
C Do the actual computations.
C
      CALL MB03VD( N, P, ILO, IHI, %VAL( A ), LDA1, LDA2, %VAL( TAU ),
     $             LDTAU, %VAL( DWORK ), INFO )
C
      IF ( IJOB.EQ.0 .AND. ICOMP.EQ.0 ) THEN
         IF ( N.GT.2 )
     $      CALL MA02LD( 'Lower', N-2, N-2, 3, 1, ZERO, ZERO, %VAL( A ),
     $                   LDA1 )
         IF ( N.GT.1 ) THEN
C
            DO 10 I = 2, P
               CALL MA02LD( 'Lower', N-1, N-1, 2+N*N*(I-1), 1, ZERO,
     $                      ZERO, %VAL( A ), LDA1 )
   10       CONTINUE
C
         END IF
      END IF
C
      IF ( ICOMP.GT.0 ) THEN
         CALL DLACPY( 'Full', N, N*P, %VAL( A ), LDA1, %VAL( Z ), LDZ1 )
         IF ( N.GT.2 )
     $      CALL MA02LD( 'Lower', N-2, N-2, 3, 1, ZERO, ZERO, %VAL( A ),
     $                   LDA1 )
         IF ( N.GT.1 ) THEN
C
            DO 20 I = 2, P
               CALL MA02LD( 'Lower', N-1, N-1, 2+N*N*(I-1), 1, ZERO,
     $                      ZERO, %VAL( A ), LDA1 )
   20       CONTINUE
C
         END IF
         CALL MB03VY( N, P, ILO, IHI, %VAL( Z ), LDZ1, LDZ2,
     $                %VAL( TAU ), LDTAU, %VAL( DWORK ), LDWORK, INFO )
      END IF
C
      IF ( IJOB.GT.0 )
     $   CALL MB03WD( JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ, %VAL( A ),
     $                LDA1, LDA2, %VAL( Z ), LDZ1, LDZ2, %VAL( WR ),
     $                %VAL( WI ), %VAL( DWORK ), LDWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      IF ( IJOB.NE.1 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( N, P*N, 0 )
         CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS(1) ),
     $                          N*N*P )
      END IF
      IF ( IJOB.GT.0 ) THEN
         IF ( NLHS.GE.2 ) THEN
            PLHS(IJOB)   = mxCreateDoubleMatrix( N, MIN( N, 1 ), 0 )
            CALL mxCopyReal8ToPtr( %VAL( WR ),
     $                             mxGetPr( PLHS(IJOB) ), N )
            PLHS(IJOB+1) = mxCreateDoubleMatrix( N, MIN( N, 1 ), 0 )
            CALL mxCopyReal8ToPtr( %VAL( WI ),
     $                             mxGetPr( PLHS(IJOB+1) ), N )
         END IF
      END IF
      IF ( ICOMP.GT.0 .AND. NLHS.GE.IJOB+2 ) THEN
         PLHS(IJOB+2) = mxCreateDoubleMatrix( N, P*N, 0 )
         CALL mxCopyReal8ToPtr( %VAL( Z ),
     $                          mxGetPr( PLHS(IJOB+2) ), N*N*P )
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      CALL mxFree( A     )
      CALL mxFree( DWORK )
      CALL mxFree( TAU   )
      CALL mxFree( WI    )
      CALL mxFree( WR    )
      CALL mxFree( Z     )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB03WD'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of PERSCHUR ***
      END
