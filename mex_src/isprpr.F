#include "fintrf.h"
C ISPRPR.F - Gateway function to check whether the transfer function of
C            a descriptor system with given generalized state space
C            realization is proper or not, using the SLICOT Library
C            routine AB13ID.
C
C RELEASE 2.0 of SLICOT Basic Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [prop(,Aout,Eout,Bout,Cout,warn)] =
C         isprpr(A,E,B,C(,jobsys,jobeig,equil,cksing,restor,update,tol))
C   [prop(,warn)] =
C         isprpr(A,E,B,C(,jobsys,0,equil,cksing,restor,0,tol))
C
C Purpose:
C   To check whether the transfer function
C
C                                   -1
C     G(lambda) := C*( lambda*E - A ) *B
C
C   of a given linear time-invariant descriptor system with generalized
C   state space realization (lambda*E-A,B,C) is proper.
C   Optionally, the system (lambda*E-A,B,C) is reduced to an equivalent
C   one (lambda*Er-Ar,Br,Cr) with only controllable and observable
C   eigenvalues in order to use it for a subsequent L_inf-norm
C   computation; or the system is reduced to an equivalent one
C   (lambda*Er-Ar,Br,Cr) without uncontrollable and unobservable
C   infinite eigenvalues. In this case, if update = 0, the system is
C   only checked for properness, the reduced system is not returned.
C
C Input parameters:
C   A       - the n-by-n system state matrix A.
C   E       - the n-by-n system descriptor matrix E.
C   B       - the n-by-m system input matrix B.
C   C       - the p-by-n system output matrix C.
C   jobsys  - (optional) specifies if the system is already in the
C             reduced form as stated in jobeig
C             = 0 : system is not in reduced form;
C             = 1 : system is in reduced form.
C             Default: jobsys = 0.
C   jobeig  - (optional) specifies which kind of eigenvalues of the
C             system pencil lambda*E-A should be removed if jobsys = 0.
C             = 0 : remove only infinite uncontrollable and unobservable
C                   eigenvalues;
C             = 1 : remove all uncontrollable and unobservable
C                   eigenvalues.
C             Default: jobeig = 0.
C   equil   - (optional) specifies whether the user wishes to
C             preliminarily equilibrate the system (lambda*E-A,B,C):
C             = 0 : do not perform equilibration;
C             = 1 : perform equilibration (scaling).
C             Default: equil = 0.
C   cksing  - (optional) specifies whether the user wishes to check if
C             the pencil (A-lambda*E) is singular:
C             = 0 : do not check singularity;
C             = 1 : check singularity.
C             If the pencil is singular, the reduced system computed for
C             cksing = 0 can have completely different eigenvalues than
C             the given system.
C             The test is performed only if jobsys = 0.
C             Default: cksing = 1.
C   restor  - (optional) specifies whether the user wishes to save the
C             system matrices before each reduction phase if jobsys = 0,
C             and restore them if no order reduction took place:
C             = 0 : do not save the matrices;
C             = 1 : save and restore.
C             While saving and restoring can sometimes give more
C             accurate reduced systems (with poles/zeros closer to the
C             original ones), this is not always true. Use restor = 0
C             for testing the properness.
C             Default: restor = 0.
C   update  - (optional) specifies whether the user wishes to update the
C             matrices A, B, and C if jobeig = 0:
C             = 0 : do not update the matrices A, B and C;
C             = 1 : update the matrices A, B and C.
C             Default: update = 0.
C   tol     - (optional) vector of length 3 containing the tolerance
C             used to set the accuracy in determining ranks (in tol(1)),
C             the tolerance used for checking pencil singularity,
C             when cksing = 1, and/or the singularity of the matrices A
C             and E, when cksing = 0 (in tol(2)), and the threshold
C             value for magnitude of the matrix elements, if equil = 1:
C             elements with magnitude less than or equal to tol(3) are
C             ignored for scaling.
C             Default: tol(1) = n*n*(epsilon_machine),
C                      tol(2) =  10*(epsilon_machine),
C                      tol(3) = max(norm_1(A,E,B,C))*(epsilon_machine),
C             where epsilon_machine is the relative machine precision.
C
C Output parameters:
C   prop    - indicates whether the transfer function of the system is
C             proper or not.
C             = 0 : the transfer function is improper;
C             = 1 : the transfer function is proper.
C   Aout    - (optional) if jobeig = 1 or update = 1, the nr-by-nr
C             reduced or transformed system state matrix Aout.
C   Eout    - (optional) if jobeig = 1 or update = 1, the nr-by-nr
C             reduced or transformed system descriptor matrix Eout.
C   Bout    - (optional) if jobeig = 1 or update = 1, the nr-by-m
C             reduced or transformed system input matrix Bout.
C   Cout    - (optional) if jobeig = 1 or update = 1, the p-by-nr
C             reduced or transformed system output matrix Cout.
C   warn    - (optional) warning indicator which indicates the quality
C             of rank decisions
C             = 0 : rank decisions can be assumed to be correct;
C             = 1 : the tolerance tol(1) is very close to a singular
C                   value of a matrix whose rank should be computed;
C                   rank decisions might be incorrect.
C
C Contributor:
C   M. Voigt, Chemnitz University of Technology, Department of
C   Mathematics, Germany, Apr. 2010.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Romania,
C   Jan. 2012, Apr. 2012, Jul. 2016.
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
      mwPointer         PLHS( * ), PRHS( * )
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      mwPointer         mxCalloc, mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsComplex, mxIsNumeric
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         CKSING, EQUIL, JOBEIG, JOBSYS, RESTOR, UPDATE
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDE, LDWORK, M,
     $                  N, NR, P, RANKE
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: DW(:)
      mwPointer         A, B, C, DWORK, E
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           LISPRP
      CHARACTER*120     TEXT
      INTEGER           ICK, IDUM( 1 ), IEIG, IP, IRST, ISYS, ITMP,
     $                  IUPD, LDW, LDWC, LIWORK, MAXMP, NC, NR1, NT
      DOUBLE PRECISION  DUM( 1 ), TEMP, TOL( 3 )
C
C .. External functions ..
      LOGICAL           AB13ID
      EXTERNAL          AB13ID
C
C     .. External Subroutines ..
      EXTERNAL          DLACPY
C
C ..Intrinsic functions..
      INTRINSIC         MAX
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ISPRPR requires at least 4 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ISPRPR requires at least 1 output argument' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   A(nxn), E(nxn), B(nxm), C(pxn)
C
      N = mxGetM( PRHS(1) )
      M = mxGetN( PRHS(3) )
      P = mxGetM( PRHS(4) )
C
      IF ( mxGetN( PRHS(1) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a numeric matrix' )
      END IF
C
      IF ( mxGetM( PRHS(2) ).NE.N .OR. mxGetN( PRHS(2) ).NE.N )
     $      THEN
         CALL mexErrMsgTxt
     $       ( 'E must be a square matrix of the same order as A' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'E must be a numeric matrix' )
      END IF
C
      IF ( mxGetM( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same row dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a numeric matrix' )
      END IF
C
      IF ( mxGetN( PRHS(4) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $       ( 'C must have the same column dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a numeric matrix' )
      END IF
C
C   jobsys, jobeig, equil, cksing, restor, update, tol
C
      IP = 5
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'JOBSYS must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'JOBSYS must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ISYS = TEMP
         IF ( ISYS.NE.0 .AND. ISYS.NE.1 ) THEN
              CALL mexErrMsgTxt( 'JOBSYS must be 0 or 1' )
         END IF
         IF ( ISYS.EQ.0 ) THEN
            JOBSYS = 'R'
         ELSE
            JOBSYS = 'N'
         END IF
         IP = IP + 1
      ELSE
         ISYS   = 0
         JOBSYS = 'R'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'JOBEIG must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'JOBEIG must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IEIG = TEMP
         IF ( IEIG.NE.0 .AND. IEIG.NE.1 ) THEN
              CALL mexErrMsgTxt( 'JOBEIG must be 0 or 1' )
         END IF
         IF ( IEIG.EQ.0 ) THEN
            JOBEIG = 'I'
         ELSE
            JOBEIG = 'A'
         END IF
         IP = IP + 1
      ELSE
         IEIG   = 0
         JOBEIG = 'I'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'EQUIL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'EQUIL must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ITMP = TEMP
         IF ( ITMP.NE.0 .AND. ITMP.NE.1 ) THEN
              CALL mexErrMsgTxt( 'EQUIL must be 0 or 1' )
         END IF
         IF ( ITMP.EQ.1 ) THEN
            EQUIL = 'S'
         ELSE
            EQUIL = 'N'
         END IF
         IP = IP + 1
      ELSE
         ITMP  = 0
         EQUIL = 'N'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'CKSING must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'CKSING must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ICK = TEMP
         IF ( ICK.NE.0 .AND. ICK.NE.1 ) THEN
              CALL mexErrMsgTxt( 'CKSING must be 0 or 1' )
         END IF
         IF ( ICK.EQ.0 ) THEN
            CKSING = 'N'
         ELSE
            CKSING = 'C'
         END IF
         IP = IP + 1
      ELSE IF ( ISYS.EQ.0 ) THEN
         ICK    = 1
         CKSING = 'C'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'RESTOR must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'RESTOR must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IRST = TEMP
         IF ( IRST.NE.0 .AND. IRST.NE.1 ) THEN
              CALL mexErrMsgTxt( 'RESTOR must be 0 or 1' )
         END IF
         IF ( IRST.EQ.0 ) THEN
            RESTOR = 'N'
         ELSE
            RESTOR = 'R'
         END IF
         IP = IP + 1
      ELSE
         IRST   = 0
         RESTOR = 'N'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'UPDATE must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'UPDATE must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IUPD = TEMP
         IF ( IUPD.NE.0 .AND. IUPD.NE.1 ) THEN
              CALL mexErrMsgTxt( 'UPDATE must be 0 or 1' )
         END IF
         IF ( IUPD.EQ.0 ) THEN
            UPDATE = 'N'
         ELSE
            UPDATE = 'U'
         END IF
         IP = IP + 1
      ELSE
         IUPD   = 0
         UPDATE = 'N'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         NT = mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) )
         IF ( NT.GT.3 ) THEN
              CALL mexErrMsgTxt
     $             ( 'TOL must be a vector of length at most 3' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'TOL must be a numeric vector' )
         END IF
         IF ( NT.GT.0 ) THEN
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TOL, NT )
            IF ( TOL( 1 ).GE.ONE ) THEN
               CALL mexErrMsgTxt( 'TOL(1) must be less than one' )
            ELSE IF ( NT.EQ.2 ) THEN
               IF ( TOL( 2 ).GE.ONE )
     $            CALL mexErrMsgTxt
     $                 ( 'TOL(2) must be less than one' )
               TOL( 3 ) = -ONE
            ELSE IF ( NT.EQ.3 ) THEN
               IF ( TOL( 3 ).GE.ONE )
     $            CALL mexErrMsgTxt
     $                 ( 'TOL(3) must be less than one' )
            ELSE
               TOL( 2 ) = ZERO
               TOL( 3 ) = -ONE
            END IF
         ELSE
            TOL( 1 ) = ZERO
            TOL( 2 ) = ZERO
            TOL( 3 ) = -ONE
         END IF
      ELSE
         TOL( 1 ) = ZERO
         TOL( 2 ) = ZERO
         TOL( 3 ) = -ONE
      END IF
C
C Determine the lengths of working arrays.
C Use optimal value for LDWORK.
C
      MAXMP = MAX( M, P )
      LDA   = MAX( 1, N )
      LDE   = LDA
      LDB   = LDA
      LDWC  = MAX( 1, P )
      IF ( N.GT.0 ) THEN
         LDC = MAX( 1, MAXMP )
      ELSE
         LDC = 1
      END IF
      IF ( IEIG.EQ.1 .OR. IUPD.EQ.1 ) THEN
         LDW = MAX( LDA, LDWC )
         NC  = MAX( N, M )
      ELSE
         LDW = LDWC
         NC  = N
      END IF
      IF ( ISYS.EQ.0 ) THEN
         LIWORK = 2*N + MAXMP + 7
      ELSE
         LIWORK = N
      END IF
C
      LISPRP = AB13ID( JOBSYS, JOBEIG, EQUIL, CKSING, RESTOR, UPDATE, N,
     $                 M, P, DUM, LDA, DUM, LDE, DUM, LDB, DUM, LDC, NR,
     $                 RANKE, TOL, IDUM, DUM, -1, IWARN, INFO )
      LDWORK = DUM( 1 )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( DW( LDW*NC ), IWORK( LIWORK ) )
      A     = mxCalloc( LDA*N,     8 )
      B     = mxCalloc( LDB*MAXMP, 8 )
      C     = mxCalloc( LDC*N,     8 )
      DWORK = mxCalloc( LDWORK,    8 )
      E     = mxCalloc( LDE*N,     8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), %VAL( A ),   N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( E ),   N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL( B ),   N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), DW, P*N )
      CALL DLACPY( 'Full', P, N, DW, LDWC, %VAL( C ), LDC )
C
C Do the actual computations.
C
      LISPRP = AB13ID( JOBSYS, JOBEIG, EQUIL, CKSING, RESTOR, UPDATE, N,
     $                 M, P, %VAL( A ), LDA, %VAL( E ), LDE, %VAL( B ),
     $                 LDB, %VAL( C ), LDC, NR, RANKE, TOL, IWORK,
     $                 %VAL( DWORK ), LDWORK, IWARN, INFO )
C
C Copy output to MATLAB workspace.
C The original system is returned if the pencil lambda*E-A is singular.
C
      IF ( INFO.EQ.0 .OR. INFO.EQ.1 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( 1, 1, 0 )
         IF ( LISPRP ) THEN
            CALL mxCopyReal8ToPtr( ONE, mxGetPr( PLHS(1) ), 1 )
         ELSE
            CALL mxCopyReal8ToPtr( ZERO, mxGetPr( PLHS(1) ), 1 )
         END IF
C
         IP = 2
         IF ( NLHS.GE.IP .AND. ( IEIG.EQ.1 .OR. IUPD.EQ.1 ) ) THEN
            NR1 = MAX( 1, NR )
            CALL DLACPY( 'Full', NR, NR, %VAL( A ), LDA, DW, NR1 )
            PLHS(IP) = mxCreateDoubleMatrix( NR, NR, 0 )
            CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(IP) ),
     $                             NR*NR )
            IP = IP + 1
C
            CALL DLACPY( 'Full', NR, NR, %VAL( E ), LDE, DW, NR1 )
            PLHS(IP) = mxCreateDoubleMatrix( NR, NR, 0 )
            CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(IP) ),
     $                             NR*NR )
            IP = IP + 1
C
            CALL DLACPY( 'Full', NR, M, %VAL( B ), LDB, DW, NR1 )
            PLHS(IP) = mxCreateDoubleMatrix( NR, M, 0 )
            CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(IP) ),
     $                             NR*M )
            IP = IP + 1
C
            CALL DLACPY( 'Full', P, NR, %VAL( C ), LDC, DW, LDWC )
            PLHS(IP) = mxCreateDoubleMatrix( P, NR, 0 )
            CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(IP) ),
     $                             P*NR )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
            IF ( IWARN.EQ.1 ) THEN
               CALL mxCopyReal8ToPtr( ONE, mxGetPr( PLHS(IP) ), 1 )
            ELSE
               CALL mxCopyReal8ToPtr( ZERO, mxGetPr( PLHS(IP) ), 1 )
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( DW, IWORK )
      CALL mxFree( A     )
      CALL mxFree( B     )
      CALL mxFree( C     )
      CALL mxFree( DWORK )
      CALL mxFree( E     )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' The pencil A - lambda*E is singular.'',
     $                  '' The original system is returned.'')' )
         CALL mexPrintf( TEXT )
      END IF
C
      RETURN
C *** Last line of ISPRPR ***
      END
 
