C HAMILEIG.F - Gateway function to compute the eigenvalues of a
C              Hamiltonian matrix using the square-reduced approach
C              implemented in SLICOT routines MB03SD and MB04ZD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [Ao,QGo(,U)]       = Hamileig(A,QG,job(,compu,S)),        job = -1;
C   [WR,WI(,Ao,QGo,U)] = Hamileig(A,QG(,job,jobscl,compu,S)), job =  0,
C                                                          or job =  1.
C
C   [Ao,QGo(,U)]       = Hamileig(A,QG,-1(,compu,S))
C   [WR,WI(,Ao,QGo,U)] = Hamileig(A,QG(,0, jobscl,compu,S))
C   [WR,WI]            = Hamileig(A,QG, 1(,jobscl))
C
C Purpose:
C   To transform a Hamiltonian matrix
C
C             ( A   G  )
C         H = (      T )                                           (1)
C             ( Q  -A  )
C
C   (with G and Q symmetric matrices) into a square-reduced Hamiltonian
C   matrix
C
C              ( A'   G'  )
C         H' = (        T ),                                       (2)
C              ( Q'  -A'  )
C                                                               T
C   by an orthogonal symplectic similarity transformation H' = U H U,
C   where
C
C               (  U1   U2 )
C           U = (          ),                                      (3)
C               ( -U2   U1 )
C
C   and to compute the eigenvalues of H.  Therefore, H' is such that
C
C         2    ( A''   G'' )
C       H'  =  (         T ),                                      (4)
C              ( 0    A''  )
C
C   with A'' upper Hessenberg and G'' skew symmetric.  The square roots
C   of the eigenvalues of A'' = A'*A' + G'*Q' are the eigenvalues of H.
C
C Input parameters:
C   A      - the n-by-n matrix A.
C   QG     - an  n-by-(n+1) matrix containing the triangles of the
C            symmetric matrices Q and G, as follows:
C            the leading n-by-n lower triangular part contains the lower
C            triangle of the matrix Q, and the n-by-n upper triangular
C            part of the submatrix in the columns 2 to n+1 contains the
C            upper triangle of the matrix G of H in (1).
C            So, if i >= j, then Q(i,j) = Q(j,i) is stored in QG(i,j)
C            and G(i,j) = G(j,i) is stored in QG(j,i+1).
C            QG is an empty matrix if n = 0.
C   job    - (optional) scalar indicating the computation to be
C            performed, as follows:
C            =-1 :  compute the square-reduced matrix H' in (2);
C            = 0 :  compute the eigenvalues of H (default);
C            = 1 :  compute the eigenvalues of H, assuming that the
C                   given Hamiltonian matrix is already in the reduced
C                   form (2).
C   jobscl - (optional) if job >= 0, scalar specifying whether or not
C            balancing operations should be performed when computing the
C            eigenvalues of A'', as follows:
C            = 0 :  do not use balancing;
C            = 1 :  do scaling in order to equilibrate the rows
C                   and columns of A'' (default).
C            If job = -1, jobscl is not used.
C   compu  - (optional) if job <= 0, scalar indicating whether the
C            orthogonal symplectic similarity transformation matrix U is
C            returned or accumulated into an orthogonal symplectic
C            matrix, or if the transformation matrix is not required,
C            as follows:
C            = 0 :  U is not required (default);
C            = 1 :  on entry, U need not be set;
C                   on exit, U contains the orthogonal symplectic
C                   matrix U;
C            = 2 :  the orthogonal symplectic similarity transformations
C                   are accumulated into U;
C                   on input, U must contain an orthogonal symplectic
C                   matrix S;
C                   on exit, U contains S*U.
C            If job = 1, compu is not used.
C   S      - (optional) if job <= 0 and compu = 2, an n-by-2*n matrix
C            containing the first n rows of the given orthogonal
C            symplectic matrix S.
C
C Output parameters:
C   WR,WI  - if job >= 0, the n-vectors of real parts and imaginary
C            parts, respectively, of the n computed eigenvalues of H'
C            with non-negative real part. The remaining n eigenvalues
C            are the negatives of these eigenvalues.
C            Eigenvalues are stored in WR and WI in decreasing order of
C            magnitude of the real parts, i.e., WR(I) >= WR(I+1).
C            (In particular, an eigenvalue closest to the imaginary
C             axis is WR(N)+WI(N)i.)
C            In addition, eigenvalues with zero real part are sorted in
C            decreasing order of magnitude of imaginary parts.  Note
C            that non-real eigenvalues with non-zero real part appear
C            in complex conjugate pairs, but eigenvalues with zero real
C            part do not, in general, appear in complex conjugate
C            pairs.
C   Ao     - if job <= 0, the computed n-by-n submatrix A' of H' in (2).
C   QGo    - if job <= 0, the computed n-by-(n+1) matrix containing the
C            triangles of the symmetric matrices Q' and G' of H' in (2),
C            stored in the same way as the initial matrices Q and G.
C   U      - if job <= 0 and compu > 0, an n-by-2*n matrix containing
C            the first n rows of the computed orthogonal symplectic
C            matrix U.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2002.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009,
C   Dec. 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Mex-file interface parameters ..
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         COMPU, JOBSCL
      INTEGER           INFO, LDA, LDQG, LDU, LDWORK, N
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), DWORK(:), QG(:,:),
     $                                 U(:,:), WI(:), WR(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           ICOMP, IJOBS, IP, ITMP, JOB
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          MB03SD, MB04ZD
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HAMILEIG requires at least 2 input arguments' )
      ELSE IF ( NLHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HAMILEIG requires at least 2 output arguments' )
      END IF
C
C   A(nxn), QG(nx(n+1)) (, job(, jobscl), compu, S).
C
      N = mxGetM( PRHS(1) )
      IF ( mxGetN( PRHS(1) ).NE.N ) THEN
         WRITE( TEXT, '(''A must have '',I5,'' rows and columns'')' ) N
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(2) ).NE.N ) THEN
         WRITE( TEXT, '(''QG must have '',I5,'' rows'')' ) N
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( mxGetN( PRHS(2) ).NE.N + MIN( N, 1 ) ) THEN
         WRITE( TEXT, '(''QG must have '',I5,'' columns'')' )
     $          N + MIN( N, 1 )
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'QG must be a real matrix' )
      END IF
C
C   job
C
      IP = 3
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
         JOB = TEMP
         IF ( JOB.LT.-1 .OR. JOB.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'JOB has -1, 0, or 1 the only admissible values' )
         END IF
         IP = IP + 1
      ELSE
         JOB = 0
      END IF
C
C   jobscl
C
      IF ( JOB.GE.0 ) THEN
         IF ( NRHS.GE.IP ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'JOBSCL must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'JOBSCL must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            IJOBS = TEMP
            IF ( IJOBS.LT.0 .OR. IJOBS.GT.1 ) THEN
               CALL mexErrMsgTxt
     $           ( 'JOBSCL has 0 or 1 the only admissible values' )
            END IF
            IP = IP + 1
         ELSE
            IJOBS = 1
         END IF
C
         IF ( IJOBS.EQ.0 ) THEN
            JOBSCL = 'N'
         ELSE
            JOBSCL = 'S'
         END IF
      END IF
C
C   compu
C
      IF ( JOB.LE.0 ) THEN
         IF ( NRHS.GE.IP ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'COMPU must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'COMPU must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            ICOMP = TEMP
            IF ( ICOMP.LT.0 .OR. ICOMP.GT.2 ) THEN
               CALL mexErrMsgTxt
     $              ( 'COMPU has 0, 1 or 2 the only admissible values' )
            END IF
            IP = IP + 1
         ELSE
            ICOMP = 0
         END IF
C
         IF ( ICOMP.EQ.0 ) THEN
            COMPU = 'N'
         ELSE IF ( ICOMP.EQ.1 ) THEN
            COMPU = 'I'
         ELSE
            COMPU = 'V'
         END IF
      END IF
C
      IF ( JOB.LE.0 .AND .ICOMP.EQ.2 ) THEN
         IF ( NRHS.LT.IP ) THEN
            WRITE( TEXT, '(''HAMILEIG requires '',I2,
     $                     '' input arguments'')' ) IP
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxGetM( PRHS(IP) ).NE.N ) THEN
            WRITE( TEXT, '(''S must have '', I5, '' rows'')' ) N
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxGetN( PRHS(IP) ).NE.2*N ) THEN
            WRITE( TEXT, '(''S must have '', I5, '' columns'')' ) 2*N
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'S must be a real matrix' )
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA  = MAX( 1, N )
      LDQG = LDA
      IF ( JOB.LE.0 .AND. ICOMP.GT.0 ) THEN
         LDU  = LDA
         ITMP = 2*N
      ELSE
         LDU  = 1
         ITMP = 1
      END IF
C
C   ldwork
C
      IF ( JOB.LT.0 ) THEN
         LDWORK = 2*N
      ELSE
         LDWORK = MAX( 1, N*( N + 1 ) )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( A( LDA, N ), DWORK( LDWORK ), QG( LDQG, N+1 ) )
      IF ( JOB.LE.0 )
     $   ALLOCATE ( U( LDU, ITMP ) )
      IF ( JOB.GE.0 )
     $   ALLOCATE ( WI( N ), WR( N ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), A,  N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), QG, N*( N + 1 ) )
      IF ( JOB.LE.0 .AND. ICOMP.EQ.2 )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), U, 2*N*N )
C
C Do the actual computations.
C
      IF ( JOB.LE.0 )
     $   CALL MB04ZD( COMPU, N, A, LDA, QG, LDQG, U, LDU, DWORK, INFO )
C
      IF ( JOB.GE.0 )
     $   CALL MB03SD( JOBSCL, N, A, LDA, QG, LDQG, WR, WI, DWORK,
     $                LDWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      IF ( JOB.LT.0 ) THEN
         IP = 1
      ELSE
         IP = 3
      END IF
C
      IF ( NLHS.GE.IP ) THEN
         PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
         CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(IP) ), N*N )
         IP = IP + 1
         IF ( NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N + MIN( N, 1 ), 0 )
            CALL mxCopyReal8ToPtr( QG, mxGetPr( PLHS(IP) ),
     $                             N*( N + MIN( N, 1 ) ) )
            IP = IP + 1
         END IF
         IF ( ICOMP.GT.0 .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, 2*N, 0 )
            CALL mxCopyReal8ToPtr( U, mxGetPr( PLHS(IP) ), 2*N*N )
         END IF
      END IF
C
      IF ( JOB.GE.0 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( N, MIN( N, 1 ), 0 )
         CALL mxCopyReal8ToPtr( WR, mxGetPr( PLHS(1) ), N )
         PLHS(2) = mxCreateDoubleMatrix( N, MIN( N, 1 ), 0 )
         CALL mxCopyReal8ToPtr( WI, mxGetPr( PLHS(2) ), N )
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( A, DWORK, QG )
      IF ( JOB.LE.0 )
     $   DEALLOCATE ( U )
      IF ( JOB.GE.0 )
     $   DEALLOCATE ( WI, WR )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM MB03SD'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of HAMILEIG ***
      END
