#include "fintrf.h"
C HAPACK_HAEIG.F
C   Gateway function for Hamiltonian eigenvalue computations.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C                    [...] = hapack_haeig(task,...)
C
C   PVL decomposition
C   task = 1 :    [Ar,QGr] = hapack_haeig(1,A,QG)
C           [U1,U2,Ar,QGr] = hapack_haeig(1,A,QG)
C
C   URV decomposition
C   task = 2 :     [S,T,G] = hapack_haeig(2,A,QG)
C            [U1,U2,S,T,G] = hapack_haeig(2,A,QG)
C      [U1,U2,V1,V2,S,T,G] = hapack_haeig(2,A,QG)
C                     [Br] = hapack_haeig(2,B)
C               [U1,U2,Br] = hapack_haeig(2,B)
C         [U1,U2,V1,V2,Br] = hapack_haeig(2,B)
C
C   URV/PS decomposition
C   task = 3 : [S,T,wr,wi] = hapack_haeig(3,A,QG)
C            [S,T,G,wr,wi] = hapack_haeig(3,A,QG)
C        [U1,U2,S,T,wr,wi] = hapack_haeig(3,A,QG)
C      [U1,U2,S,T,G,wr,wi] = hapack_haeig(3,A,QG)
C        [U1,U2,V1,V2,...
C               S,T,wr,wi] = hapack_haeig(3,A,QG)
C       [U1,U2,V1,V2,...
C             S,T,G,wr,wi] = hapack_haeig(3,A,QG)
C
C   Balancing
C   task = 4 : [scale,ilo] = hapack_haeig(4,A,QG,balanc)
C       [Ar,QGr,scale,ilo] = hapack_haeig(4,A,QG,balanc)
C
C   Simply eigenvalues:
C   task = 5:      [wr,wi] = hapack_haeig(5,A,QG,balanc)
C
C   Reordering Hamiltonian Schur form:
C   task = 6 :    [Ar,QGr] = hapack_haeig(6,A,QG,select,lower)
C           [U1,U2,Ar,QGr] = hapack_haeig(6,A,QG,select,lower)
C           [U1,U2,Ar,QGr] = hapack_haeig(6,U1,U2,A,QG,select,lower)
C
C   Selected stable/unstable invariant subspaces
C   task = 7 :        [US] = hapack_haeig(7,S,T,U1,U2,V1,V2,select)
C                  [US,UU] = hapack_haeig(7,S,T,U1,U2,V1,V2,select)
C
C   task = 8 :        [UU] = hapack_haeig(8,S,T,U1,U2,V1,V2,select)
C
C   Complete stable/unstable invariant subspaces:
C   task = 9 :        [US] = hapack_haeig(9,S,T,U1,U2,V1,V2)
C                     [US] = hapack_haeig(9,S,T,G,U1,U2,V1,V2)
C                  [US,UU] = hapack_haeig(9,S,T,U1,U2,V1,V2)
C                  [US,UU] = hapack_haeig(9,S,T,G,U1,U2,V1,V2)
C
C   task = 10 :       [UU] = hapack_haeig(10,S,T,U1,U2,V1,V2)
C                     [UU] = hapack_haeig(10,S,T,G,U1,U2,V1,V2)
C Purpose:
C   To compute eigenvalues and invariant subspaces of a 2*N-by-2*N
C   Hamiltoninan matrix H = [A, G; Q, -A']. On input as well as on
C   output Hamiltonian matrices are stored in compressed format in
C   an N-by-N array Ar and an N-by-(N+1) array QGr.
C
C   task =  1: Computes the PVL form Hr = [Ar, Gr; Qr, -Ar'] of a
C              Hamiltonian H with diagonal Qr and Ar in upper Hessenberg
C              form. Optionally, N-by-N matrices U1 and U2 are computed
C              so that U = [U1, U2; -U2, U1] satisfies U'*U = I and
C              H = U*Hr*U'.
C
C   task =  2: Computes the symplectic URV form Br = [T, G; 0, S'] of a
C              Hamiltonian H or a general matrix B so that T is upper
C              triangular and S is in upper Hessenberg form. Optionally,
C              N-by-N matrices U1, U2, V1 and V2 are computed so that
C              U = [U1, U2; -U2, U1], V = [V1', V2'; -V2', V1'] (or
C              V = [V1, V2'; -V2', V1] if B is Hamiltonian and stored in
C              compressed format) satisfy U'*U = I, V'*V = I and
C              B = U*Br*V'.
C
C   task =  3: Computes the symplectic URV/periodic Schur form
C              Hr = [T, G; 0, S'] of a Hamiltonian matrix H so that T is
C              upper triangular and S is in real Schur form. Optionally,
C              N-by-N matrices U1, U2, V1 and V2 are computed so that
C              U = [U1, U2; -U2, U1], V = [V1, V2; -V2, V1] satisfy
C              U'*U = I, V'*V = I and H = U*Hr*V'. Also, the nonpositive
C              square roots of the eigenvalues of -S*T, which are equal
C              to the stable eigenvalues of H, are computed and returned
C              in vectors wr and wi containing the real and imaginary
C              parts, respectively.
C
C   task =  4: Finds an accurate similarity transformation T such that
C              Hr = T\H*T has, as nearly as possible, approximately
C              equal row and column norms. T is a permutation of a
C              diagonal matrix and symplectic. T is stored in an
C              n-vector SCALE as described in MB04DD.
C
C   task =  5: Computes the stable eigenvalues of H and returns vectors
C              wr and wi containing the real and imaginary parts,
C              respectively.
C
C   task =  6: Reorders the eigenvalues selected by the arrays SELECT
C              and LOWER to the top left part of a matrix H in
C              Hamiltonian Schur form. Optionally, the corresponding
C              orthogonal symplectic transformation matrix
C              U = [U1, U2; -U2, U1] is computed or a given matrix
C              U = [U1, U2; -U2, U1] is post-multiplied by this matrix.
C
C   task =  7
C   task =  8: Computes the invariant subspace US (UU) belonging to
C              a set of selected stable (unstable) eigenvalues of H.
C              The user must provide the matrices U1, U2, V1, V2, S
C              and T as computed by MB03XD (task 3).
C
C   task =  9
C   task = 10: Computes the invariant subspace US (UU) belonging to
C              all stable (unstable) eigenvalues of H. The user must
C              provide the matrices U1, U2, V1, V2, S, T and
C              possibly G as computed by MB03XD (task 3). Note that
C              if G is provided then US and UU are computed from a
C              larger set of vectors which may imply more accurate
C              bases.
C
C Input parameters:
C   task   - integer option to determine the computation to perform as
C            described above.
C   A      - real N-by-N matrix. If task = 6 this matrix must be in real
C            Schur form.
C   B      - real 2*N-by-2*N matrix.
C   G      - real N-by-N matrix.
C   QG     - real N-by-(N+1) matrix containing the symmetric matrices
C            Q and G in compressed format. If task = 6 the part
C            containing the matrix Q must be zero.
C   S      - a real N-by-N matrix in real Schur form.
C   T      - a real, upper triangular N-by-N matrix.
C   U1,U2  - real N-by-N matrices.
C   V1,V2  - real N-by-N matrices.
C   balanc - determines whether H should be permuted (balanc = 'p'),
C            scaled (balanc = 's') or permuted and scaled (balanc = 'b')
C            prior to eigenvalue computations. Otherwise balanc = 'n'.
C   lower  - N-vector controlling whether lambda or -lambda of each
C            selected eigenvalue pair (lambda,-lambda) is to be
C            reordered to the top part of a Hamiltonian Schur form.
C   select - if task = 6: N-vector specifying the eigenvalues to be
C            reordered to the top part of a Hamiltonian Schur form; and
C            if task > 6: N-vector specifying the eigenvalue pairs to
C            which stable/unstable invariant subspaces are to be
C            computed.
C
C Output parameters:
C   Ar     - real N-by-N matrix containing the (1,1) block of the
C            returned Hamiltonian matrix Hr.
C   QGr    - real N-by-(N+1) matrix containing the symmetric (2,1) and
C            (1,2) blocks of the returned (skew-)Hamiltonian matrix Hr
C            in compressed format.
C   Br     - real 2*N-by-2*N matrix in symplectic URV form.
C   T,S,G  - submatrices of a symplectic URV(/periodic Schur) form of H.
C   U1,U2  - real N-by-N matrices so that U = [U1, U2; -U2, U1]
C            satisfies U'*U = I.
C   V1,V2  - real N-by-N matrices so that V = [V1', V2'; -V2', V1'] (if
C            task = 2) or V = [V1, V2; -V2, V1] (otherwise) satisfies
C            V'*V = I.
C   US     - real 2*N-by-M matrix containing the orthonormal basis of a
C            stable invariant subspace.
C   UU     - real 2*N-by-M matrix containing the orthonormal basis of an
C            unstable invariant subspace.
C   ilo    - ilo-1 is the number of deflated eigenvalues in the balanced
C            Hamiltonian matrix.
C   scale  - real N-vector containing the scaling factors as returned by
C            MB04DD.
C   wr     - real N-vector containing the real parts of the eigenvalues
C            of H.
C   wi     - real N-vector containing the imaginary parts of the
C            eigenvalues of H.
C
C Contributor:
C   D. Kressner, UCL CESAME, Louvain-la-Neuve, June 2003.
C
C Revisions:
C   V. Sima, June 2008 (SLICOT version), Apr. 2009, Aug. 2011,
C   Oct. 2011, Sep. 2012, Aug. 2015.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C
C     .. Mex-file interface parameters ..
      mwPointer         PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C     .. Mex-file integer functions ..
      mwPointer         mxCalloc, mxCreateDoubleMatrix,
     $                  mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxGetString, mxIsChar,
     $                  mxIsComplex, mxIsNumeric
C
C     .. Scalar parameters used by HAPACK subroutines ..
      CHARACTER         BALANC, JOB, JOBU, JOBV, METH, TRAN
      INTEGER           ILO, INFO, LDA, LDB, LDG, LDQ, LDQG, LDS, LDT,
     $                  LDWORK, LDU1, LDU2, LDUS, LDUU, LDV1, LDV2, M,
     $                  N, USCOL
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: LOWER(:), LWORK(:), SELECT(:)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: DSEL(:), TW(:,:)
      mwPointer         A, B, B11, B12, B21, B22, CSL, CSR, DW, DWORK,
     $                  G, Q, QG, SCALE, T, TAUL, TAUR, U1, U2, US, UU,
     $                  V1, V2, WI, WR
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           COMPG, COMPQ, COMPT, COMPU, COMPUS, COMPUU,
     $                  COMPV, LW(1), PAIR, READA, READB, READG, READQG,
     $                  READT, UPDU, UPDV
      INTEGER           I, IERR, IP, IW(1), J, JG, N1, NLM, NN, TASK
      DOUBLE PRECISION  DM(1), TEMP
C
C     .. External Functions ..
      INTEGER           UE01MD
      EXTERNAL          UE01MD
C
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DLACPY, DLASCL, DLASET, MA02ED, MA02KD,
     $                  MA02KV, MA02LD, MB03TD, MB03XD, MB03ZD, MB04DD,
     $                  MB04PB, MB04TB, MB04WP, MB04WR
C
C     ..Intrinsic Functions..
      INTRINSIC         DBLE, INT, MAX, MIN
C
C     Check for proper number of arguments.
C
      IF ( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HAPACK_HAEIG requires at least 2 input arguments' )
      ELSE IF ( NLHS.GT.9 ) THEN
         CALL mexErrMsgTxt
     $        ( 'HAPACK_HAEIG requires at most 9 output arguments' )
      END IF
C
C     Check dimensions of input parameters and read/set scalar
C     parameters.
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 )
     $   CALL mexErrMsgTxt( 'TASK must be a scalar' )
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 )
     $   CALL mexErrMsgTxt( 'TASK must be an integer scalar' )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      TASK = TEMP
C
      IF ( TASK.LT.1 .OR. TASK.GT.10 )
     $   CALL mexErrMsgTxt
     $        ( 'The only admissible values of TASK are 1, 2, ..., 10')
C
      COMPG  = .FALSE.
      COMPQ  = .FALSE.
      COMPT  = .FALSE.
      COMPU  = .FALSE.
      COMPUS = .FALSE.
      COMPUU = .FALSE.
      COMPV  = .FALSE.
      UPDU   = .FALSE.
      UPDV   = .FALSE.
      READA  = .TRUE.
      READB  = .FALSE.
      READG  = .FALSE.
      READQG = .TRUE.
      READT  = .FALSE.
      IF ( TASK.EQ.1 ) THEN
         IF ( NLHS.GT.2 )
     $      COMPU = .TRUE.
         NLM = 4
      ELSE IF ( TASK.EQ.2 ) THEN
         IF ( NRHS.GT.2 ) THEN
            NLM = 7
            IF ( NLHS.GT.3 )
     $         COMPU = .TRUE.
            IF ( NLHS.GT.5 )
     $         COMPV = .TRUE.
            COMPQ = .TRUE.
            COMPT = .TRUE.
         ELSE
            READA  = .FALSE.
            READB  = .TRUE.
            READQG = .FALSE.
            IF ( NLHS.GT.1 )
     $         COMPU = .TRUE.
            IF ( NLHS.GT.3 )
     $         COMPV = .TRUE.
            NLM = 5
         END IF
      ELSE IF ( TASK.EQ.3 ) THEN
         NLM = 9
         IF ( NLHS.GT.5 )
     $      COMPU = .TRUE.
         IF ( NLHS.GT.7 )
     $      COMPV = .TRUE.
         IF ( NLHS.EQ.5 .OR. NLHS.EQ.7 .OR. NLHS.EQ.9 )
     $      COMPG = .TRUE.
         COMPT = .TRUE.
      ELSE IF ( TASK.EQ.4 ) THEN
         NLM = 4
      ELSE IF ( TASK.EQ.5 ) THEN
         COMPT = .TRUE.
         NLM = 2
      ELSE IF ( TASK.EQ.6 ) THEN
         IF ( NLHS.GT.2 )
     $      COMPU = .TRUE.
         IF ( NRHS.GE.7 )
     $      UPDU = .TRUE.
         NLM = 4
      ELSE IF ( TASK.EQ.7 .OR. TASK.EQ.8 ) THEN
         UPDU = .TRUE.
         UPDV = .TRUE.
         READQG = .FALSE.
         READT  = .TRUE.
         IF ( TASK.EQ.7 )
     $      COMPUS = .TRUE.
         IF ( TASK.EQ.8 .OR. NLHS.GT.1 )
     $      COMPUU = .TRUE.
         NLM = 2
      ELSE IF ( TASK.EQ.9 .OR. TASK.EQ.10 ) THEN
         UPDU = .TRUE.
         UPDV = .TRUE.
         READQG = .FALSE.
         READT  = .TRUE.
         IF ( NRHS.GT.7 )
     $      READG = .TRUE.
         IF ( TASK.EQ.9 )
     $      COMPUS = .TRUE.
         IF ( TASK.EQ.10 .OR. NLHS.GT.1 )
     $      COMPUU = .TRUE.
         NLM = 2
      END IF
      IF ( TASK.EQ.8 .OR. TASK.EQ.10 )
     $   NLM = 1
C
      IF ( NLHS.GT.NLM ) THEN
         WRITE( TEXT, '('' HAPACK_HAEIG requires at most '', I4,
     $            '' output arguments'')' ) NLM
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      IF ( ( TASK.LE.3 ) .AND. NRHS.GT.3 ) THEN
         CALL mexErrMsgTxt
     $           ( 'HAPACK_HAEIG requires at most 3 input arguments' )
      ELSE IF ( ( TASK.EQ.4 .OR. TASK.EQ.5 ) .AND. NRHS.GT.4 ) THEN
         CALL mexErrMsgTxt
     $           ( 'HAPACK_HAEIG requires at most 4 input arguments' )
      ELSE IF ( TASK.EQ.6 .AND. .NOT.COMPU .AND. NRHS.GT.5 ) THEN
         CALL mexErrMsgTxt
     $           ( 'HAPACK_HAEIG requires at most 5 input arguments' )
      ELSE IF ( TASK.EQ.6 .AND. NRHS.GT.7 ) THEN
         CALL mexErrMsgTxt
     $           ( 'HAPACK_HAEIG requires at most 7 input arguments' )
      ELSE IF ( TASK.GE.7 .AND. NRHS.GT.8 ) THEN
         CALL mexErrMsgTxt
     $           ( 'HAPACK_HAEIG requires at most 8 input arguments' )
      END IF
C
      IF ( ( TASK.EQ.1 .OR. TASK.EQ.3 ) .AND. NRHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $           ( 'HAPACK_HAEIG requires at least 3 input arguments' )
      ELSE IF ( ( TASK.EQ.4 .OR. TASK.EQ.5 ) .AND. NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $           ( 'HAPACK_HAEIG requires at least 4 input arguments' )
      ELSE IF ( TASK.EQ.6 .AND. NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $           ( 'HAPACK_HAEIG requires at least 5 input arguments' )
      ELSE IF ( ( TASK.EQ.9 .OR. TASK.EQ.10 ) .AND. NRHS.LT.7 ) THEN
         CALL mexErrMsgTxt
     $           ( 'HAPACK_HAEIG requires at least 7 input arguments' )
      ELSE IF ( ( TASK.EQ.7 .OR. TASK.EQ.8  ) .AND. NRHS.LT.8 ) THEN
         CALL mexErrMsgTxt
     $           ( 'HAPACK_HAEIG requires at least 8 input arguments' )
      END IF
C
C     Check BALANC.
C
      IF ( TASK.EQ.4 .OR. TASK.EQ.5 ) THEN
         IP = 4
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 .OR.
     $        mxIsChar( PRHS(IP) ).NE.1 )
     $      CALL mexErrMsgTxt( 'BALANC must be a character' )
         IERR = mxGetString( PRHS(IP), BALANC, 1 )
      END IF
C
      IF ( UPDU .AND. TASK.EQ.6 ) THEN
         IP = 4
      ELSE
         IP = 2
      END IF
C
      IF ( READA ) THEN
C
C        Check A.
C
         M = mxGetM( PRHS(IP) )
         N = mxGetN( PRHS(IP) )
         IF ( M.LT.0 .OR. N.LT.0 .OR. M.NE.N )
     $      CALL mexErrMsgTxt( 'A must be a square matrix' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'A must be a real matrix' )
         IP = IP + 1
      END IF
C
      IF ( READB ) THEN
C
C        Check B.
C
         M = mxGetM( PRHS(IP) )
         N = mxGetN( PRHS(IP) ) / 2
         IF ( M.LT.0 .OR. N.LT.0 .OR. M.NE.2*N .OR.
     $        N*2.NE.( mxGetN( PRHS(IP) ) ) )
     $      CALL mexErrMsgTxt( 'B must be a 2n-by-2n matrix' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'B must be a real matrix' )
         IP = IP + 1
      END IF
      N1 = N + 1
C
      IF ( READQG ) THEN
C
C        Check QG.
C
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( M.NE.N .OR. NN.NE.N1 )
     $      CALL mexErrMsgTxt( 'QG must be an n-by-n+1 matrix' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'QG must be a real matrix' )
         IP = IP + 1
      END IF
C
      IF ( READT ) THEN
C
C        Check T.
C
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( M.NE.N .OR. NN.NE.N )
     $      CALL mexErrMsgTxt( 'T must be an n-by-n matrix' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'T must be a real matrix' )
         IP = IP + 1
      END IF
C
      IF ( READG ) THEN
C
C        Check G.
C
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( M.NE.N .OR. NN.NE.N )
     $      CALL mexErrMsgTxt( 'G must be an n-by-n matrix' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'G must be a real matrix' )
         IP = IP + 1
      END IF
C
      IF ( UPDU ) THEN
C
C        Check U1 and U2.
C
         IF ( TASK.EQ.6 )
     $      IP = 2
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( M.NE.N .OR. NN.NE.N )
     $      CALL mexErrMsgTxt( 'U1 must be an n-by-n matrix' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $   CALL mexErrMsgTxt( 'U1 must be a real matrix' )
         IP = IP + 1
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( M.NE.N .OR. NN.NE.N )
     $      CALL mexErrMsgTxt( 'U2 must be an n-by-n matrix' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'U2 must be a real matrix' )
         IP = IP + 1
      END IF
      IF ( UPDV ) THEN
C
C        Check V1 and V2.
C
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( M.NE.N .OR. NN.NE.N )
     $      CALL mexErrMsgTxt( 'V1 must be an n-by-n matrix' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $   CALL mexErrMsgTxt( 'V1 must be a real matrix' )
         IP = IP + 1
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( M.NE.N .OR. NN.NE.N )
     $      CALL mexErrMsgTxt( 'V2 must be an n-by-n matrix' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'V2 must be a real matrix' )
         IP = IP + 1
      END IF
C
C     Check SELECT AND LOWER.
C
      IF ( TASK.EQ.6 ) THEN
         IF ( UPDU ) THEN
            IP = 6
         ELSE
            IP = 4
         END IF
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( .NOT.( ( M.EQ.1 .AND. NN.EQ.N ) .OR.
     $               ( M.EQ.N .AND. NN.EQ.1 ) ) )
     $      CALL mexErrMsgTxt( 'SELECT must be a vector of length n' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'SELECT must be a real vector' )
         IP = IP + 1
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( .NOT.( ( M.EQ.1 .AND. NN.EQ.N ) .OR.
     $               ( M.EQ.N .AND. NN.EQ.1 ) ) )
     $      CALL mexErrMsgTxt( 'LOWER must be a vector of length n' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'LOWER must be a real vector' )
      END IF
      IF ( TASK.EQ.7 .OR. TASK.EQ.8 ) THEN
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( .NOT.( ( M.EQ.1 .AND. NN.EQ.N ) .OR.
     $               ( M.EQ.N .AND. NN.EQ.1 ) ) )
     $      CALL mexErrMsgTxt( 'SELECT must be a vector of length n' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'SELECT must be a real vector' )
      END IF
C
C     Determine the lenghts of working arrays.
C
      LDA  = 1
      LDB  = 1
      LDG  = 1
      LDQG = 1
      LDU1 = 1
      LDU2 = 1
      LDV1 = 1
      LDV2 = 1
      LDUS = 1
      LDUU = 1
      IF ( READA )
     $   LDA = MAX( 1, N )
      IF ( READB ) THEN
         LDB = MAX( 1, 2*N )
         LDS = MAX( 1,   N )
      END IF
      IF ( READG )
     $   LDG = MAX( 1, N )
      IF ( READQG )
     $   LDQG = MAX( 1, N )
      IF ( COMPQ )
     $   LDQ = MAX( 1, N )
      IF ( COMPT .OR. READT )
     $   LDT = MAX( 1, N )
      IF ( COMPU .OR. UPDU ) THEN
         LDU1 = MAX( 1, N )
         LDU2 = LDU1
      END IF
      IF ( COMPV .OR. UPDU ) THEN
         LDV1 = MAX( 1, N )
         LDV2 = LDV1
      END IF
      IF ( COMPUS )
     $   LDUS = MAX( 1, 2*N )
      IF ( COMPUU )
     $   LDUU = MAX( 1, 2*N )
C
C     Workspace requirements to enable blocked algorithms.
C
      LDWORK = 1
      IF ( TASK.EQ.1 ) THEN
         CALL MB04PB( N, 1, DM, LDA, DM, LDQG, DM, DM, DM, -1, INFO )
         LDWORK = MAX( LDWORK, INT( DM(1) ) )
         IF ( COMPU ) THEN
            CALL MB04WP( N, 1, DM, LDU1, DM, LDU2, DM, DM, DM, -1,
     $                   INFO )
            LDWORK = MAX( LDWORK, INT( DM(1) ) )
         END IF
C
      ELSE IF ( TASK.EQ.2 ) THEN
C
         IF ( READB ) THEN
            CALL MB04TB( 'No Transpose', 'No Transpose', N, 1, DM, LDB,
     $                   DM, LDB, DM, LDB, DM, LDB, DM, DM, DM, DM, DM,
     $                   -1, INFO )
         ELSE
            CALL MB04TB( 'No Transpose', 'Transpose', N, 1, DM, LDT, DM,
     $                   LDA, DM, LDQG, DM, LDQ, DM, DM, DM, DM, DM, -1,
     $                   INFO )
         END IF
         LDWORK = MAX( LDWORK, INT( DM(1) ) )
         IF ( COMPU.OR.COMPV ) THEN
            CALL MB04WR( 'U', 'No Transpose', N, 1, DM, LDU1, DM, LDU2,
     $                   DM, DM, DM, -1, INFO )
            LDWORK = MAX( LDWORK, INT( DM(1) ) )
         END IF
C
      ELSE IF ( TASK.EQ.3 ) THEN
C
         IF ( COMPG ) THEN
            JOB = 'G'
         ELSE
            JOB = 'S'
         END IF
         IF ( COMPU ) THEN
            JOBU = 'U'
         ELSE
            JOBU = 'N'
         END IF
         IF ( COMPV ) THEN
            JOBV = 'V'
         ELSE
            JOBV = 'N'
         END IF
         CALL MB03XD( 'N', JOB, JOBU, JOBV, N, DM, LDA, DM, LDQG, DM,
     $                LDT, DM, LDU1, DM, LDU2, DM, LDV1, DM, LDV2, DM,
     $                DM, 1, DM, DM, -1, INFO )
         LDWORK = MAX( LDWORK, INT( DM(1) ) )
C
      ELSE IF ( TASK.EQ.5 ) THEN
C
         CALL MB03XD( BALANC, 'E', 'N', 'N', N, A, LDA, QG, LDQG, T,
     $                LDT, DM, 1, DM, 1, DM, 1, DM, 1, DM, DM, 1, DM,
     $                DM, -1, INFO )
         LDWORK = MAX( LDWORK, INT( DM(1) ) )
C
      ELSE IF ( TASK.EQ.6 ) THEN
C
         LDWORK = MAX( LDWORK, N )
C
      ELSE IF ( TASK.EQ.9 .OR. TASK.EQ.10 ) THEN
C
         IF ( COMPUS .AND. COMPUU ) THEN
            JOB = 'B'
         ELSE IF ( COMPUS ) THEN
            JOB = 'S'
         ELSE
            JOB = 'U'
         END IF
         IF ( READG ) THEN
            METH  = 'L'
            USCOL = 2*N
         ELSE
            METH  = 'S'
            USCOL = N
         END IF
         CALL MB03ZD( 'All', METH, JOB, 'N', 'B', LW, N, USCOL, 1,
     $                DM, DM, LDA, DM, LDT, DM, LDG, DM, LDU1, DM, LDU2,
     $                DM, LDV1, DM, LDV2, M, DM, DM, DM, LDUS, DM, LDUU,
     $                LW, IW, DM, -1, INFO )
         LDWORK = MAX( LDWORK, INT( DM(1) ) )
C
      END IF
C
C     Allocate variable dimension local arrays.
C     !Fortran 90/95
C
      IF ( TASK.NE.7 .AND. TASK.NE.8 )
     $   DWORK = mxCalloc( LDWORK, 8 )
      IF ( READA )
     $   A = mxCalloc( LDA*N, 8 )
      IF ( READB )
     $   B = mxCalloc( 4*N*N, 8 )
      IF ( READG )
     $   G = mxCalloc( LDG*N, 8 )
      IF ( READQG )
     $   QG = mxCalloc( LDQG*N1, 8 )
      IF ( COMPQ )
     $   Q = mxCalloc( LDQ*N, 8 )
      IF ( COMPT .OR. READT )
     $   T = mxCalloc( LDT*N, 8 )
      IF ( TASK.EQ.1 .OR. TASK.EQ.2 ) THEN
         CSL  = mxCalloc( 2*N, 8 )
         TAUL = mxCalloc(   N, 8 )
         IF ( TASK.EQ.2 ) THEN
            CSR  = mxCalloc( 2*N, 8 )
            TAUR = mxCalloc(   N, 8 )
         END IF
      END IF
      IF ( TASK.EQ.3 .OR. TASK.GE.5 ) THEN
         WR = mxCalloc( N, 8 )
         WI = mxCalloc( N, 8 )
      END IF
      IF ( COMPU .OR. UPDU ) THEN
         U1 = mxCalloc( LDU1*N, 8 )
         U2 = mxCalloc( LDU2*N, 8 )
      END IF
      IF ( COMPV .OR. UPDV ) THEN
         V1 = mxCalloc( LDV1*N, 8 )
         V2 = mxCalloc( LDV2*N, 8 )
      END IF
      IF ( TASK.EQ.4 .OR. TASK.EQ.5 )
     $   SCALE = mxCalloc( N, 8 )
      IF ( TASK.EQ.6 .OR. TASK.EQ.7 .OR. TASK.EQ.8 ) THEN
         ALLOCATE ( SELECT(N) )
         IF ( TASK.NE.6 )
     $      ALLOCATE ( IWORK(M) )
      END IF
      IF ( TASK.EQ.6 )
     $   ALLOCATE ( LOWER(N) )
      IF ( TASK.EQ.9 .OR. TASK.EQ.10 ) THEN
         IF ( READG ) THEN
            ALLOCATE ( IWORK(2*N), LWORK(2*N) )
         ELSE
            ALLOCATE ( IWORK(N), LWORK(0) )
         END IF
      END IF
      JG = MIN( N1, 2 )
C
C     Copy inputs from MATLAB workspace to locally allocated arrays.
C
      IP = 1
      IF ( TASK.EQ.6 ) THEN
         IF ( UPDU ) THEN
            IP = IP + 1
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(U1),
     $                             N*N )
            IP = IP + 1
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(U2),
     $                             N*N )
         ELSE IF ( COMPU ) THEN
            CALL DLASET( 'All', N, N, ZERO, ONE,  %VAL(U1), LDU1 )
            CALL DLASET( 'All', N, N, ZERO, ZERO, %VAL(U2), LDU2 )
         END IF
      END IF
      IF ( READA ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(A), N*N )
      END IF
      IF ( READQG ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(QG),
     $                          N*N1 )
      END IF
      IF ( READB ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(B),
     $                          4*N*N )
      END IF
      IF ( READT ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(T), N*N )
      END IF
      IF ( READG ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(G), N*N )
      END IF
      IF ( UPDU .AND. TASK.GE.7 ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(U1), N*N )
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(U2), N*N )
      END IF
      IF ( UPDV .AND. TASK.GE.7 ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(V1), N*N )
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL(V2), N*N )
      END IF
      IF ( TASK.EQ.6 .OR. TASK.EQ.7 .OR. TASK.EQ.8 ) THEN
         ALLOCATE ( DSEL(N) )
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), DSEL, N )
         DO 10  I = 1, N
            SELECT(I) = DSEL(I).EQ.ONE
   10    CONTINUE
         IF ( TASK.EQ.6 ) THEN
            IP = IP + 1
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), DSEL, N )
            DO 20  I = 1, N
               LOWER(I) = DSEL(I).EQ.ONE
   20       CONTINUE
         END IF
         DEALLOCATE ( DSEL )
      END IF
C
C     Determine dimensions of DWORK, US and UU if TASK = 7,8.
C
      IF ( TASK.EQ.7 .OR. TASK.EQ.8 ) THEN
         ALLOCATE ( TW(LDA,N) )
         CALL DLACPY( 'Lower', N, N, %VAL(A), LDA, TW, LDA )
         M = 0
         PAIR = .FALSE.
         DO 30 I = 1, N
            IF ( PAIR ) THEN
               PAIR = .FALSE.
            ELSE
               IF ( I.LT.N ) THEN
                  IF ( TW(I+1,I).EQ.ZERO ) THEN
                     IF ( SELECT(I) )
     $                  M = M + 1
                  ELSE
                     PAIR = .TRUE.
                     IF ( SELECT(I) .OR. SELECT(I+1) )
     $                  M = M + 2
                  END IF
               ELSE
                  IF ( SELECT(N) )
     $               M = M + 1
               END IF
            END IF
   30    CONTINUE
         DEALLOCATE ( TW )
         USCOL  = M
         IF ( COMPUS .AND. COMPUU ) THEN
            JOB = 'B'
         ELSE IF ( COMPUS ) THEN
            JOB = 'S'
         ELSE
            JOB = 'U'
         END IF
         CALL MB03ZD( 'Selected', 'S', JOB, 'N', 'B', SELECT, N, USCOL,
     $                1, DM, %VAL(A), LDA, DM, LDT, DM, 1, DM, LDU1, DM,
     $                LDU2, DM, LDV1, DM, LDV2, M, DM, DM, DM, LDUS,
     $                DM, LDUU, LW, IW, DM, -1, INFO )
         LDWORK = MAX( LDWORK, INT( DM(1) ) )
         DWORK = mxCalloc( LDWORK, 8 )
      END IF
      IF ( COMPUS )
     $   US = mxCalloc( LDUS*USCOL, 8 )
      IF ( COMPUU )
     $   UU = mxCalloc( LDUU*USCOL, 8 )
C
C     Do the actual computations.
C
      IF ( TASK.EQ.1 ) THEN
         CALL MB04PB( N, 1, %VAL(A), LDA, %VAL(QG), LDQG, %VAL(CSL),
     $                %VAL(TAUL), %VAL(DWORK), LDWORK, INFO )
         IF ( INFO.EQ.0 ) THEN
            IF ( COMPU ) THEN
               CALL DLACPY( 'Lower', N, N, %VAL(A), LDA, %VAL(U1), LDU1
     $                    )
               CALL DLACPY( 'Lower', N, N, %VAL(QG), LDQG, %VAL(U2),
     $                      LDU2 )
               CALL MB04WP( N, 1, %VAL(U1), LDU1, %VAL(U2), LDU2,
     $                      %VAL(CSL), %VAL(TAUL), %VAL(DWORK), LDWORK,
     $                      INFO )
            END IF
            IF ( N.GT.2 ) THEN
               CALL MA02LD( 'Lower', N-2, N-2, 3, 1, ZERO, ZERO,
     $                      %VAL(A), LDA )
            END IF
            IF ( N.GT.1 ) THEN
               CALL MA02LD( 'Lower', N-1, N-1, 2, 1, ZERO, ZERO,
     $                      %VAL(QG), LDQG )
            END IF
         END IF
C
      ELSE IF ( TASK.EQ.2 ) THEN
C
         IF ( READB ) THEN
            B11 = mxCalloc( LDS*N, 8 )
            B12 = mxCalloc( LDS*N, 8 )
            B21 = mxCalloc( LDS*N, 8 )
            B22 = mxCalloc( LDS*N, 8 )
            CALL MA02KD( 'All', N, N,  1,  1, 1, 1, %VAL(B), LDB,
     $                   %VAL(B11), LDS )
            CALL MA02KD( 'All', N, N, N1,  1, 1, 1, %VAL(B), LDB,
     $                   %VAL(B21), LDS )
            CALL MA02KD( 'All', N, N, 1,  N1, 1, 1, %VAL(B), LDB,
     $                   %VAL(B12), LDS )
            CALL MA02KD( 'All', N, N, N1, N1, 1, 1, %VAL(B), LDB,
     $                   %VAL(B22), LDS )
            CALL MB04TB( 'No Transpose', 'No Transpose', N, 1,
     $                   %VAL(B11), LDS, %VAL(B22), LDS, %VAL(B12), LDS,
     $                   %VAL(B21), LDS, %VAL(CSL), %VAL(CSR),
     $                   %VAL(TAUL), %VAL(TAUR), %VAL(DWORK), LDWORK,
     $                   INFO )
            IF ( INFO.EQ.0 ) THEN
               IF ( COMPU ) THEN
                  CALL DLACPY( 'Lower', N, N, %VAL(B11), LDS, %VAL(U1),
     $                         LDU1 )
                  CALL DLACPY( 'Lower', N, N, %VAL(B21), LDS, %VAL(U2),
     $                         LDU2 )
               END IF
               IF ( COMPV ) THEN
                  CALL DLACPY( 'Upper', N, N, %VAL(B21), LDS, %VAL(V2),
     $                         LDV2 )
                  CALL DLACPY( 'Upper', N, N, %VAL(B22), LDS, %VAL(V1),
     $                         LDV1 )
               END IF
               CALL MA02KD( 'Upper', N, N, 1, 1,  1,  1, %VAL(B11), LDS,
     $                      %VAL(B), LDB )
               CALL MA02LD( 'Lower', 2*N-1, N, 2, 1, ZERO, ZERO,
     $                      %VAL(B), LDB )
               CALL MA02KD( 'All',   N, N, 1, 1,  1, N1, %VAL(B12), LDS,
     $                      %VAL(B), LDB )
               CALL MA02KV( N, 1, N*(LDB+1)+1, %VAL(B22), %VAL(B) )
               CALL MA02KD( 'Lower', N, N-1, 1, 2, N1, N1+1, %VAL(B22),
     $                      LDS, %VAL(B), LDB )
               IF ( N.GT.2 )
     $            CALL MA02LD( 'Upper', N-2, N-2, N1, N+3, ZERO, ZERO,
     $                         %VAL(B), LDB )
               TRAN = 'N'
            END IF
C
         ELSE
C
            CALL DLACPY( 'All', N, N, %VAL(A), LDA, %VAL(T), LDT )
            CALL DLASCL( 'General', 0, 0, ONE, -ONE, N, N, %VAL(A), LDA,
     $                   IERR )
            CALL DLACPY( 'Lower', N, N, %VAL(QG), LDQG, %VAL(Q), LDQ )
            CALL MA02ED( 'Lower', N, %VAL(Q), LDQ )
            DW = mxCalloc( LDQG*N, 8 )
            CALL MA02KD( 'Upper', N, N, 1, JG, 1, 1, %VAL(QG), LDQG,
     $                   %VAL(DW), LDQG )
            CALL MA02ED( 'Upper',  N, %VAL(DW), LDQG )
            CALL DLACPY( 'All', N, N, %VAL(DW), LDQG, %VAL(QG), LDQG )
            CALL mxFree( DW )
C
            CALL MB04TB( 'No Transpose', 'Transpose', N, 1, %VAL(T),
     $                   LDT, %VAL(A), LDA, %VAL(QG), LDQG, %VAL(Q),
     $                   LDQ, %VAL(CSL), %VAL(CSR), %VAL(TAUL),
     $                   %VAL(TAUR), %VAL(DWORK), LDWORK, INFO )
            IF ( INFO.EQ.0 ) THEN
               IF ( COMPU ) THEN
                  CALL DLACPY( 'Lower', N, N, %VAL(T), LDT, %VAL(U1),
     $                         LDU1 )
                  CALL DLACPY( 'Lower', N, N, %VAL(Q), LDQ, %VAL(U2),
     $                         LDU2 )
               END IF
               IF ( COMPV ) THEN
                  CALL DLACPY( 'Lower', N, N, %VAL(A), LDA, %VAL(V1),
     $                         LDV1 )
                  CALL DLACPY( 'Upper', N, N, %VAL(Q), LDQ, %VAL(V2),
     $                         LDV2 )
               END IF
               IF ( N.GT.1 ) THEN
                  CALL MA02LD( 'Lower', N-1, N-1, 2, 1, ZERO, ZERO,
     $                         %VAL(T), LDT )
                  IF ( N.GT.2 ) THEN
                     CALL MA02LD( 'Lower', N-2, N-2, 3, 1, ZERO, ZERO,
     $                            %VAL(A), LDA )
                  END IF
               END IF
            END IF
            TRAN = 'T'
C
         END IF
         IF ( COMPU .AND. INFO.EQ.0 )
     $      CALL MB04WR( 'U', 'No Transpose', N, 1, %VAL(U1), LDU1,
     $                   %VAL(U2), LDU2, %VAL(CSL), %VAL(TAUL),
     $                   %VAL(DWORK), LDWORK, INFO )
         IF ( COMPV .AND. INFO.EQ.0 )
     $      CALL MB04WR( 'V', TRAN, N, 1, %VAL(V1), LDV1, %VAL(V2),
     $                   LDV2, %VAL(CSR), %VAL(TAUR), %VAL(DWORK),
     $                   LDWORK, INFO )
C
      ELSE IF ( TASK.EQ.3 ) THEN
C
         CALL MB03XD( 'N', JOB, JOBU, JOBV, N, %VAL(A), LDA, %VAL(QG),
     $                LDQG, %VAL(T), LDT, %VAL(U1), LDU1, %VAL(U2),
     $                LDU2, %VAL(V1), LDV1, %VAL(V2), LDV2, %VAL(WR),
     $                %VAL(WI), ILO, %VAL(DWORK), %VAL(DWORK), LDWORK,
     $                INFO )
         DO 40 J = 2, N1
            CALL MA02KV( N, (J-1)*LDQG+1, (J-2)*LDQG+1, %VAL(QG),
     $                   %VAL(QG) )
   40    CONTINUE
         ALLOCATE ( TW(LDQG,2) )
         CALL DCOPY( N, %VAL(WR), 1, TW, 1 )
         CALL DCOPY( N, %VAL(WI), 1, TW(1,2), 1 )
         I = 1
   50    CONTINUE
C        WHILE( I <= N ) DO
         IF ( I.LE.N ) THEN
            TW(I,1) = -TW(I,1)
            IF ( TW(I,1).NE.ZERO .AND. TW(I,2).NE.ZERO ) THEN
               I = I + 1
               TW(I,2) = -TW(I,2)
            END IF
            I = I + 1
            GO TO 50
         END IF
C        END WHILE 50
         CALL DCOPY( N, TW, 1, %VAL(WR), 1 )
         CALL DCOPY( N, TW(1,2), 1, %VAL(WI), 1 )
         DEALLOCATE ( TW )
C
      ELSE IF ( TASK.EQ.4 ) THEN
C
         CALL MB04DD( BALANC, N, %VAL(A), LDA, %VAL(QG), LDQG, ILO,
     $                %VAL(SCALE), INFO )
C
      ELSE IF ( TASK.EQ.5 ) THEN
C
         CALL MB03XD( BALANC, 'E', 'N', 'N', N, %VAL(A), LDA, %VAL(QG),
     $                LDQG, %VAL(T), LDT, %VAL(DWORK), 1, %VAL(DWORK),
     $                1, %VAL(DWORK), 1, %VAL(DWORK), 1, %VAL(WR),
     $                %VAL(WI), ILO, %VAL(SCALE), %VAL(DWORK), LDWORK,
     $                INFO )
         ALLOCATE ( TW(LDQG,2) )
         CALL DCOPY( N, %VAL(WR), 1, TW, 1 )
         CALL DCOPY( N, %VAL(WI), 1, TW(1,2), 1 )
         I = 1
   60    CONTINUE
C        WHILE( I <= N ) DO
         IF ( I.LE.N ) THEN
            TW(I,1) = -TW(I,1)
            IF ( TW(I,1).NE.ZERO .AND. TW(I,2).NE.ZERO ) THEN
               I = I + 1
               TW(I,2) = -TW(I,2)
            END IF
            I = I + 1
            GO TO 60
         END IF
C        END WHILE 60
         CALL DCOPY( N, TW, 1, %VAL(WR), 1 )
         CALL DCOPY( N, TW(1,2), 1, %VAL(WI), 1 )
         DEALLOCATE ( TW )
C
      ELSE IF ( TASK.EQ.6 ) THEN
C
         DW = mxCalloc( N, 8 )
         CALL MA02KV( N, 1, 1, %VAL(QG), %VAL(DW) )
         DO 70 J = 2, N1
            CALL MA02KV( N, (J-1)*LDQG+1, (J-2)*LDQG+1, %VAL(QG),
     $                   %VAL(QG) )
   70    CONTINUE
         IF ( COMPU ) THEN
            CALL MB03TD( 'Hamiltonian', 'Update', SELECT, LOWER, N,
     $                   %VAL(A), LDA, %VAL(QG), LDQG, %VAL(U1), LDU1,
     $                   %VAL(U2), LDU2, %VAL(WR), %VAL(WI), M,
     $                   %VAL(DWORK), LDWORK, INFO )
         ELSE
            CALL MB03TD( 'Hamiltonian', 'No Update', SELECT, LOWER,
     $                   N, %VAL(A), LDA, %VAL(QG), LDQG, %VAL(DWORK),
     $                   1, %VAL(DWORK), 1, %VAL(WR), %VAL(WI), M,
     $                   %VAL(DWORK), LDWORK, INFO )
         END IF
         DO 80 J = N1, 2, -1
            CALL MA02KV( N, (J-2)*LDQG+1, (J-1)*LDQG+1, %VAL(QG),
     $                   %VAL(QG) )
   80    CONTINUE
         CALL MA02KV( N, 1, 1, %VAL(DW), %VAL(QG) )
         CALL mxFree( DW )
C
      ELSE IF ( TASK.EQ.7 .OR. TASK.EQ.8 ) THEN
C
         CALL MB03ZD( 'Selected', 'S', JOB, 'N', 'B', SELECT, N, USCOL,
     $                1, %VAL(DWORK), %VAL(A), LDA, %VAL(T), LDT,
     $                %VAL(DWORK), 1, %VAL(U1), LDU1, %VAL(U2), LDU2,
     $                %VAL(V1), LDV1, %VAL(V2), LDV2, M, %VAL(WR),
     $                %VAL(WI), %VAL(US), LDUS, %VAL(UU), LDUU, LWORK,
     $                IWORK, %VAL(DWORK), LDWORK, INFO )
C
      ELSE IF ( TASK.EQ.9 .OR. TASK.EQ.10 ) THEN
C
         CALL MB03ZD( 'All', METH, JOB, 'N', 'B', SELECT, N, USCOL, 1,
     $                %VAL(DWORK), %VAL(A), LDA, %VAL(T), LDT, %VAL(G),
     $                LDG, %VAL(U1), LDU1, %VAL(U2), LDU2, %VAL(V1),
     $                LDV1, %VAL(V2), LDV2, M, %VAL(WR), %VAL(WI),
     $                %VAL(US), LDUS, %VAL(UU), LDUU, LWORK, IWORK,
     $                %VAL(DWORK), LDWORK, INFO )
      END IF
C
C     Copy output to MATLAB workspace.
C
      IP = 1
      IF ( INFO.EQ.0 ) THEN
         IF ( COMPUS ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( 2*N, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL(US), mxGetPr( PLHS(IP) ),
     $                             2*N*M )
            IP = IP + 1
         END IF
         IF ( COMPUU ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( 2*N, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL(UU), mxGetPr( PLHS(IP) ),
     $                             2*N*M )
            IP = IP + 1
         END IF
         IF ( COMPU ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(U1), mxGetPr( PLHS(IP) ),
     $                             N*N )
            IP = IP + 1
         END IF
         IF ( COMPU .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(U2), mxGetPr( PLHS(IP) ),
     $                             N*N )
            IP = IP + 1
         END IF
         IF ( COMPV .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(V1), mxGetPr( PLHS(IP) ),
     $                             N*N )
            IP = IP + 1
         END IF
         IF ( COMPV .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(V2), mxGetPr( PLHS(IP) ),
     $                             N*N )
            IP = IP + 1
         END IF
         IF ( (   TASK.EQ.1 .OR.
     $          ( TASK.EQ.2 .AND. READA ) .OR.
     $            TASK.EQ.3 .OR.
     $          ( TASK.EQ.4 .AND.   NLHS.GT.2 ) .OR.
     $            TASK.EQ.6 ) .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(A), mxGetPr( PLHS(IP) ),
     $                             N*N )
            IP = IP + 1
         END IF
         IF ( ( ( TASK.EQ.2 .AND. READA ) .OR.
     $            TASK.EQ.3 ) .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(T), mxGetPr( PLHS(IP) ),
     $                             N*N )
            IP = IP + 1
         END IF
         IF ( ( TASK.EQ.1 .OR.
     $        ( TASK.EQ.4 .AND. NLHS.GT.2 ) .OR.
     $          TASK.EQ.6 ) .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N1, 0 )
            CALL mxCopyReal8ToPtr( %VAL(QG), mxGetPr( PLHS(IP) ),
     $                             N*N1 )
            IP = IP + 1
         END IF
         IF ( ( ( TASK.EQ.2 .AND. READA ) .OR.
     $            TASK.EQ.3 .AND. COMPG ) .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(QG), mxGetPr( PLHS(IP) ),
     $                             N*N )
            IP = IP + 1
         END IF
         IF ( READB .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( 2*N, 2*N, 0 )
            CALL mxCopyReal8ToPtr( %VAL(B), mxGetPr( PLHS(IP) ),
     $                             4*N*N )
            IP = IP + 1
         END IF
         IF ( ( TASK.EQ.3 .OR. TASK.EQ.5 ) .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, 1, 0 )
            CALL mxCopyReal8ToPtr( %VAL(WR), mxGetPr( PLHS(IP) ),
     $                             N )
            IP = IP + 1
         END IF
         IF ( ( TASK.EQ.3 .OR. TASK.EQ.5 ) .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, 1, 0 )
            CALL mxCopyReal8ToPtr( %VAL(WI), mxGetPr( PLHS(IP) ),
     $                             N )
            IP = IP + 1
         END IF
         IF ( TASK.EQ.4 .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( N, 1, 0 )
            CALL mxCopyReal8ToPtr( %VAL(SCALE), mxGetPr( PLHS(IP) ),
     $                             N )
            IP = IP + 1
         END IF
         IF ( TASK.EQ.4 .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
            TEMP = DBLE( ILO )
            CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS(IP) ), 1 )
         END IF
      END IF
C
C     Deallocate local arrays.
C     !Fortran 90/95
C
      CALL mxFree( DWORK )
      IF ( READA )
     $   CALL mxFree( A )
      IF ( READB )
     $   CALL mxFree( B )
      IF ( READG )
     $   CALL mxFree( G )
      IF ( READQG )
     $   CALL mxFree( QG )
      IF ( COMPQ )
     $   CALL mxFree( Q )
      IF ( COMPT .OR. READT )
     $   CALL mxFree( T )
      IF ( TASK.EQ.1 .OR. TASK.EQ.2 ) THEN
         CALL mxFree( CSL  )
         CALL mxFree( TAUL )
         IF ( TASK.EQ.2 ) THEN
            CALL mxFree( CSR  )
            CALL mxFree( TAUR )
         END IF
      END IF
      IF ( TASK.EQ.3 .OR. TASK.GE.5 ) THEN
         CALL mxFree( WI )
         CALL mxFree( WR )
      END IF
      IF ( COMPU .OR. UPDU ) THEN
         CALL mxFree( U1 )
         CALL mxFree( U2 )
      END IF
      IF ( COMPV .OR. UPDV ) THEN
         CALL mxFree( V1 )
         CALL mxFree( V2 )
      END IF
      IF ( COMPUS )
     $   CALL mxFree( US )
      IF ( COMPUU )
     $   CALL mxFree( UU )
      IF ( TASK.EQ.4 .OR. TASK.EQ.5 )
     $   CALL mxFree( SCALE )
      IF ( TASK.EQ.6 .OR. TASK.EQ.7 .OR. TASK.EQ.8 ) THEN
         DEALLOCATE ( SELECT )
         IF ( TASK.NE.6 )
     $      DEALLOCATE ( IWORK )
      END IF
      IF ( TASK.EQ.6 )
     $   DEALLOCATE ( LOWER )
      IF ( TASK.EQ.9 .OR. TASK.EQ.10  )
     $   DEALLOCATE ( IWORK, LWORK )
C
C     Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( TASK.EQ.1 ) THEN
            IF ( COMPU ) THEN
               WRITE( TEXT, '('' INFO = '', I4,
     $                '' ON EXIT FROM MB04PB/MB04WP'')' )  INFO
            ELSE
               WRITE( TEXT,
     $                '('' INFO = '', I4, '' ON EXIT FROM MB04PB'')' )
     $                INFO
            END IF
         ELSE IF ( TASK.EQ.2 ) THEN
            IF ( COMPU.OR.COMPV ) THEN
               WRITE( TEXT, '('' INFO = '', I4,
     $                '' ON EXIT FROM MB04TB/MB04WR'')' )  INFO
            ELSE
               WRITE( TEXT,
     $                '('' INFO = '', I4, '' ON EXIT FROM MB04TB'')' )
     $                INFO
            END IF
         ELSE IF ( TASK.EQ.3 .OR. TASK.EQ.5 ) THEN
            IF ( INFO.GT.0 ) THEN
               WRITE( TEXT,
     $         '('' The periodic QR algorithm failed to converge. '')' )
            ELSE
               WRITE( TEXT,
     $                '('' INFO = '', I4, '' ON EXIT FROM MB03XD'')' )
     $                INFO
            END IF
         ELSE IF ( TASK.EQ.4 ) THEN
            WRITE( TEXT, '('' INFO = '', I4,
     $                '' ON EXIT FROM MB04DD'')' )  INFO
         ELSE IF ( TASK.EQ.6 ) THEN
            IF ( INFO.GT.0 ) THEN
               WRITE( TEXT,
     $         '('' Reordering of eigenvalues failed. '')' )
            ELSE
               WRITE( TEXT,
     $                '('' INFO = '', I4, '' ON EXIT FROM MB03TD'')' )
     $                INFO
            END IF
         ELSE IF ( TASK.GE.7 ) THEN
            IF ( INFO.EQ.1 ) THEN
               WRITE( TEXT,
     $            '('' Some of the eigenvalues are too close to the '',
     $            ''imaginary axis. '')' )
            ELSE IF ( INFO.EQ.2 .OR. INFO.EQ.4 ) THEN
               WRITE( TEXT,
     $         '('' Reordering of eigenvalues failed. '')' )
            ELSE IF ( INFO.EQ.3 ) THEN
               WRITE( TEXT,
     $         '('' The QR algorithm failed to converge. '')' )
            ELSE
               WRITE( TEXT,
     $                '('' INFO = '', I4, '' ON EXIT FROM MB03ZD'')' )
     $                INFO
            END IF
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of HAPACK_HAEIG ***
      END
