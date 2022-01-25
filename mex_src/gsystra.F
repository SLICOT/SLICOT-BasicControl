#include "fintrf.h"
C GSYSTRA.F  - Gateway function to perform various equivalence
C              transformations for descriptor systems with scaling,
C              generalized Schur form, etc., using SLICOT routines
C              TG01AD, TG01BD, TG01CD, TG01DD, TG01ED, TG01FD, and
C              TG01WD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [(Ao,Eo(,Bo)(,Co)(lscal,rscal)(,Q)(,Z)(,ranks)(,ev,db)] =
C                           GSYSTRA(task,A,E(,B)(,C)(,flag)(,Q1,Z1))
C
C   [Ao,Eo,Bo,Co(,lscal,rscal)] = GSYSTRA(0,A,E,B,C(,flag))
C   [Ao,Eo,Bo,Co(,Q,Z)]         = GSYSTRA(1,A,E,B,C(,flag)(,Q1,Z1))
C   [Ao,Eo,Bo(,Q)]              = GSYSTRA(2,A,E,B(,flag)(,Q1))
C   [Ao,Eo,Co(,Z)]              = GSYSTRA(3,A,E,C(,flag)(,Z1))
C   [Ao,Eo,Bo,Co(,Q,Z)(,ranks)] = GSYSTRA(4,A,E,B,C(,flag))
C   [Ao,Eo,Bo,Co(,Q,Z)(,ranks)] = GSYSTRA(5,A,E,B,C(,flag)(,Q1,Z1))
C   [Ao,Eo,Bo,Co(,Q,Z)(,ev,db)] = GSYSTRA(6,A,E,B,C)
C
C Purpose:
C   To perform one of the equivalence transformations specified by the
C   value of the parameter task (task = 0, 1, ..., 6), for a descriptor
C   triple (A-lambda E,B,C), or its parts:
C
C   0) To balance the matrices of the system pencil
C
C      S =  ( A  B ) - lambda ( E  0 ) .                             (1)
C           ( C  0 )          ( 0  0 )
C
C   Balancing involves diagonal equivalence transformations
C   (Dl*A*Dr - lambda Dl*E*Dr, Dl*B, C*Dr) applied to the system
C   (A-lambda E,B,C) to make the rows and columns of system pencil
C   matrices  diag(Dl,I)*S*diag(Dr,I)  as close in norm as possible.
C   Balancing may reduce the 1-norms of the matrices of the system
C   pencil S. The balancing can be performed optionally on the
C   following particular system pencils
C
C      S = A-lambda E,
C      S = ( A-lambda E  B ),    or
C      S = ( A-lambda E ).
C          (     C      )
C
C   1) To reduce the matrices A and E of the system pencil (1) to
C   generalized upper Hessenberg form using orthogonal transformations,
C
C      Q'*A*Z = H,   Q'*E*Z = T,
C
C   where H is upper Hessenberg, T is upper triangular, Q and Z
C   are orthogonal, and ' means transpose. The corresponding
C   transformations, written compactly as diag(Q',I)*S*diag(Z,I),
C   are also applied to B and C, getting Q'*B and C*Z.
C   The orthogonal matrices Q and Z are determined as products of
C   Givens rotations. They may either be formed explicitly, or they
C   may be postmultiplied into input matrices Q1 and Z1, so that
C
C      Q1*A*Z1' = (Q1*Q)*H*(Z1*Z)'
C      Q1*E*Z1' = (Q1*Q)*T*(Z1*Z)'.
C
C   2) To reduce the descriptor system pair (A-lambda E,B) to the
C   QR-coordinate form by computing an orthogonal transformation
C   matrix Q such that the transformed descriptor system pair
C   (Q'*A-lambda Q'*E, Q'*B) has the descriptor matrix Q'*E in an
C   upper trapezoidal form. The left orthogonal transformations
C   performed to reduce E can be optionally accumulated.
C
C   3) To reduce the descriptor system pair (C,A-lambda E) to the
C   RQ-coordinate form by computing an orthogonal transformation
C   matrix Z such that the transformed descriptor system pair
C   (C*Z,A*Z-lambda E*Z) has the descriptor matrix E*Z in an upper
C   trapezoidal form. The right orthogonal transformations performed
C   to reduce E can be optionally accumulated.
C
C   4) To compute for the descriptor system (A-lambda E,B,C) the
C   orthogonal transformation matrices Q and Z such that the
C   transformed system (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z) is in an
C   SVD (singular value decomposition) coordinate form with
C   the system matrices Q'*A*Z and Q'*E*Z in the form
C
C               ( A11  A12 )             ( Er  0 )
C      Q'*A*Z = (          ) ,  Q'*E*Z = (       ) ,                 (2)
C               ( A21  A22 )             (  0  0 )
C
C   where Er is an invertible diagonal matrix having on the diagonal
C   the decreasingly ordered nonzero singular values of E. Optionally,
C   the A22 matrix can be further reduced to the SVD form
C
C               ( Ar  0 )
C               (       ) ,                                          (3)
C               (  0  0 )
C
C   where Ar is an invertible diagonal matrix having on the diagonal
C   the decreasingly ordered nonzero singular values of A22. The left
C   and/or right orthogonal transformations performed to reduce E
C   and A22 are accumulated.
C
C   5) To compute for the descriptor system (A-lambda E,B,C) the
C   orthogonal transformation matrices Q and Z such that the
C   transformed system (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z) is in an
C   SVD-like coordinate form (2), where Er is an upper triangular
C   invertible matrix. Optionally, the A22 matrix can be further
C   reduced to the form
C
C               ( Ar  X )
C               (       ) ,
C               (  0  0 )
C
C   with Ar an upper triangular invertible matrix, and X either a full
C   or a zero matrix. The left and/or right orthogonal transformations
C   performed to reduce E and A22 can be optionally accumulated.
C
C   6) To reduce the pair (A,E) to a real generalized Schur form
C   by using an orthogonal equivalence transformation
C   (A,E) <-- (Q'*A*Z,Q'*E*Z) and to apply the transformation
C   to the matrices B and C: B <-- Q'*B and C <-- C*Z.
C
C Input parameters:
C   task   - integer specifying the computations to be performed.
C            task = 0 :  balance the descriptor system;
C            task = 1 :  compute the generalized upper Hessenberg form;
C            task = 2 :  compute QR-coordinate form;
C            task = 3 :  compute RQ-coordinate form;
C            task = 4 :  compute SVD coordinate form;
C            task = 5 :  compute SVD-like coordinate form;
C            task = 6 :  compute real generalized Schur form.
C   A      - the nl-by-n state dynamics matrix A, where nl = n, if
C            task = 1 or task = 6, and nl = l, otherwise.
C   E      - the nl-by-n descriptor matrix E. If task = 1 and jobe = 1,
C            the matrix E is assumed upper triangular.
C   B      - the nl-by-m input/state matrix B.
C   C      - the p-by-n state/output matrix C.
C   flag   - (optional) real vector specifying various options,
C            depending on task.
C
C            For task = 0, flag has length 2, and contains:
C            flag(1) = job : indicates which matrices are involved in
C            balancing, as follows:
C               job = 0 :  All matrices are involved in balancing;
C               job = 1 :  B, A and E are involved in balancing;
C               job = 2 :  C, A and E are involved in balancing;
C               job = 3 :  B and C are not involved in balancing.
C            flag(2) = thresh >= 0 : threshold value for magnitude of
C            elements: elements with magnitude less than or equal to
C            thresh are ignored for balancing.
C            Default :  flag = [ 0; 0 ].
C
C            For task = 1, flag has length 3 (at most), and contains:
C            flag(1) = jobe : specifies whether E is a general square
C            or an upper triangular matrix, as follows:
C               jobe = 0 :  E is a general square matrix;
C               jobe = 1 :  E is an upper triangular matrix.
C            flag(2) = compq : indicates what should be done with
C            matrix Q, as follows:
C               compq = 0 :  do not compute Q;
C               compq = 1 :  Q is initialized to the unit matrix, and
C                            the orthogonal matrix Q is returned;
C               compq = 2 :  Q is initialized to an orthogonal matrix Q1
C                            and the product Q1*Q is returned.
C            flag(3) = compz : indicates what should be done with
C            matrix Z, as follows:
C               compz = 0 :  do not compute Z;
C               compz = 1 :  Z is initialized to the unit matrix, and
C                            the orthogonal matrix Z is returned;
C               compz = 2 :  Z is initialized to an orthogonal matrix Z1
C                            and the product Z1*Z is returned.
C            Default :  flag = [ 0; 0; 0 ].
C
C            For task = 2, flag has length 1, and contains compq above.
C            Default :  flag = 0.
C
C            For task = 3, flag has length 1, and contains compz above.
C            Default :  flag = 0.
C
C            For task = 4, flag has length 2, and contains:
C            flag(1) = joba : specifies whether or not A22 should be
C            reduced, as follows:
C               joba = 0 :  do not reduce A22;
C               joba = 1 :  reduce A22 to an SVD form.
C            flag(2) = tol  : the tolerance to be used in determining
C            the rank of E and of A22. If tol > 0, then singular values
C            less than tol*svmax are treated as zero, where svmax is the
C            maximum singular value of E or of its estimate for A.
C            If tol <= 0, the default tolerance toldef = eps*l*n is used
C            instead, where eps is the machine precision. tol < 1.
C            Default :  flag = [ 0; 0 ].
C
C            For task = 5, flag has length 4, and contains:
C            flag(1) = joba : specifies whether or not A22 should be
C            reduced, as follows:
C               joba = 0 :  do not reduce A22;
C               joba = 1 :  reduce A22 to an SVD-like form.
C               joba = 2 :  reduce A22 to an upper trapezoidal form.
C            flag(2) = compq : see above.
C            flag(3) = compz : see above.
C            flag(4) = tol  : the tolerance to be used in determining
C            the rank of E and of A22. If tol > 0, then the given value
C            of tol is used as a lower bound for the reciprocal
C            condition numbers of leading submatrices of R or R22 in
C            the QR decompositions E*P = Q*R of E or A22*P22 = Q22*R22
C            of A22. A submatrix whose estimated condition number is
C            less than 1/tol is considered to be of full rank.
C            If tol <= 0, the default tolerance toldef = eps*l*n is
C            used instead, where eps is the machine precision. tol < 1.
C            Default :  flag = [ 0; 0; 0; 0 ].
C
C   Q1     - (optional) if compq = 2 the nl-by-nl given orthogonal
C            matrix Q1.
C   Z1     - (optional) if compz = 2 the n-by-n given orthogonal
C            matrix Z1.
C
C Output parameters:
C   Ao     - If task = 0, the l-by-n balanced matrix Dl*A*Dr.
C            If task = 1, the n-by-n upper Hessenberg matrix H = Q'*A*Z.
C            The elements below the first subdiagonal are set to zero.
C            If task = 2, the l-by-n transformed matrix Q'*A.
C            If task = 3, the l-by-n transformed matrix A*Z.
C            If task = 4, the l-by-n transformed matrix Q'*A*Z. If
C            joba = 1, this matrix is in the form
C
C                         ( A11  *   *  )
C                Q'*A*Z = (  *   Ar  0  ) ,
C                         (  *   0   0  )
C
C            where A11 is a ranke-by-ranke (see below) matrix and Ar
C            is a rnka22-by-rnka22 invertible diagonal matrix, with
C            decreasingly ordered positive diagonal elements.
C            If task = 5, the l-by-n transformed matrix Q'*A*Z. If
C            joba = 2, this matrix is in the form
C
C                         ( A11  *   *  )
C                Q'*A*Z = (  *   Ar  Z  ) ,
C                         (  *   0   0  )
C
C            where A11 is a ranke-by-ranke matrix and Ar is a
C            rnka22-by-rnka22 invertible upper triangular matrix.
C            If joba = 1 then A has the above form with Z = 0.
C            If task = 6, the n-by-n matrix Q'*A*Z in an upper
C            quasi-triangular form. The elements below the first
C            subdiagonal are set to zero.
C   Eo     - If task = 0, the l-by-n balanced matrix Dl*E*Dr.
C            If task = 1, the n-by-n upper triangular matrix T = Q'*E*Z.
C            The elements below the diagonal are set to zero.
C            If task = 2, the l-by-n transformed matrix Q'*E in upper
C            trapezoidal form, i.e.,
C
C                      ( E11 )
C               Q'*E = (     ) ,     if l >= n ,
C                      (  0  )
C            or
C
C               Q'*E = ( E11 E12 ),  if l < n ,
C
C            where E11 is a min(l,n)-by-min(l,n) upper triangular
C            matrix.
C            If task = 3, the l-by-n transformed matrix E*Z in upper
C            trapezoidal form, i.e.,
C
C                     ( E11 )
C               E*Z = (     ) ,      if l >= n ,
C                     (  R  )
C            or
C
C               E*Z = ( 0  R ),      if l < n ,
C
C            where R is a min(l,n)-by-min(l,n) upper triangular matrix.
C            If task = 4, the l-by-n transformed matrix Q'*E*Z in (2),
C            where Er is a ranke-by-ranke invertible diagonal matrix
C            having on the diagonal the decreasingly ordered positive
C            singular values of E.
C            If task = 5, the l-by-n transformed matrix Q'*E*Z in (2),
C            where Er is a ranke-by-ranke upper triangular invertible
C            matrix.
C            If task = 6, the n-by-n matrix Q'*E*Z in an upper
C            triangular form. The elements below the diagonal are set
C            to zero.
C   Bo     - If task = 0, the l-by-m balanced matrix Dl*B.
C            If task = 1, task = 2, or task >= 4, the nl-by-m
C            transformed matrix Q'*B.
C   Co     - If task = 0, the p-by-n balanced matrix C*Dr.
C            If task = 1, or task >= 3, the p-by-n transformed
C            matrix C*Z.
C   lscal  - the l-vector of scaling factors applied to A, E (and B)
C            from the left. If Dl(j) is the scaling factor applied to
C            row j, then lscal(j) = Dl(j), for j = 1,...,l.
C   rscal  - the n-vector of scaling factors applied to A, E (and C)
C            from the right. If Dr(j) is the scaling factor applied to
C            column j, then rscal(j) = Dr(j), for j = 1,...,n.
C   Q      - the nl-by-nl orthogonal matrix Q.
C            If compq = 1, Q' is the product of the transformations
C            which are applied to A, E, (and B) on the left.
C            If compq = 2, Q is the orthogonal matrix product Q1*Q.
C            If task = 6, Q contains the left orthogonal transformation
C            matrix used to reduce (A,E) to the real generalized Schur
C            form. (The columns of Q are the left generalized Schur
C            vectors of the pair (A,E).)
C   Z      - the n-by-n orthogonal matrix Z.
C            If compz = 1, Z is the product of the transformations
C            which are applied to A, E, (and C) on the right.
C            If compz = 2, Z is the orthogonal matrix product Z1*Z.
C            If task = 6, Z contains the right orthogonal transformation
C            matrix used to reduce (A,E) to the real generalized Schur
C            form. (The columns of Z are the right generalized Schur
C            vectors of the pair (A,E).)
C   ranks  - (optional) real vector with dimension at most 2 containing
C            the ranks of Er and Ar (if joba > 0). Specifically:
C            ranks(1) contains ranke, the effective (if task = 4), or
C            estimated (if task = 5) rank of matrix E, and thus also the
C            order of the invertible diagonal (if task = 4), or upper
C            triangular (if task = 5) submatrix Er. If task = 4, it is
C            computed as the number of singular values of E greater than
C            tol*svemax, with svemax the maximum singular value of E.
C            If joba > 0, ranks(2) contains rnka22, the effective
C            (if task = 4), or estimated (if task = 5) rank of matrix
C            A22, and thus also the order of the invertible diagonal
C            (if task = 4), or upper triangular (if task = 5) submatrix
C            Ar. If task = 4, it is computed as the number of singular
C            values of A22 greater than tol*svamax, where svamax is an
C            estimate of the maximum singular value of A.
C            If joba = 0, then rnka22 is not returned.
C   ev     - (optional) if task = 6, complex vector of length n.
C            The "alpha" part of the generalized eigenvalues.
C   db     - (optional) if task = 6, real vector of length n containing
C            the "beta" part of the generalized eigenvalues, i.e., the
C            generalized eigenvalues are ev(k)/db(k), k=1,...,n.
C            ev(k) and db(k), k=1,...,n, are the diagonals of the
C            complex Schur form that would result if the 2-by-2 diagonal
C            blocks of the real Schur form of (A,E) were further reduced
C            to triangular form using 2-by-2 complex unitary
C            transformations.
C
C Method
C   For task = 1, a QR factorization of E is computed and the
C   transformations are applied to A, B, and possibly Q1. Then, A is
C   reduced to upper Hessenberg form, preserving E triangular, by
C   an unblocked reduction, using two sequences of plane rotations
C   applied alternately from the left and from the right. The
C   corresponding transformations may be accumulated and/or applied
C   to the matrices B, C, and possibly Z1. If jobe = 1, the initial
C   reduction of E to upper triangular form is skipped.
C
C   For task = 2, a QR factorization of E is computed to reduce it
C   to the upper trapezoidal form. The transformations are also applied
C   to the rest of system matrices, A <- Q'*A ,  B <- Q'*B.
C
C   For task = 3, an RQ factorization of E is computed to reduce it
C   the upper trapezoidal form. The transformations are also applied
C   to the rest of system matrices, A <- A*Z,  C <- C*Z.
C
C   For task = 4, the singular value decomposition (SVD) of E is
C   computed, in the form
C
C                  ( Er  0 )
C         E  = Q * (       ) * Z'
C                  (  0  0 )
C
C   and the largest ranke-by-ranke leading diagonal submatrix Er is
C   found, whose condition number is less than 1/tol. Hence, ranke
C   defines the effective rank of matrix E.
C   If joba = 1, the same reduction is performed on A22 in the
C   partitioned matrix Q'*A*Z in (2), to obtain the transformed A22 in
C   the form (3), with Ar an invertible diagonal matrix. The accumulated
C   transformations are also applied to the rest of matrices,
C   B <- Q'*B,  C <- C*Z.
C
C   For task = 5, a truncated QR factorization of E, with column
C   pivoting, is computed in the form
C
C                     ( E11 E12 )
C         E * P = Q * (         )
C                     (  0  E22 )
C
C   and the largest ranke-by-ranke leading submatrix E11 whose
C   estimated condition number is less than 1/tol is found. Hence,
C   ranke defines the effective rank of matrix E. Further E22, being
C   negligible, is set to zero, and an orthogonal matrix Y is determined
C   such that
C
C         ( E11 E12 ) = ( Er  0 ) * Y .
C
C   The overal transformation matrix Z results as Z = P*Y and the
C   resulting transformed matrices Q'*A*Z and Q'*E*Z have the form (2),
C   where Er is an upper triangular invertible matrix.
C   If joba = 1 the same reduction is performed on A22 to obtain it
C   in the form (3), with Ar an upper triangular invertible matrix.
C   If joba = 2 then A22 is row compressed using a QR factorization
C   with column pivoting to the form (3), with Ar an upper triangular
C   invertible matrix. The transformations are also applied to the rest
C   of system matrices  B <- Q'*B, C <- C*Z.
C
C   If task = 6, the pair (A,E) is reduced to a real generalized Schur
C   form using an orthogonal equivalence transformation
C   (A,E) <- (Q'*A*Z,Q'*E*Z) and the transformation is applied to the
C   matrices B and C: B <- Q'*B and C <- C*Z.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2005,
C   Apr. 2009, Dec. 2010, July 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
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
      CHARACTER         COMPQ, COMPZ, JOB, JOBA, JOBE
      INTEGER           IHI, ILO, INFO, L, LDA, LDB, LDC, LDE, LDQ,
     $                  LDWORK, LDZ, M, N, P, RANKE, RNKA22
      DOUBLE PRECISION  THRESH, TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      mwPointer         A, ALPHAI, ALPHAR, B, BETA, C, DWORK, E, EV,
     $                  LSCALE, Q, RSCALE, Z
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           WITHQ, WITHQ1, WITHZ, WITHZ1
      INTEGER           FLAG(4), I, ICOMPQ, ICOMPZ, IJOB, IP, ISIZE,
     $                  ITASK, LFLAG, LIWORK, NC, NQ, NZ
      DOUBLE PRECISION  FLAGR(4), TEMP
C
C .. External subroutines ..
      EXTERNAL          DGGBAK, DGGBAL, TG01AD, TG01BD, TG01CD, TG01DD,
     $                  TG01ED, TG01FD, TG01WD
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'GSYSTRA requires at least 3 input arguments' )
      ELSE IF ( NLHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'GSYSTRA requires at least 2 output arguments' )
      END IF
C
C   task, A(nlxn), E(nlxn)(, B(nl,m), C(p,n), flag, Q1(nl,nl), Z1(n,n))
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
      IF ( ITASK.LT.0 .OR. ITASK.GT.6 ) THEN
         CALL mexErrMsgTxt
     $     ( 'TASK has 0, 1, ..., or 6 the only admissible values' )
      END IF
C
      IF ( ITASK.EQ.0 .OR. ITASK.EQ.4 ) THEN
         LFLAG = 2
      ELSE IF ( ITASK.EQ.1 ) THEN
         LFLAG = 3
      ELSE IF ( ITASK.EQ.2 .OR. ITASK.EQ.3 ) THEN
         LFLAG = 1
      ELSE IF ( ITASK.EQ.5 ) THEN
         LFLAG = 4
      ELSE
         LFLAG = 0
      END IF
C
      L = mxGetM( PRHS(2) )
      N = mxGetN( PRHS(2) )
      IF ( ( ITASK.EQ.1 .OR. ITASK.EQ.6 ) .AND. L.NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(3) ).NE.L .OR.
     $     mxGetN( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'E must have the same size as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'E must be a real matrix' )
      END IF
C
      IP = 4
      M  = 0
      IF ( NRHS.GE.IP ) THEN
         IF ( ITASK.NE.3 ) THEN
            M =  mxGetN( PRHS(IP) )
            IF ( mxGetM( PRHS(IP) ).NE.L ) THEN
               CALL mexErrMsgTxt( 'B must have the same rowsize as A' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'B must be a real matrix' )
            END IF
            IP = IP + 1
         END IF
      END IF
C
      P = 0
      IF ( NRHS.GE.IP ) THEN
         IF ( ITASK.NE.2 ) THEN
            P =  mxGetM( PRHS(IP) )
            IF ( mxGetN( PRHS(IP) ).NE.N ) THEN
               CALL mexErrMsgTxt
     $              ( 'C must have the same columnsize as A' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'C must be a real matrix' )
            END IF
            IP = IP + 1
         END IF
      END IF
C
C   flag
C
      DO 10 I = 1, LFLAG
         FLAG(I)  = 0
         FLAGR(I) = ZERO
   10 CONTINUE
C
      IF ( NRHS.GE.IP ) THEN
         ISIZE = mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) )
         IF ( ISIZE.GT.LFLAG ) THEN
            WRITE( TEXT, '(''FLAG must have at most '', I2,
     $                     '' elements'')' ) LFLAG
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'FLAG must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), FLAGR, ISIZE )
         IF ( ISIZE.GT.0 ) FLAG(1) = FLAGR(1)
         IF ( ISIZE.GT.1 ) FLAG(2) = FLAGR(2)
         IF ( ISIZE.GT.2 ) FLAG(3) = FLAGR(3)
         IP = IP + 1
      END IF
C
      ICOMPQ = 0
      ICOMPZ = 0
      TOL    = ZERO
C
      IF ( ITASK.NE.6 )
     $   IJOB = FLAG(1)
      IF ( ITASK.EQ.0 ) THEN
         IF ( IJOB.EQ.0 ) THEN
            JOB = 'A'
         ELSE IF ( IJOB.EQ.1 ) THEN
            JOB = 'B'
         ELSE IF ( IJOB.EQ.2 ) THEN
            JOB = 'C'
         ELSE IF ( IJOB.EQ.3 ) THEN
            JOB = 'N'
         ELSE
            CALL mexErrMsgTxt
     $        ( 'JOB has 0, 1, 2, or 3 the only admissible values' )
         END IF
         IF ( FLAGR(2).GE.ZERO ) THEN
            THRESH = FLAGR(2)
         ELSE
            CALL mexErrMsgTxt( 'THRESH must be non-negative' )
         END IF
C
      ELSE IF ( ITASK.EQ.1 ) THEN
         IF ( IJOB.EQ.0 ) THEN
            JOBE = 'G'
         ELSE IF ( IJOB.EQ.1 ) THEN
            JOBE = 'U'
         ELSE
            CALL mexErrMsgTxt
     $        ( 'JOBE has 0 or 1 the only admissible values' )
         END IF
         ICOMPQ = FLAG(2)
         ICOMPZ = FLAG(3)
C
      ELSE IF ( ITASK.EQ.2 ) THEN
         ICOMPQ = IJOB
C
      ELSE IF ( ITASK.EQ.3 ) THEN
         ICOMPZ = IJOB
C
      ELSE IF ( ITASK.EQ.4 ) THEN
         IF ( IJOB.EQ.0 ) THEN
            JOBA = 'N'
         ELSE IF ( IJOB.EQ.1 ) THEN
            JOBA = 'R'
         ELSE
            CALL mexErrMsgTxt
     $        ( 'JOBA has 0 or 1 the only admissible values' )
         END IF
         TOL = FLAGR(2)
C
      ELSE IF ( ITASK.EQ.5 ) THEN
         IF ( IJOB.EQ.0 ) THEN
            JOBA = 'N'
         ELSE IF ( IJOB.EQ.1 ) THEN
            JOBA = 'R'
         ELSE IF ( IJOB.EQ.2 ) THEN
            JOBA = 'T'
         ELSE
            CALL mexErrMsgTxt
     $        ( 'JOBA has 0, 1, or 2 the only admissible values' )
         END IF
         ICOMPQ = FLAG(2)
         ICOMPZ = FLAG(3)
         TOL    = FLAGR(4)
      END IF
C
      WITHQ1 = ITASK.EQ.1 .OR. ITASK.EQ.2 .OR. ITASK.EQ.5
      WITHZ1 = ITASK.EQ.1 .OR. ITASK.EQ.3 .OR. ITASK.EQ.5
C
      IF ( WITHQ1 ) THEN
         IF ( ICOMPQ.EQ.0 ) THEN
            COMPQ = 'N'
         ELSE IF ( ICOMPQ.EQ.1 ) THEN
            COMPQ = 'I'
         ELSE IF ( ICOMPQ.EQ.2 ) THEN
            IF ( ITASK.EQ.1 ) THEN
               COMPQ = 'V'
            ELSE
               COMPQ = 'U'
            END IF
         ELSE
            CALL mexErrMsgTxt
     $        ( 'COMPQ has 0, 1, or 2 the only admissible values' )
         END IF
      END IF
C
      IF ( WITHZ1 ) THEN
         IF ( ICOMPZ.EQ.0 ) THEN
            COMPZ = 'N'
         ELSE IF ( ICOMPZ.EQ.1 ) THEN
            COMPZ = 'I'
         ELSE IF ( ICOMPZ.EQ.2 ) THEN
            IF ( ITASK.EQ.1 ) THEN
               COMPZ = 'V'
            ELSE
               COMPZ = 'U'
            END IF
         ELSE
            CALL mexErrMsgTxt
     $        ( 'COMPZ has 0, 1, or 2 the only admissible values' )
         END IF
      END IF
C
      WITHQ  = ( ITASK.EQ.2 .AND. ICOMPQ.GE.1 .AND. NLHS.GE.4 ) .OR.
     $         ( ITASK.NE.2 .AND. ICOMPQ.GE.1 .AND. WITHQ1 .AND.
     $           NLHS.GE.5 ) .OR. ( ITASK.EQ.4 .OR. ITASK.EQ.6 )
      WITHZ  = ( ITASK.EQ.3 .AND. ICOMPZ.GE.1 .AND. NLHS.GE.4 ) .OR.
     $         ( ITASK.NE.3 .AND. ICOMPZ.GE.1 .AND. WITHZ1 .AND.
     $           NLHS.GE.4 + MIN( ICOMPQ, 1 ) ) .OR.
     $         ( ITASK.EQ.4 .OR. ITASK.EQ.6 )
C
      WITHQ1 = WITHQ1 .AND. ICOMPQ.EQ.2
      WITHZ1 = WITHZ1 .AND. ICOMPZ.EQ.2
C
C   Q1
C
      IF ( NRHS.GE.IP .AND. WITHQ1 ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.L .OR.
     $        mxGetN( PRHS(IP) ).NE.L ) THEN
            WRITE( TEXT,
     $          '(''Q1 must have '',I5,'' rows and columns'')' ) L
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'Q1 must be a real matrix' )
         END IF
         IP = IP + 1
      END IF
C
C   Z1
C
      IF ( NRHS.GE.IP .AND. WITHZ1 ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.N .OR.
     $        mxGetN( PRHS(IP) ).NE.N ) THEN
            WRITE( TEXT,
     $          '(''Z1 must have '',I5,'' rows and columns'')' ) N
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'Z1 must be a real matrix' )
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA = MAX( 1, L )
      LDE = LDA
      IF ( M.GT.0 .OR. ITASK.EQ.6 ) THEN
         LDB = LDA
      ELSE
         LDB = 1
      END IF
      LDC = MAX( 1, P )
      IF ( ITASK.NE.2 ) THEN
         NC = N
      ELSE
         NC = 0
      END IF
      IF ( WITHQ1 .OR. WITHQ ) THEN
         LDQ = LDA
         NQ  = L
      ELSE
         LDQ = 1
         NQ  = 0
      END IF
      IF ( WITHZ1 .OR. WITHZ ) THEN
         LDZ = MAX( 1, N )
         NZ  = N
      ELSE
         LDZ = 1
         NZ  = 0
      END IF
C
C   ldwork
C   Use larger values, for optimum performance in TG01BD, TG01CD,
C   TG01DD, and TG01WD.
C
      IF ( ITASK.EQ.0 ) THEN
         LDWORK = 3*( L + N )
      ELSE IF ( ITASK.EQ.1 ) THEN
         LDWORK = MAX( 1, 8*N )
         IF ( IJOB.EQ.0 )
     $      LDWORK = MAX( LDWORK, 3*N + MAX( N, M ) )
      ELSE IF ( ITASK.EQ.2 ) THEN
         LDWORK = MAX( 1, MIN( L, N ) + MAX( L, N, M ) )
      ELSE IF ( ITASK.EQ.3 ) THEN
         LDWORK = MAX( 1, MIN( L, N ) + MAX( L, N, P ) )
      ELSE IF ( ITASK.EQ.4 ) THEN
         LDWORK = MAX( 1, MIN( L, N ) +
     $                 MAX( 3*MIN( L, N ) + MAX( L, N ),
     $                      5*MIN( L, N ), M, P ) )
      ELSE IF ( ITASK.EQ.5 ) THEN
         LIWORK = N
         LDWORK = MAX( 1, P, MIN( L, N ) + MAX( 3*N, M, L ) )
      ELSE
         LDWORK = MAX( 7*( N + 1 ), 8*N ) + 16
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      A     = mxCalloc(  LDA*N,  8 )
      B     = mxCalloc(  LDB*M,  8 )
      C     = mxCalloc(  LDC*NC, 8 )
      DWORK = mxCalloc( LDWORK,  8 )
      E     = mxCalloc(  LDE*N,  8 )
      Q     = mxCalloc(  LDQ*NQ, 8 )
      Z     = mxCalloc(  LDZ*NZ, 8 )
      IF ( ITASK.EQ.0 ) THEN
         LSCALE = mxCalloc( L, 8 )
         RSCALE = mxCalloc( N, 8 )
      ELSE IF ( ITASK.EQ.1 ) THEN
         ALPHAI = mxCalloc( N, 8 )
         ALPHAR = mxCalloc( N, 8 )
      ELSE IF ( ITASK.EQ.5 ) THEN
         ALLOCATE ( IWORK( LIWORK ) )
      ELSE IF ( ITASK.EQ.6 ) THEN
         EV     = mxCalloc( N, 16 )
         ALPHAI = mxCalloc( N, 8 )
         ALPHAR = mxCalloc( N, 8 )
         BETA   = mxCalloc( N, 8 )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( A ), L*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL( E ), L*N )
      IP = 4
      IF ( NRHS.GE.IP .AND. ITASK.NE.3 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL( B ), L*M )
         IP = IP + 1
      END IF
      IF ( NRHS.GE.IP .AND. ITASK.NE.2 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL( C ), P*N )
         IP = IP + 1
      END IF
      IP = IP + 1
      IF ( NRHS.GE.IP .AND. WITHQ1 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL( Q ), L*L )
         IP = IP + 1
      END IF
      IF ( NRHS.GE.IP .AND. WITHZ1 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL( Z ), N*N )
      END IF
C
C Do the actual computations.
C
      IF ( ITASK.EQ.0 ) THEN
         CALL TG01AD( JOB, L, N, M, P, THRESH, %VAL( A ), LDA,
     $                %VAL( E ), LDE, %VAL( B ), LDB, %VAL( C ), LDC,
     $                %VAL( LSCALE ), %VAL( RSCALE ), %VAL( DWORK ),
     $                INFO )
      ELSE IF ( ITASK.EQ.1 ) THEN
         IF ( N.GT.0 ) THEN
            CALL DGGBAL( 'Permute', N, %VAL( A ), LDA, %VAL( E ), LDE,
     $                   ILO, IHI, %VAL( ALPHAI ), %VAL( ALPHAR ),
     $                   %VAL( DWORK ), INFO )
         ELSE
            ILO = 1
            IHI = 0
         END IF
         CALL TG01BD( JOBE, COMPQ, COMPZ, N, M, P, ILO, IHI, %VAL( A ),
     $                LDA, %VAL( E ), LDE, %VAL( B ), LDB, %VAL( C ),
     $                LDC, %VAL( Q ), LDQ, %VAL( Z ), LDZ,
     $                %VAL( DWORK ), LDWORK, INFO )
         IF ( N.GT.0 ) THEN
            IF( ICOMPQ.GT.0 )
     $         CALL DGGBAK( 'Permute', 'Left', N, ILO, IHI,
     $                      %VAL( ALPHAI ), %VAL( ALPHAR ), N,
     $                      %VAL( Q ), LDQ, INFO )
            IF( ICOMPZ.GT.0 )
     $         CALL DGGBAK( 'Permute', 'Right', N, ILO, IHI,
     $                      %VAL( ALPHAI ), %VAL( ALPHAR ), N,
     $                      %VAL( Z ), LDZ, INFO )
         END IF
      ELSE IF ( ITASK.EQ.2 ) THEN
         CALL TG01CD( COMPQ, L, N, M, %VAL( A ), LDA, %VAL( E ), LDE,
     $                %VAL( B ), LDB, %VAL( Q ), LDQ, %VAL( DWORK ),
     $                LDWORK, INFO )
      ELSE IF ( ITASK.EQ.3 ) THEN
         CALL TG01DD( COMPZ, L, N, P, %VAL( A ), LDA, %VAL( E ), LDE,
     $                %VAL( C ), LDC, %VAL( Z ), LDZ, %VAL( DWORK ),
     $                LDWORK, INFO )
      ELSE IF ( ITASK.EQ.4 ) THEN
         CALL TG01ED( JOBA, L, N, M, P, %VAL( A ), LDA, %VAL( E ), LDE,
     $                %VAL( B ), LDB, %VAL( C ), LDC, %VAL( Q ), LDQ,
     $                %VAL( Z ), LDZ, RANKE, RNKA22, TOL, %VAL( DWORK ),
     $                LDWORK, INFO )
      ELSE IF ( ITASK.EQ.5 ) THEN
         CALL TG01FD( COMPQ, COMPZ, JOBA, L, N, M, P, %VAL( A ), LDA,
     $                %VAL( E ), LDE, %VAL( B ), LDB, %VAL( C ), LDC,
     $                %VAL( Q ), LDQ, %VAL( Z ), LDZ, RANKE, RNKA22,
     $                TOL, IWORK, %VAL( DWORK ), LDWORK, INFO )
      ELSE IF ( ITASK.EQ.6 ) THEN
         CALL TG01WD( N, M, P, %VAL( A ), LDA, %VAL( E ), LDE,
     $                %VAL( B ), LDB, %VAL( C ), LDC, %VAL( Q ), LDQ,
     $                %VAL( Z ), LDZ, %VAL( ALPHAR ), %VAL( ALPHAI ),
     $                %VAL( BETA ), %VAL( DWORK ), LDWORK, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         IP = 3
         PLHS(1) = mxCreateDoubleMatrix( L, N, 0 )
         CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS(1) ), L*N )
         PLHS(2) = mxCreateDoubleMatrix( L, N, 0 )
         CALL mxCopyReal8ToPtr( %VAL( E ), mxGetPr( PLHS(2) ), L*N )
         IF ( NLHS.GE.IP .AND. ITASK.NE.3 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( L, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( B ), mxGetPr( PLHS(IP) ),
     $                             L*M )
            IP = IP + 1
         END IF
         IF ( NLHS.GE.IP .AND. ITASK.NE.2 ) THEN
            PLHS(IP) = mxCreateDoubleMatrix( P, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( C ), mxGetPr( PLHS(IP) ),
     $                             P*N )
            IP = IP + 1
         END IF
         IF ( NLHS.GE.IP ) THEN
            IF ( ITASK.EQ.0 ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( L, MIN( L, 1 ), 0 )
               CALL mxCopyReal8ToPtr( %VAL( LSCALE ),
     $                                mxGetPr( PLHS(IP) ), L*MIN( L, 1 )
     $                              )
               IP = IP + 1
               IF ( NLHS.GE.IP ) THEN
                  PLHS(IP) = mxCreateDoubleMatrix( N, MIN( N, 1 ),
     $                                             0 )
                  CALL mxCopyReal8ToPtr( %VAL( RSCALE ),
     $                                   mxGetPr( PLHS(IP) ),
     $                                   N*MIN( N, 1 ) )
               END IF
            ELSE
               IF ( WITHQ ) THEN
                  PLHS(IP) = mxCreateDoubleMatrix( L, L, 0 )
                  CALL mxCopyReal8ToPtr( %VAL( Q ),
     $                                   mxGetPr( PLHS(IP) ), L*L )
                  IP = IP + 1
               END IF
               IF ( NLHS.GE.IP .AND. WITHZ ) THEN
                  PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
                  CALL mxCopyReal8ToPtr( %VAL( Z ),
     $                                   mxGetPr( PLHS(IP) ), N*N )
                  IP = IP + 1
               END IF
               IF ( NLHS.GE.IP ) THEN
                  IF ( ITASK.EQ.4 .OR. ITASK.EQ.5 ) THEN
                     FLAGR(1) = RANKE
                     IF ( IJOB.EQ.0 ) THEN
                        ISIZE = 1
                     ELSE
                        ISIZE = 2
                        FLAGR(2) = RNKA22
                     END IF
                     PLHS(IP) = mxCreateDoubleMatrix( ISIZE, 1, 0 )
                     CALL mxCopyReal8ToPtr( FLAGR,
     $                                      mxGetPr( PLHS(IP) ),
     $                                      ISIZE )
                  ELSE IF ( ITASK.EQ.6 ) THEN
                     CALL mxCopyPtrToComplex16( ALPHAR, ALPHAI,
     $                                          %VAL( EV ),
     $                                          N*MIN( N, 1 ) )
                     PLHS(IP) = mxCreateDoubleMatrix( N, MIN( 1, N ),
     $                                                1 )
                     CALL mxCopyComplex16ToPtr( %VAL( EV ),
     $                                          mxGetPr( PLHS(IP) ),
     $                                          mxGetPi( PLHS(IP) ),
     $                                          N*MIN( N, 1 ) )
                     IP = IP + 1
                     PLHS(IP) = mxCreateDoubleMatrix( N, MIN( 1, N ),
     $                                                0 )
                     CALL mxCopyReal8ToPtr( %VAL( BETA ),
     $                                      mxGetPr( PLHS(IP) ),
     $                                      N*MIN( N, 1 ) )
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      CALL mxFree( A     )
      CALL mxFree( B     )
      CALL mxFree( C     )
      CALL mxFree( DWORK )
      CALL mxFree( E     )
      CALL mxFree( Q     )
      CALL mxFree( Z     )
      IF ( ITASK.EQ.0 ) THEN
         CALL mxFree( LSCALE )
         CALL mxFree( RSCALE )
      ELSE IF ( ITASK.EQ.1 ) THEN
         CALL mxFree( ALPHAI )
         CALL mxFree( ALPHAR )
      ELSE IF ( ITASK.EQ.5 ) THEN
         DEALLOCATE ( IWORK )
      ELSE IF ( ITASK.EQ.6 ) THEN
         CALL mxFree( EV )
         CALL mxFree( ALPHAI )
         CALL mxFree( ALPHAR )
         CALL mxFree( BETA   )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( ITASK.EQ.4 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM TG01ED'')' )
     $         INFO
         ELSE IF ( ITASK.EQ.6 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM TG01WD'')' )
     $         INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of GSYSTRA ***
      END
