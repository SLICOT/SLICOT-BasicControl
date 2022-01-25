% GSYSTRA.F  - MEX-function to perform various equivalence
%              transformations for descriptor systems with scaling,
%              generalized Schur form, etc., using SLICOT routines
%              TG01AD, TG01BD, TG01CD, TG01DD, TG01ED, TG01FD, and
%              TG01WD.
%
%   [(Ao,Eo(,Bo)(,Co)(lscal,rscal)(,Q)(,Z)(,ranks)(,ev,db)] =
%                           GSYSTRA(task,A,E(,B)(,C)(,flag)(,Q1,Z1))
%
%   [Ao,Eo,Bo,Co(,lscal,rscal)] = GSYSTRA(0,A,E,B,C(,flag))
%   [Ao,Eo,Bo,Co(,Q,Z)]         = GSYSTRA(1,A,E,B,C(,flag)(,Q1,Z1))
%   [Ao,Eo,Bo(,Q)]              = GSYSTRA(2,A,E,B(,flag)(,Q1))
%   [Ao,Eo,Co(,Z)]              = GSYSTRA(3,A,E,C(,flag)(,Z1))
%   [Ao,Eo,Bo,Co(,Q,Z)(,ranks)] = GSYSTRA(4,A,E,B,C(,flag))
%   [Ao,Eo,Bo,Co(,Q,Z)(,ranks)] = GSYSTRA(5,A,E,B,C(,flag)(,Q1,Z1))
%   [Ao,Eo,Bo,Co(,Q,Z)(,ev,db)] = GSYSTRA(6,A,E,B,C)
%
%   GSYSTRA performs one of the equivalence transformations specified by
%   the value of the parameter task (task = 0, 1, ..., 6), for a descriptor
%   triple (A-lambda E,B,C), or its parts:
%
%   0) To balance the matrices of the system pencil
%
%      S =  ( A  B ) - lambda ( E  0 ) .                             (1)
%           ( C  0 )          ( 0  0 )
%
%   Balancing involves diagonal equivalence transformations
%   (Dl*A*Dr - lambda Dl*E*Dr, Dl*B, C*Dr) applied to the system
%   (A-lambda E,B,C) to make the rows and columns of system pencil
%   matrices  diag(Dl,I)*S*diag(Dr,I)  as close in norm as possible.
%   Balancing may reduce the 1-norms of the matrices of the system
%   pencil S. The balancing can be performed optionally on the
%   following particular system pencils
%
%      S = A-lambda E,
%      S = ( A-lambda E  B ),    or
%      S = ( A-lambda E ).
%          (     C      )
%
%   1) To reduce the matrices A and E of the system pencil (1) to
%   generalized upper Hessenberg form using orthogonal transformations,
%
%      Q'*A*Z = H,   Q'*E*Z = T,
%
%   where H is upper Hessenberg, T is upper triangular, Q and Z
%   are orthogonal, and ' means transpose. The corresponding
%   transformations, written compactly as diag(Q',I)*S*diag(Z,I),
%   are also applied to B and C, getting Q'*B and C*Z.
%   The orthogonal matrices Q and Z are determined as products of
%   Givens rotations. They may either be formed explicitly, or they
%   may be postmultiplied into input matrices Q1 and Z1, so that
%
%      Q1*A*Z1' = (Q1*Q)*H*(Z1*Z)'
%      Q1*E*Z1' = (Q1*Q)*T*(Z1*Z)'.
%
%   2) To reduce the descriptor system pair (A-lambda E,B) to the
%   QR-coordinate form by computing an orthogonal transformation
%   matrix Q such that the transformed descriptor system pair
%   (Q'*A-lambda Q'*E, Q'*B) has the descriptor matrix Q'*E in an
%   upper trapezoidal form. The left orthogonal transformations
%   performed to reduce E can be optionally accumulated.
%
%   3) To reduce the descriptor system pair (C,A-lambda E) to the
%   RQ-coordinate form by computing an orthogonal transformation
%   matrix Z such that the transformed descriptor system pair
%   (C*Z,A*Z-lambda E*Z) has the descriptor matrix E*Z in an upper
%   trapezoidal form. The right orthogonal transformations performed
%   to reduce E can be optionally accumulated.
%
%   4) To compute for the descriptor system (A-lambda E,B,C) the
%   orthogonal transformation matrices Q and Z such that the
%   transformed system (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z) is in an
%   SVD (singular value decomposition) coordinate form with
%   the system matrices Q'*A*Z and Q'*E*Z in the form
%
%               ( A11  A12 )             ( Er  0 )
%      Q'*A*Z = (          ) ,  Q'*E*Z = (       ) ,                 (2)
%               ( A21  A22 )             (  0  0 )
%
%   where Er is an invertible diagonal matrix having on the diagonal
%   the decreasingly ordered nonzero singular values of E. Optionally,
%   the A22 matrix can be further reduced to the SVD form
%
%               ( Ar  0 )
%               (       ) ,                                          (3)
%               (  0  0 )
%
%   where Ar is an invertible diagonal matrix having on the diagonal
%   the decreasingly ordered nonzero singular values of A22. The left
%   and/or right orthogonal transformations performed to reduce E
%   and A22 are accumulated.
%
%   5) To compute for the descriptor system (A-lambda E,B,C) the
%   orthogonal transformation matrices Q and Z such that the
%   transformed system (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z) is in an
%   SVD-like coordinate form (2), where Er is an upper triangular
%   invertible matrix. Optionally, the A22 matrix can be further
%   reduced to the form
%
%               ( Ar  X )
%               (       ) ,
%               (  0  0 )
%
%   with Ar an upper triangular invertible matrix, and X either a full
%   or a zero matrix. The left and/or right orthogonal transformations
%   performed to reduce E and A22 can be optionally accumulated.
%
%   6) To reduce the pair (A,E) to a real generalized Schur form
%   by using an orthogonal equivalence transformation
%   (A,E) <-- (Q'*A*Z,Q'*E*Z) and to apply the transformation
%   to the matrices B and C: B <-- Q'*B and C <-- C*Z.
%
%   Description of input parameters:
%   task   - integer specifying the computations to be performed.
%            task = 0 :  balance the descriptor system;
%            task = 1 :  compute the generalized upper Hessenberg form;
%            task = 2 :  compute QR-coordinate form;
%            task = 3 :  compute RQ-coordinate form;
%            task = 4 :  compute SVD coordinate form;
%            task = 5 :  compute SVD-like coordinate form;
%            task = 6 :  compute real generalized Schur form.
%   A      - the nl-by-n state dynamics matrix A, where nl = n, if
%            task = 1 or task = 6, and nl = l, otherwise.
%   E      - the nl-by-n descriptor matrix E. If task = 1 and jobe = 1,
%            the matrix E is assumed upper triangular.
%   B      - the nl-by-m input/state matrix B.
%   C      - the p-by-n state/output matrix C.
%   flag   - (optional) real vector specifying various options,
%            depending on task.
%
%            For task = 0, flag has length 2, and contains:
%            flag(1) = job : indicates which matrices are involved in
%            balancing, as follows:
%               job = 0 :  All matrices are involved in balancing;
%               job = 1 :  B, A and E are involved in balancing;
%               job = 2 :  C, A and E are involved in balancing;
%               job = 3 :  B and C are not involved in balancing.
%            flag(2) = thresh >= 0 : threshold value for magnitude of
%            elements: elements with magnitude less than or equal to
%            thresh are ignored for balancing.
%            Default :  flag = [ 0; 0 ].
%
%            For task = 1, flag has length 3 (at most), and contains:
%            flag(1) = jobe : specifies whether E is a general square
%            or an upper triangular matrix, as follows:
%               jobe = 0 :  E is a general square matrix;
%               jobe = 1 :  E is an upper triangular matrix.
%            flag(2) = compq : indicates what should be done with
%            matrix Q, as follows:
%               compq = 0 :  do not compute Q;
%               compq = 1 :  Q is initialized to the unit matrix, and
%                            the orthogonal matrix Q is returned;
%               compq = 2 :  Q is initialized to an orthogonal matrix Q1
%                            and the product Q1*Q is returned.
%            flag(3) = compz : indicates what should be done with
%            matrix Z, as follows:
%               compz = 0 :  do not compute Z;
%               compz = 1 :  Z is initialized to the unit matrix, and
%                            the orthogonal matrix Z is returned;
%               compz = 2 :  Z is initialized to an orthogonal matrix Z1
%                            and the product Z1*Z is returned.
%            Default :  flag = [ 0; 0; 0 ].
%
%            For task = 2, flag has length 1, and contains compq above.
%            Default :  flag = 0.
%
%            For task = 3, flag has length 1, and contains compz above.
%            Default :  flag = 0.
%
%            For task = 4, flag has length 2, and contains:
%            flag(1) = joba : specifies whether or not A22 should be
%            reduced, as follows:
%               joba = 0 :  do not reduce A22;
%               joba = 1 :  reduce A22 to an SVD form.
%            flag(2) = tol  : the tolerance to be used in determining
%            the rank of E and of A22. If tol > 0, then singular values
%            less than tol*svmax are treated as zero, where svmax is the
%            maximum singular value of E or of its estimate for A.
%            If tol <= 0, the default tolerance toldef = eps*l*n is used
%            instead, where eps is the machine precision. tol < 1.
%            Default :  flag = [ 0; 0 ].
%
%            For task = 5, flag has length 4, and contains:
%            flag(1) = joba : specifies whether or not A22 should be
%            reduced, as follows:
%               joba = 0 :  do not reduce A22;
%               joba = 1 :  reduce A22 to an SVD-like form.
%               joba = 2 :  reduce A22 to an upper trapezoidal form.
%            flag(2) = compq : see above.
%            flag(3) = compz : see above.
%            flag(4) = tol  : the tolerance to be used in determining
%            the rank of E and of A22. If tol > 0, then the given value
%            of tol is used as a lower bound for the reciprocal
%            condition numbers of leading submatrices of R or R22 in
%            the QR decompositions E*P = Q*R of E or A22*P22 = Q22*R22
%            of A22. A submatrix whose estimated condition number is
%            less than 1/tol is considered to be of full rank.
%            If tol <= 0, the default tolerance toldef = eps*l*n is
%            used instead, where eps is the machine precision. tol < 1.
%            Default :  flag = [ 0; 0; 0; 0 ].
%
%   Q1     - (optional) if compq = 1 the nl-by-nl given orthogonal
%            matrix Q1.
%   Z1     - (optional) if compz = 1 the n-by-n given orthogonal
%            matrix Z1.
%
%   Description of output parameters:
%   Ao     - If task = 0, the l-by-n balanced matrix Dl*A*Dr.
%            If task = 1, the n-by-n upper Hessenberg matrix H = Q'*A*Z.
%            The elements below the first subdiagonal are set to zero.
%            If task = 2, the l-by-n transformed matrix Q'*A.
%            If task = 3, the l-by-n transformed matrix A*Z.
%            If task = 4, the l-by-n transformed matrix Q'*A*Z. If
%            joba = 1, this matrix is in the form
%
%                         ( A11  *   *  )
%                Q'*A*Z = (  *   Ar  0  ) ,
%                         (  *   0   0  )
%
%            where A11 is a ranke-by-ranke (see below) matrix and Ar
%            is a rnka22-by-rnka22 invertible diagonal matrix, with
%            decreasingly ordered positive diagonal elements.
%            If task = 5, the l-by-n transformed matrix Q'*A*Z. If
%            joba = 2, this matrix is in the form
%
%                         ( A11  *   *  )
%                Q'*A*Z = (  *   Ar  Z  ) ,
%                         (  *   0   0  )
%
%            where A11 is a ranke-by-ranke matrix and Ar is a
%            rnka22-by-rnka22 invertible upper triangular matrix.
%            If joba = 1 then A has the above form with Z = 0.
%            If task = 6, the n-by-n matrix Q'*A*Z in an upper
%            quasi-triangular form. The elements below the first
%            subdiagonal are set to zero.
%   Eo     - If task = 0, the l-by-n balanced matrix Dl*E*Dr.
%            If task = 1, the n-by-n upper triangular matrix T = Q'*E*Z.
%            The elements below the diagonal are set to zero.
%            If task = 2, the l-by-n transformed matrix Q'*E in upper
%            trapezoidal form, i.e.,
%
%                      ( E11 )
%               Q'*E = (     ) ,     if l >= n ,
%                      (  0  )
%            or
%
%               Q'*E = ( E11 E12 ),  if l < n ,
%
%            where E11 is a min(l,n)-by-min(l,n) upper triangular
%            matrix.
%            If task = 3, the l-by-n transformed matrix E*Z in upper
%            trapezoidal form, i.e.,
%
%                     ( E11 )
%               E*Z = (     ) ,      if l >= n ,
%                     (  R  )
%            or
%
%               E*Z = ( 0  R ),      if l < n ,
%
%            where R is a min(l,n)-by-min(l,n) upper triangular matrix.
%            If task = 4, the l-by-n transformed matrix Q'*E*Z in (2),
%            where Er is a ranke-by-ranke invertible diagonal matrix
%            having on the diagonal the decreasingly ordered positive
%            singular values of E.
%            If task = 5, the l-by-n transformed matrix Q'*E*Z in (2),
%            where Er is a ranke-by-ranke upper triangular invertible
%            matrix.
%            If task = 6, the n-by-n matrix Q'*E*Z in an upper
%            triangular form. The elements below the diagonal are set
%            to zero.
%   Bo     - If task = 0, the l-by-m balanced matrix Dl*B.
%            If task = 1, task = 2, or task >= 4, the nl-by-m
%            transformed matrix Q'*B.
%   Co     - If task = 0, the p-by-n balanced matrix C*Dr.
%            If task = 1, or task >= 3, the p-by-n transformed
%            matrix C*Z.
%   lscal  - the l-vector of scaling factors applied to A, E (and B)
%            from the left. If Dl(j) is the scaling factor applied to
%            row j, then lscal(j) = Dl(j), for j = 1,...,l.
%   rscal  - the n-vector of scaling factors applied to A, E (and C)
%            from the right. If Dr(j) is the scaling factor applied to
%            column j, then rscal(j) = Dr(j), for j = 1,...,n.
%   Q      - the nl-by-nl orthogonal matrix Q.
%            If compq = 1, Q' is the product of the transformations
%            which are applied to A, E, (and B) on the left.
%            If compq = 2, Q is the orthogonal matrix product Q1*Q.
%            If task = 6, Q contains the left orthogonal transformation
%            matrix used to reduce (A,E) to the real generalized Schur
%            form. (The columns of Q are the left generalized Schur
%            vectors of the pair (A,E).)
%   Z      - the n-by-n orthogonal matrix Z.
%            If compz = 1, Z is the product of the transformations
%            which are applied to A, E, (and C) on the right.
%            If compz = 2, Z is the orthogonal matrix product Z1*Z.
%            If task = 6, Z contains the right orthogonal transformation
%            matrix used to reduce (A,E) to the real generalized Schur
%            form. (The columns of Z are the right generalized Schur
%            vectors of the pair (A,E).)
%   ranks  - (optional) real vector with dimension at most 2 containing
%            the ranks of Er and Ar (if joba > 0). Specifically:
%            ranks(1) contains ranke, the effective (if task = 4), or
%            estimated (if task = 5) rank of matrix E, and thus also the
%            order of the invertible diagonal (if task = 4), or upper
%            triangular (if task = 5) submatrix Er. If task = 4, it is
%            computed as the number of singular values of E greater than
%            tol*svemax, with svemax the maximum singular value of E.
%            If joba > 0, ranks(2) contains rnka22, the effective
%            (if task = 4), or estimated (if task = 5) rank of matrix
%            A22, and thus also the order of the invertible diagonal
%            (if task = 4), or upper triangular (if task = 5) submatrix
%            Ar. If task = 4, it is computed as the number of singular
%            values of A22 greater than tol*svamax, where svamax is an
%            estimate of the maximum singular value of A.
%            If joba = 0, then rnka22 is not returned.
%   ev     - (optional) if task = 6, complex vector of length n.
%            The "alpha" part of the generalized eigenvalues.
%   db     - (optional) if task = 6, real vector of length n containing
%            the "beta" part of the generalized eigenvalues, i.e., the
%            generalized eigenvalues are ev(k)/db(k), k=1,...,n.
%            ev(k) and db(k), k=1,...,n, are the diagonals of the
%            complex Schur form that would result if the 2-by-2 diagonal
%            blocks of the real Schur form of (A,E) were further reduced
%            to triangular form using 2-by-2 complex unitary
%            transformations.
%
%   Method
%   For task = 1, a QR factorization of E is computed and the
%   transformations are applied to A, B, and possibly Q1. Then, A is
%   reduced to upper Hessenberg form, preserving E triangular, by
%   an unblocked reduction, using two sequences of plane rotations
%   applied alternately from the left and from the right. The
%   corresponding transformations may be accumulated and/or applied
%   to the matrices B, C, and possibly Z1. If jobe = 1, the initial
%   reduction of E to upper triangular form is skipped.
%
%   For task = 2, a QR factorization of E is computed to reduce it
%   to the upper trapezoidal form. The transformations are also applied
%   to the rest of system matrices, A <- Q'*A ,  B <- Q'*B.
%
%   For task = 3, an RQ factorization of E is computed to reduce it
%   the upper trapezoidal form. The transformations are also applied
%   to the rest of system matrices, A <- A*Z,  C <- C*Z.
%
%   For task = 4, the singular value decomposition (SVD) of E is
%   computed, in the form
%
%                  ( Er  0 )
%         E  = Q * (       ) * Z'
%                  (  0  0 )
%
%   and the largest ranke-by-ranke leading diagonal submatrix Er is
%   found, whose condition number is less than 1/tol. Hence, ranke
%   defines the effective rank of matrix E.
%   If joba = 1, the same reduction is performed on A22 in the
%   partitioned matrix Q'*A*Z in (2), to obtain the transformed A22 in
%   the form (3), with Ar an invertible diagonal matrix. The accumulated
%   transformations are also applied to the rest of matrices,
%   B <- Q'*B,  C <- C*Z.
%
%   For task = 5, a truncated QR factorization of E, with column
%   pivoting, is computed in the form
%
%                     ( E11 E12 )
%         E * P = Q * (         )
%                     (  0  E22 )
%
%   and the largest ranke-by-ranke leading submatrix E11 whose
%   estimated condition number is less than 1/tol is found. Hence,
%   ranke defines the effective rank of matrix E. Further E22, being
%   negligible, is set to zero, and an orthogonal matrix Y is determined
%   such that
%
%         ( E11 E12 ) = ( Er  0 ) * Y .
%
%   The overal transformation matrix Z results as Z = P*Y and the
%   resulting transformed matrices Q'*A*Z and Q'*E*Z have the form (2),
%   where Er is an upper triangular invertible matrix.
%   If joba = 1 the same reduction is performed on A22 to obtain it
%   in the form (3), with Ar an upper triangular invertible matrix.
%   If joba = 2 then A22 is row compressed using a QR factorization
%   with column pivoting to the form (3), with Ar an upper triangular
%   invertible matrix. The transformations are also applied to the rest
%   of system matrices  B <- Q'*B, C <- C*Z.
%
%   If task = 6, the pair (A,E) is reduced to a real generalized Schur
%   form using an orthogonal equivalence transformation
%   (A,E) <- (Q'*A*Z,Q'*E*Z) and the transformation is applied to the
%   matrices B and C: B <- Q'*B and C <- C*Z.
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
