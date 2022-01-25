#include "fintrf.h"
C GSYSCOM.F  - Gateway function to transform a descriptor system, by
C              equivalence transformations, to a controllable or
C              observable staircase form, or to a reduced (controllable,
C              observable, or irreducible) form, using SLICOT routines
C              TG01HD, TG01ID, and TG01JD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [Ao,Eo,Bo,Co(,Q,Z)(,orders,sizes)] =
C                       GSYSCOM(task,A,E,B,C(,flag)(,Q1,Z1))
C
C   [Ao,Eo,Bo,Co(,Q,Z)(,orders,sizes)] =
C                       GSYSCOM(1,A,E,B,C(,flag)(,Q1,Z1))
C   [Ao,Eo,Bo,Co(,Q,Z)(,orders,sizes)] =
C                       GSYSCOM(2,A,E,B,C(,flag)(,Q1,Z1))
C   [Ao,Eo,Bo,Co(,infred,sizes)] =
C                       GSYSCOM(3,A,E,B,C(,flag))
C
C Purpose:
C   To perform one of the equivalence transformations specified by the
C   value of the parameter task (task = 1, 2, or 3), for a descriptor
C   triple (A-lambda E,B,C), to reduce it to a staircase form, or to a
C   controllable, observable, or irreducible form:
C
C   1) To compute orthogonal transformation matrices Q and Z which
C   reduce the n-th order descriptor system (A-lambda*E,B,C) to the form
C
C               ( Ac  *  )             ( Ec  *  )           ( Bc )
C      Q'*A*Z = (        ) ,  Q'*E*Z = (        ) ,  Q'*B = (    ) ,
C               ( 0  Anc )             ( 0  Enc )           ( 0  )
C
C         C*Z = ( Cc Cnc ) ,                                         (1)
C
C   where the c-th order descriptor system (Ac-lambda*Ec,Bc,Cc) is
C   finite and/or infinite controllable. The pencil Anc - lambda*Enc
C   is regular of order n-c and contains the uncontrollable finite
C   and/or infinite eigenvalues of the pencil A-lambda*E. The reduced
C   order descriptor system (Ac-lambda*Ec,Bc,Cc) has the same
C   transfer-function matrix as the original system (A-lambda*E,B,C).
C   The left and/or right orthogonal transformations Q and Z performed
C   to reduce the system matrices can be optionally accumulated.
C
C   For jobcon = 0 or 2 (see parameter flag below), the pencil
C   ( Bc Ec-lambda*Ac ) has full row rank c for all finite lambda and
C   is in a staircase form with
C                    _      _          _        _
C                  ( E1,0   E1,1  ...  E1,k-1   E1,k  )
C                  (        _          _        _     )
C      ( Bc Ec ) = (  0     E2,1  ...  E2,k-1   E2,k  ) ,            (2)
C                  (              ...  _        _     )
C                  (  0       0   ...  Ek,k-1   Ek,k  )
C
C                    _          _        _
C                  ( A1,1  ...  A1,k-1   A1,k  )
C                  (            _        _     )
C        Ac      = (   0   ...  A2,k-1   A2,k  ) ,                   (3)
C                  (       ...           _     )
C                  (   0   ...    0      Ak,k  )
C         _
C   where Ei,i-1 is an rtau(i)-by-rtau(i-1) full row rank matrix
C                          _
C   (with rtau(0) = m) and Ai,i is an rtau(i)-by-rtau(i) upper
C   triangular matrix.
C
C   For jobcon = 1, the pencil ( Bc Ac-lambda*Ec ) has full row
C   rank c for all finite lambda and is in a staircase form with
C                    _     _          _        _
C                  ( A1,0  A1,1  ...  A1,k-1   A1,k  )
C                  (       _          _        _     )
C      ( Bc Ac ) = (  0    A2,1  ...  A2,k-1   A2,k  ) ,             (4)
C                  (             ...  _        _     )
C                  (  0      0   ...  Ak,k-1   Ak,k  )
C
C                    _          _        _
C                  ( E1,1  ...  E1,k-1   E1,k  )
C                  (            _        _     )
C        Ec      = (   0   ...  E2,k-1   E2,k  ) ,                   (5)
C                  (       ...           _     )
C                  (   0   ...    0      Ek,k  )
C         _
C   where Ai,i-1 is an rtau(i)-by-rtau(i-1) full row rank matrix
C                          _
C   (with rtau(0) = m) and Ei,i is an rtau(i)-by-rtau(i) upper
C   triangular matrix.
C
C   For jobcon = 0, the (n-c)-by-(n-c) regular pencil Anc - lambda*Enc
C   has the form
C
C                         ( Ainc - lambda*Einc         *          )
C      Anc - lambda*Enc = (                                       ) ,
C                         (        0           Afnc - lambda*Efnc )
C
C   where:
C     1) the inc-by-inc regular pencil Ainc - lambda*Einc, with Ainc
C        upper triangular and nonsingular, contains the uncontrollable
C        infinite eigenvalues of A - lambda*E;
C     2) the (n-c-inc)-by-(n-c-inc) regular pencil Afnc - lambda*Efnc,
C        with Efnc upper triangular and nonsingular, contains the
C        uncontrollable finite eigenvalues of A - lambda*E.
C
C   Note: The significance of the two diagonal blocks can be
C         interchanged by calling the gateway with the arguments A and E
C         interchanged. In this case, Ainc - lambda*Einc contains the
C         uncontrollable zero eigenvalues of A - lambda*E, while
C         Afnc - lambda*Efnc contains the uncontrollable nonzero finite
C         and infinite eigenvalues of A - lambda*E.
C
C   For jobcon = 1, the pencil Anc - lambda*Enc has the form
C
C      Anc - lambda*Enc = Afnc - lambda*Efnc ,
C
C   where the regular pencil Afnc - lambda*Efnc, with Efnc upper
C   triangular and nonsingular, contains the uncontrollable finite
C   eigenvalues of A - lambda*E.
C
C   For jobcon = 2, the pencil Anc - lambda*Enc has the form
C
C      Anc - lambda*Enc = Ainc - lambda*Einc ,
C
C   where the regular pencil Ainc - lambda*Einc, with Ainc upper
C   triangular and nonsingular, contains the uncontrollable nonzero
C   finite and infinite eigenvalues of A - lambda*E.
C
C   2) To compute orthogonal transformation matrices Q and Z which
C   reduce the n-th order descriptor system (A-lambda*E,B,C) to the form
C
C               ( Ano  * )             ( Eno  * )           ( Bno )
C      Q'*A*Z = (        ) ,  Q'*E*Z = (        ) ,  Q'*B = (     ) ,
C               ( 0   Ao )             ( 0   Eo )           ( Bo  )
C
C         C*Z = ( 0   Co ) ,                                         (6)
C
C   where the o-th order descriptor system (Ao-lambda*Eo,Bo,Co) is a
C   finite and/or infinite observable. The pencil Ano - lambda*Eno is
C   regular of order n-o and contains the unobservable finite and/or
C   infinite eigenvalues of the pencil A-lambda*E. The reduced order
C   descriptor system (Ao-lambda*Eo,Bo,Co) has the same
C   transfer-function matrix as the original system (A-lambda*E,B,C).
C   The left and/or right orthogonal transformations Q and Z performed
C   to reduce the system matrices can be optionally accumulated.
C
C   For jobobs = 0 or 2, the pencil ( Eo-lambda*Ao ) has full column
C                                   (      Co      )
C   rank o for all finite lambda and is in a staircase form with
C                    _      _            _      _
C                  ( Ek,k   Ek,k-1   ... Ek,2   Ek,1   )
C                  ( _      _            _      _      )
C      ( Eo )  =   ( Ek-1,k Ek-1,k-1 ... Ek-1,2 Ek-1,1 ) ,           (7)
C      ( Co )      (     ...         ... _      _      )
C                  (  0       0      ... E1,2   E1,1   )
C                  (                            _      )
C                  (  0       0      ... 0      E0,1   )
C                    _          _      _
C                  ( Ak,k  ...  Ak,2   Ak,1 )
C                  (       ...  _      _    )
C        Ao    =   (   0   ...  A2,2   A2,1 ) ,                      (8)
C                  (                   _    )
C                  (   0   ...    0    A1,1 )
C         _
C   where Ei-1,i is a ctau(i-1)-by-ctau(i) full column rank matrix
C                          _
C   (with ctau(0) = p) and Ai,i is a ctau(i)-by-ctau(i) upper
C   triangular matrix.
C
C   For jobobs = 1, the pencil ( Ao-lambda*Eo ) has full column
C                              (      Co      )
C   rank o for all finite lambda and is in a staircase form with
C                    _      _            _      _
C                  ( Ak,k   Ak,k-1   ... Ak,2   Ak,1   )
C                  ( _      _            _      _      )
C      ( Ao )  =   ( Ak-1,k Ak-1,k-1 ... Ak-1,2 Ak-1,1 ) ,           (9)
C      ( Co )      (     ...         ... _      _      )
C                  (  0       0      ... A1,2   A1,1   )
C                  (                            _      )
C                  (  0       0      ... 0      A0,1   )
C                    _          _      _
C                  ( Ek,k  ...  Ek,2   Ek,1 )
C                  (       ...  _      _    )
C        Eo    =   (   0   ...  E2,2   E2,1 ) ,                     (10)
C                  (                   _    )
C                  (   0   ...    0    E1,1 )
C         _
C   where Ai-1,i is a ctau(i-1)-by-ctau(i) full column rank matrix
C                          _
C   (with ctau(0) = p) and Ei,i is a ctau(i)-by-ctau(i) upper
C   triangular matrix.
C
C   For jobobs = 0, the (n-o)-by-(n-o) regular pencil Ano - lambda*Eno
C   has the form
C
C                         ( Afno - lambda*Efno         *          )
C      Ano - lambda*Eno = (                                       ) ,
C                         (        0           Aino - lambda*Eino )
C
C   where:
C     1) the ino-by-ino regular pencil Aino - lambda*Eino, with Aino
C        upper triangular and nonsingular, contains the unobservable
C        infinite eigenvalues of A - lambda*E;
C     2) the (n-o-ino)-by-(n-o-ino) regular pencil Afno - lambda*Efno,
C        with Efno upper triangular and nonsingular, contains the
C        unobservable finite eigenvalues of A - lambda*E.
C
C   Note: The significance of the two diagonal blocks can be
C         interchanged by calling the gateway with the
C         arguments A and E interchanged. In this case,
C         Aino - lambda*Eino contains the unobservable zero
C         eigenvalues of A - lambda*E, while Afno - lambda*Efno
C         contains the unobservable nonzero finite and infinite
C         eigenvalues of A - lambda*E.
C
C   For jobobs = 1, the pencil Ano - lambda*Eno has the form
C
C      Ano - lambda*Eno = Afno - lambda*Efno ,
C
C   where the regular pencil Afno - lambda*Efno, with Efno upper
C   triangular and nonsingular, contains the unobservable finite
C   eigenvalues of A - lambda*E.
C
C   For jobobs = 2, the pencil Ano - lambda*Eno has the form
C
C      Ano - lambda*Eno = Aino - lambda*Eino ,
C
C   where the regular pencil Aino - lambda*Eino, with Aino upper
C   triangular and nonsingular, contains the unobservable nonzero
C   finite and infinite eigenvalues of A - lambda*E.
C
C   3) To find a reduced (controllable, observable, or irreducible)
C   descriptor representation (Ar-lambda*Er,Br,Cr) for an original
C   descriptor representation (A-lambda*E,B,C). The pencil Ar-lambda*Er
C   is in an upper block Hessenberg form, with either Ar or Er upper
C   triangular.
C
C Input parameters:
C   task   - integer specifying the computations to be performed.
C            task = 1 :  compute controllability staircase form (1);
C            task = 2 :  compute observability staircase form (6);
C            task = 3 :  compute controllable, observable, or
C                        irreducible form.
C   A      - the n-by-n state dynamics matrix A.
C   E      - the n-by-n descriptor matrix E.
C   B      - the n-by-m input/state matrix B.
C   C      - the p-by-n state/output matrix C.
C   flag   - (optional) real vector of length 4 specifying various
C            options, depending on task.
C
C            For task = 1, flag contains:
C            flag(1) = jobcon : indicates what to do, as follows:
C               jobcon = 0 :  separate both finite and infinite
C                             uncontrollable eigenvalues;
C               jobcon = 1 :  separate only finite uncontrollable
C                             eigenvalues;
C               jobcon = 2 :  separate only nonzero finite and infinite
C                             uncontrollable eigenvalues.
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
C            flag(4) = tol : see below.
C            Default :  flag = [ 0; 0; 0; 0 ].
C
C            For task = 2, flag contains:
C            flag(1) = jobobs : indicates what to do, as follows:
C               jobobs = 0 :  separate both finite and infinite
C                             unobservable eigenvalues;
C               jobobs = 1 :  separate only finite unobservable
C                             eigenvalues;
C               jobobs = 2 :  separate only nonzero finite and infinite
C                             unobservable eigenvalues.
C            flag(2) = compq : see task = 1 above.
C            flag(3) = compz : see above.
C            flag(4) = tol : see below.
C            Default :  flag = [ 0; 0; 0; 0 ].
C
C            For task = 3, flag contains:
C            flag(1) = job : indicates what to do, as follows:
C               job = 0 :  remove both the uncontrollable and
C                          unobservable parts to get an irreducible
C                          descriptor representation;
C               job = 1 :  remove the uncontrollable part only to get a
C                          controllable descriptor representation;
C               job = 2 :  remove the unobservable part only to get an
C                          observable descriptor representation.
C            flag(2) = systyp : indicates the type of descriptor system
C            algorithm to be applied according to the assumed
C            transfer-function matrix as follows:
C               systyp = 0 :  rational transfer-function matrix;
C               systyp = 1 :  proper (standard) transfer-function
C                             matrix;
C               systyp = 2 :  polynomial transfer-function matrix.
C            flag(3) = equil : specifies whether the user wishes to
C            preliminarily scale the system (A-lambda*E,B,C) as follows:
C               equil = 0 :  perform scaling;
C               equil = 1 :  do not perform scaling.
C            flag(4) = tol : see below.
C            Default :  flag = [ 0; 0; 0; 0 ].
C
C            tol is a real scalar indicating the tolerance to be used
C            in rank determinations when transforming (A-lambda*E,B,C),
C            or (A-lambda*E, B) if task = 1, or (A'-lambda*E',C')', if
C            task = 2. If tol > 0, then the given value of tol is used
C            as a lower bound for reciprocal condition numbers in rank
C            determinations; a (sub)matrix whose estimated condition
C            number is less than 1/tol is considered to be of full rank.
C            If tol <= 0, the default tolerance toldef = eps*n*n is
C            used instead, where eps is the machine precision. tol < 1.
C   Q1     - (optional) if compq = 2 the n-by-n given orthogonal
C            matrix Q1.
C   Z1     - (optional) if compz = 2 the n-by-n given orthogonal
C            matrix Z1.
C
C Output parameters:
C   Ao     - If task = 1, the n-by-n transformed state matrix Q'*A*Z,
C
C                         ( Ac   *  )
C                Q'*A*Z = (         ) ,
C                         ( 0   Anc )
C
C            where Ac is c-by-c and Anc is (n-c)-by-(n-c).
C            If jobcon = 1, the matrix ( Bc Ac ) is in the
C            controllability staircase form (4).
C            If jobcon = 0 or 2, the submatrix Ac is upper triangular.
C            If jobcon = 0, the Anc matrix has the form
C
C                      ( Ainc   *  )
C                Anc = (           ) ,
C                      (  0   Afnc )
C
C            where the inc-by-inc matrix Ainc is nonsingular and upper
C            triangular.
C            If jobcon = 2, Anc is nonsingular and upper triangular.
C
C            If task = 2, the n-by-n transformed state matrix Q'*A*Z,
C
C                         ( Ano  *  )
C                Q'*A*Z = (         ) ,
C                         ( 0    Ao )
C
C            where Ao is o-by-o and Ano is (n-o)-by-(n-o).
C            If jobobs = 1, the matrix ( Ao ) is in the observability
C                                      ( Co )
C            staircase form (9).
C            If jobobs = 0 or 2, the submatrix Ao is upper triangular.
C            If jobobs = 0, the submatrix Ano has the form
C
C                      ( Afno   *  )
C                Ano = (           ) ,
C                      (  0   Aino )
C
C            where the ino-by-ino matrix Aino is nonsingular and upper
C            triangular.
C            If jobobs = 2, Ano is nonsingular and upper triangular.
C
C            If task = 3, the nr-by-nr reduced order state matrix Ar of
C            an irreducible, controllable, or observable realization for
C            the original system, depending on the value of job,
C            job = 0, job = 1, or job = 2, respectively.
C            The matrix Ar is upper triangular if systyp = 0 or 2.
C            If systyp = 1 and job = 1, the matrix [Br Ar] is in a
C            controllable staircase form.
C            If systyp = 1 and job = 0 or 2, the matrix ( Ar ) is in an
C                                                       ( Cr )
C            observable staircase form.
C            The block structure of staircase forms is contained
C            in the leading infred(7) elements of the vector sizes.
C   Eo     - If task = 1, the n-by-n transformed descriptor matrix
C            Q'*E*Z,
C
C                        ( Ec   *  )
C               Q'*E*Z = (         ) ,
C                        ( 0   Enc )
C
C            where Ec is c-by-c and Enc is (n-c)-by-(n-c).
C            If jobcon = 0 or 2, the matrix ( Bc Ec ) is in the
C            controllability staircase form (2).
C            If jobcon = 1, the submatrix Ec is upper triangular.
C            If jobcon = 0, the Enc matrix has the form
C
C                      ( Einc   *  )
C                Enc = (           ) ,
C                      (  0   Efnc )
C
C            where the inc-by-inc matrix Einc is nilpotent and the
C            (n-c-inc)-by-(n-c-inc) matrix Efnc is nonsingular and
C            upper triangular.
C            If jobcon = 1, Enc is nonsingular and upper triangular.
C
C            If task = 2, the n-by-n transformed descriptor matrix
C            Q'*E*Z,
C
C                         ( Eno  *  )
C                Q'*E*Z = (         ) ,
C                         ( 0    Eo )
C
C            where Eo is o-by-o and Eno is (n-o)-by-(n-o).
C            If jobobs = 0 or 2, the matrix ( Eo ) is in the
C                                           ( Co )
C            observability staircase form (7).
C            If jobobs = 1, the submatrix Eo is upper triangular.
C            If jobobs = 0, the Eno matrix has the form
C
C                      ( Efno   *  )
C                Eno = (           ) ,
C                      (  0   Eino )
C
C            where the ino-by-ino matrix Eino is nilpotent and the
C            (n-o-ino)-by-(n-o-ino) matrix Efno is nonsingular and
C            upper triangular.
C            If jobobs = 1, Eno is nonsingular and upper triangular.
C
C            If task = 3, the nr-by-nr reduced order descriptor matrix
C            Er of an irreducible, controllable, or observable
C            realization for the original system, depending on the value
C            of job, job = 0, job = 1, or job = 2, respectively.
C            The resulting Er has infred(6) nonzero sub-diagonals.
C            If at least for one k = 1,...,4, infred(k) >= 0, then the
C            resulting Er is structured being either upper triangular
C            or block Hessenberg, in accordance to the last
C            performed order reduction phase (see Method).
C            The block structure of staircase forms is contained
C            in the leading infred(7) elements of the vector sizes.
C   Bo     - If task = 1 or task = 2, the n-by-m transformed state/input
C            matrix Q'*B. If task = 1, this matrix has the form
C
C                       ( Bc )
C                Q'*B = (    ) ,
C                       ( 0  )
C
C            where Bc is c-by-m.
C            For jobcon = 0 or 2, the matrix ( Bc Ec ) is in the
C            controllability staircase form (2).
C            For jobcon = 1, the matrix ( Bc Ac ) is in the
C            controllability staircase form (4).
C            If task = 3, the nr-by-m reduced input matrix Br of an
C            irreducible, controllable, or observable realization for
C            the original system, depending on the value of job,
C            job = 0, job = 1, or job = 2, respectively.
C            If job = 1, only the first sizes(1) rows of B are nonzero.
C   Co     - If task = 1, the p-by-n transformed matrix C*Z.
C            If task = 2, the p-by-n transformed matrix
C
C                C*Z = (  0   Co ) ,
C
C            where Co is p-by-o.
C            If jobobs = 0 or 2, the matrix ( Eo ) is in the
C                                           ( Co )
C            observability staircase form (7).
C            If jobobs = 1, the matrix ( Ao ) is in the observability
C                                      ( Co )
C            staircase form (9).
C            If task = 3, the p-by-nr transformed state/output matrix Cr
C            of an irreducible, controllable, or observable realization
C            for the original system, depending on the value of job,
C            job = 0, job = 1, or job = 2, respectively.
C            If job = 0, or job = 2, only the last sizes(1) columns
C            (in the first nr columns) of C are nonzero.
C   Q      - the n-by-n orthogonal matrix Q.
C            If compq = 1, Q' is the product of the transformations
C            which are applied to A, E, and B on the left.
C            If compq = 2, Q is the orthogonal matrix product Q1*Q.
C   Z      - the n-by-n orthogonal matrix Z.
C            If compz = 1, Z is the product of the transformations
C            which are applied to A, E, and C on the right.
C            If compz = 2, Z is the orthogonal matrix product Z1*Z.
C   orders - (optional) if task < 3, integer vector of length 3
C            indicating the orders of the subsystems and the number of
C            full rank blocks.
C            If task = 1, orders contains c, inc, and nrb:
C               c is the order of the reduced matrices Ac and Ec, and
C            the number of rows of reduced matrix Bc; also the order
C            of the controllable part of the pair (A-lambda*E,B).
C               For jobcon = 0, inc is the order of the reduced matrices
C            Ainc and Einc, and also the number of uncontrollable
C            infinite eigenvalues of the pencil A - lambda*E.
C            For jobcon <> 0, inc has no significance and it is set to
C            zero.
C                                                            _
C               nrb is the number k, of full row rank blocks Ei,i in the
C            staircase form of the pencil (Bc Ec-lambda*Ac) in (2) and
C                                                                  _
C            (3) (for jobcon = 0 or 2), or of full row rank blocks Ai,i
C            in the staircase form of the pencil (Bc Ac-lambda*Ec) in
C            (4) and (5) (for jobcon = 1).
C
C            If task = 2, orders contains o, ino, and nrb:
C               o is the order of the reduced matrices Ao and Eo, and
C            the number of columns of reduced matrix Co; also the order
C            of the observable part of the pair (C, A-lambda*E).
C               For jobobs = 0, ino is the order of the reduced matrices
C            Aino and Eino, and also the number of unobservable
C            infinite eigenvalues of the pencil A - lambda*E.
C            For jobobs <> 0, ino has no significance and it is set to
C            zero.
C                                                               _
C               nrb is the number k, of full column rank blocks Ei-1,i
C            in the staircase form of the pencil (Eo-lambda*Ao) in (7)
C                                                (    Co      )
C            and (8) (for jobobs = 0 or 2), or of full column rank
C                   _
C            blocks Ai-1,i in the staircase form of the pencil
C            (Ao-lambda*Eo) in (9) and (10) (for jobobs = 1).
C            (     Co     )
C   infred - (optional) if task = 3, integer array of dimension 7,
C            containing information on performed reduction and on
C            structure of resulting system matrices as follows:
C            infred(k) >= 0 (k = 1, 2, 3, or 4) if Phase k of reduction
C                           (see Method) has been performed. In this
C                           case, infred(k) is the achieved order
C                           reduction in Phase k.
C            infred(k) < 0  (k = 1, 2, 3, or 4) if Phase k was not
C                           performed.
C            infred(5)  -   the number of nonzero sub-diagonals of A.
C            infred(6)  -   the number of nonzero sub-diagonals of E.
C            infred(7)  -   the number of blocks in the resulting
C                           staircase form at last performed reduction
C                           phase. The block dimensions are contained
C                           in the first infred(7) elements of sizes.
C   sizes -  (optional) if task = 1 or 2, integer nrb-vector containing
C            rtau, if task = 1, or ctau, if task = 2.
C               rtau(i), for i = 1, ..., nrb, is the row dimension of
C                                    _         _
C            the full row rank block Ei,i-1 or Ai,i-1 in the staircase
C            form (2) or (4) for jobcon = 0 or 2, or for jobcon = 1,
C            respectively.
C               ctau(i), for i = 1, ..., nrb, is the column dimension
C                                          _         _
C            of the full column rank block Ei-1,i or Ai-1,i in the
C            staircase form (7) or (9) for jobobs = 0 or 2, or
C            for jobobs = 1, respectively.
C            If task = 3, integer vector of dimension infred(7), whose
C            elements contain the orders of the diagonal blocks of
C            Ar-lambda*Er.
C
C Method
C   If task = 3, the order reduction is performed in 4 phases:
C   Phase 1: Eliminate all finite uncontrollable eigenvalues.
C            The resulting matrix ( Br Ar ) is in a controllable
C            staircase form, and Er is upper triangular.
C            This phase is performed if job = 0 or 1 and systyp = 0
C            or 1.
C   Phase 2: Eliminate all infinite and finite nonzero uncontrollable
C            eigenvalues. The resulting matrix ( Br Er ) is in a
C            controllable staircase form, and Ar is upper triangular.
C            This phase is performed if job = 0 or 1 and systyp = 0
C            or 2.
C   Phase 3: Eliminate all finite unobservable eigenvalues.
C            The resulting matrix ( Ar ) is in an observable
C                                 ( Cr )
C            staircase form, and Er is upper triangular.
C            This phase is performed if job = 0 or 2 and systyp = 0
C            or 1.
C   Phase 4: Eliminate all infinite and finite nonzero unobservable
C            eigenvalues. The resulting matrix ( Er ) is in an
C                                              ( Cr )
C            observable staircase form, and Ar is upper triangular.
C            This phase is performed if job = 0 or 2 and systyp = 0
C            or 2.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003.
C
C Revisions:
C   V. Sima, Feb. 2004, Apr. 2009, Dec. 2010, July 2012, Apr. 2017.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
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
      CHARACTER         COMPQ, COMPZ, EQUIL, JOB, JOBCON, JOBOBS,
     $                  SYSTYP
      INTEGER           INFO, LDA, LDB, LDC, LDE, LDQ, LDWORK, LDZ, M,
     $                  N, NCONT, NIUCON, NIUOBS, NLBLCK, NOBSV, NR,
     $                  NRBLCK, P
      DOUBLE PRECISION  TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: CTAU(:), INFRED(:), IWORK(:),
     $                                 RTAU(:)
      DOUBLE PRECISION, ALLOCATABLE :: DW(:), TMP(:)
      mwPointer         A, B, C, DWK, DWORK, E, Q, Z
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           WITHQ, WITHQ1, WITHZ, WITHZ1
      INTEGER           FLAG(4), I, ICOMPQ, ICOMPZ, IEQU, IJOB, IP, IPS,
     $                  ISIZE, ISYS, ITASK, LIWORK, MB, NM, NQ, NZ
      DOUBLE PRECISION  FLAGR(4), TEMP
C
C .. External subroutines ..
      EXTERNAL          DLACPY, TG01HD, TG01ID, TG01JD
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'GSYSCOM requires at least 5 input arguments' )
      ELSE IF ( NLHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'GSYSCOM requires at least 4 output arguments' )
      END IF
C
C   task, A(nxn), E(nxn), B(n,m), C(p,n) (, flag, Q1(n,n), Z1(n,n))
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
      IF ( ITASK.LT.1 .OR. ITASK.GT.3 ) THEN
         CALL mexErrMsgTxt
     $     ( 'TASK has 1, 2, or 3 the only admissible values' )
      END IF
C
      N =  mxGetM( PRHS(2) )
      IF ( mxGetN( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
C
      IF ( mxGetM( PRHS(3) ).NE.N .OR.
     $     mxGetN( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'E must have the same size as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'E must be a real matrix' )
      END IF
C
      IP = 6
      M  = mxGetN( PRHS(4) )
      IF ( mxGetM( PRHS(4) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same rowsize as A' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
C
      P =  mxGetM( PRHS(5) )
      IF ( mxGetN( PRHS(5) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'C must have the same columnsize as A' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
C
C   flag
C
      DO 10 I = 1, 4
         FLAG(I)  = 0
         FLAGR(I) = ZERO
   10 CONTINUE
C
      IF ( NRHS.GE.IP ) THEN
         ISIZE = mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) )
         IF ( ISIZE.GT.4 ) THEN
            CALL mexErrMsgTxt( 'FLAG must have at most 4 elements' )
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
      IPS = IP
C
      TOL = FLAGR(4)
      IF ( TOL.GT.ONE ) THEN
         CALL mexErrMsgTxt( 'TOL must be less than 1' )
      END IF
      ICOMPQ = 0
      ICOMPZ = 0
C
      IJOB = FLAG(1)
      IF ( ITASK.EQ.1 ) THEN
         IF ( IJOB.EQ.0 ) THEN
            JOBCON = 'C'
         ELSE IF ( IJOB.EQ.1 ) THEN
            JOBCON = 'F'
         ELSE IF ( IJOB.EQ.2 ) THEN
            JOBCON = 'I'
         ELSE
            CALL mexErrMsgTxt
     $        ( 'JOBCON has 0, 1, or 2 the only admissible values' )
         END IF
         ICOMPQ = FLAG(2)
         ICOMPZ = FLAG(3)
C
      ELSE IF ( ITASK.EQ.2 ) THEN
         IF ( IJOB.EQ.0 ) THEN
            JOBOBS = 'O'
         ELSE IF ( IJOB.EQ.1 ) THEN
            JOBOBS = 'F'
         ELSE IF ( IJOB.EQ.2 ) THEN
            JOBOBS = 'I'
         ELSE
            CALL mexErrMsgTxt
     $        ( 'JOBOBS has 0, 1, or 2 the only admissible values' )
         END IF
         ICOMPQ = FLAG(2)
         ICOMPZ = FLAG(3)
C
      ELSE IF ( ITASK.EQ.3 ) THEN
         IF ( IJOB.EQ.0 ) THEN
            JOB = 'I'
         ELSE IF ( IJOB.EQ.1 ) THEN
            JOB = 'C'
         ELSE IF ( IJOB.EQ.2 ) THEN
            JOB = 'O'
         ELSE
            CALL mexErrMsgTxt
     $        ( 'JOB has 0, 1, or 2 the only admissible values' )
         END IF
C
         ISYS = FLAG(2)
         IF ( ISYS.EQ.0 ) THEN
            SYSTYP = 'R'
         ELSE IF ( ISYS.EQ.1 ) THEN
            SYSTYP = 'S'
         ELSE IF ( ISYS.EQ.2 ) THEN
            SYSTYP = 'P'
         ELSE
            CALL mexErrMsgTxt
     $        ( 'SYSTYP has 0, 1, or 2 the only admissible values' )
         END IF
C
         IEQU = FLAG(3)
         IF ( IEQU.EQ.0 ) THEN
            EQUIL = 'S'
         ELSE IF ( IEQU.EQ.1 ) THEN
            EQUIL = 'N'
         ELSE
            CALL mexErrMsgTxt
     $        ( 'EQUIL has 0, or 1 the only admissible values' )
         END IF
      END IF
C
      WITHQ1 = ITASK.EQ.1 .OR. ITASK.EQ.2
      WITHZ1 = WITHQ1
C
      IF ( WITHQ1 ) THEN
         IF ( ICOMPQ.EQ.0 ) THEN
            COMPQ = 'N'
         ELSE IF ( ICOMPQ.EQ.1 ) THEN
            COMPQ = 'I'
         ELSE IF ( ICOMPQ.EQ.2 ) THEN
            COMPQ = 'U'
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
            COMPZ = 'U'
         ELSE
            CALL mexErrMsgTxt
     $        ( 'COMPZ has 0, 1, or 2 the only admissible values' )
         END IF
      END IF
C
      WITHQ  = ICOMPQ.GE.1 .AND. NLHS.GT.4
      WITHZ  = ICOMPZ.GE.1 .AND. NLHS.GT.4 + MIN( ICOMPQ, 1 )
C
      WITHQ1 = ICOMPQ.EQ.2
      WITHZ1 = ICOMPZ.EQ.2
C
C   Q1
C
      IF ( NRHS.GE.IP .AND. WITHQ1 ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.N .OR.
     $        mxGetN( PRHS(IP) ).NE.N ) THEN
            CALL mexErrMsgTxt( 'Q1 must have the same size as A' )
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
            CALL mexErrMsgTxt( 'Z1 must have the same size as A' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'Z1 must be a real matrix' )
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA = MAX( 1, N )
      LDE = LDA
      LDB = LDA
      LDC = MAX( 1, P )
      MB  = M
      IF ( ITASK.EQ.2 ) THEN
         MB  = MAX( M, P )
         LDC = MAX( LDC, M )
      ELSE IF ( ITASK.EQ.3 ) THEN
         IF ( N.GT.0 ) THEN
            LDC = MAX( LDC, M )
         ELSE
            LDC = 1
         END IF
         IF ( IJOB.NE.1 )
     $      MB = MAX( M, P )
      END IF
      IF ( WITHQ1 .OR. WITHQ ) THEN
         LDQ = LDA
         NQ  = N
      ELSE
         LDQ = 1
         NQ  = 1
      END IF
      IF ( WITHZ1 .OR. WITHZ ) THEN
         LDZ = LDA
         NZ  = N
      ELSE
         LDZ = 1
         NZ  = 1
      END IF
C
C   ldwork
C   Use LDWORK >= 2*N*N+N*M+N*P+MAX(N,2*M,2*P) for TG01JD, to get
C   more accurate results by performing only those reductions phases
C   (see Method), where effective order reduction occurs.
C
      IF ( ITASK.EQ.1 ) THEN
         LIWORK = M
         LDWORK = MAX( N, 2*M )
      ELSE IF ( ITASK.EQ.2 ) THEN
         LIWORK = P
         LDWORK = MAX( N, 2*P )
      ELSE IF ( ITASK.EQ.3 ) THEN
         LIWORK = N + MAX( M, P )
         IF ( IJOB.EQ.0 .OR. ISYS.EQ.0 )
     $      LIWORK = LIWORK + N
         LDWORK = 2*N*N + N*M + N*P +  MAX( N, 2*M, 2*P )
         IF ( IEQU.EQ.0 )
     $      LDWORK = MAX( LDWORK, 8*N )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( IWORK( LIWORK ) )
      A     = mxCalloc(  LDA*N,  8 )
      B     = mxCalloc(  LDB*MB, 8 )
      C     = mxCalloc(  LDC*N,  8 )
      DWORK = mxCalloc( LDWORK,  8 )
      E     = mxCalloc(  LDE*N,  8 )
      IF ( WITHQ .OR. WITHQ1 )
     $   Q = mxCalloc( LDQ*NQ, 8 )
      IF ( WITHZ .OR. WITHZ1 )
     $   Z = mxCalloc( LDZ*NZ, 8 )
      IF ( ITASK.EQ.1 ) THEN
         ALLOCATE ( RTAU( N ) )
      ELSE IF ( ITASK.EQ.2 ) THEN
         ALLOCATE ( CTAU( N ) )
      ELSE
         ALLOCATE ( INFRED( 7 ), TMP( 7 ) )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( A ), N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL( E ), N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL( B ), N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), %VAL( C ), P*N )
      IP = IPS
      IF ( NRHS.GE.IP .AND. WITHQ1 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL( Q ), N*N )
         IP = IP + 1
      END IF
      IF ( NRHS.GE.IP .AND. WITHZ1 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), %VAL( Z ), N*N )
      END IF
C
C Do the actual computations.
C
      IF ( ITASK.GE.2 ) THEN
         IF ( M.GT.P .AND. P.GT.0 ) THEN
            DWK = mxCalloc( P*N, 8 )
            CALL mxCopyPtrToReal8( C, %VAL( DWK ), P*N )
            CALL DLACPY( 'Full', P, N, %VAL( DWK ), P, %VAL( C ), LDC )
            CALL mxFree( DWK )
         END IF
      END IF
C
      IF ( ITASK.EQ.1 ) THEN
         CALL TG01HD( JOBCON, COMPQ, COMPZ, N, M, P, %VAL( A ), LDA,
     $                %VAL( E ), LDE, %VAL( B ), LDB, %VAL( C ), LDC,
     $                %VAL( Q ), LDQ, %VAL( Z ), LDZ, NCONT, NIUCON,
     $                NRBLCK, RTAU, TOL, IWORK, %VAL( DWORK ), INFO )
      ELSE IF ( ITASK.EQ.2 ) THEN
         CALL TG01ID( JOBOBS, COMPQ, COMPZ, N, M, P, %VAL( A ), LDA,
     $                %VAL( E ), LDE, %VAL( B ), LDB, %VAL( C ), LDC,
     $                %VAL( Q ), LDQ, %VAL( Z ), LDZ, NOBSV, NIUOBS,
     $                NLBLCK, CTAU, TOL, IWORK, %VAL( DWORK ), INFO )
      ELSE
         CALL TG01JD( JOB, SYSTYP, EQUIL, N, M, P, %VAL( A ), LDA,
     $                %VAL( E ), LDE, %VAL( B ), LDB, %VAL( C ), LDC,
     $                NR, INFRED, TOL, IWORK, %VAL( DWORK ), LDWORK,
     $                INFO )
      END IF
C
      IF ( ITASK.GE.2 .AND. M.GT.P .AND. P.GT.0 )
     $   CALL DLACPY( 'Full', P, N, %VAL( C ), LDC, %VAL( C ), P )
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         IP = 5
         IF ( ITASK.EQ.3 ) THEN
            IF ( N.GT.NR .AND. NR.GT.0 ) THEN
               CALL DLACPY( 'Full', N, N, %VAL( A ), LDA, %VAL( A ),
     $                      NR )
               CALL DLACPY( 'Full', N, N, %VAL( E ), LDE, %VAL( E ),
     $                      NR )
               CALL DLACPY( 'Full', N, M, %VAL( B ), LDB, %VAL( B ),
     $                      NR )
            END IF
            IF ( M.GT.P .AND. P.GT.0 .AND. N.GT.0 )
     $         CALL DLACPY( 'Full', P, N, %VAL( C ), LDC, %VAL( C ), P )
            PLHS(1) = mxCreateDoubleMatrix( NR, NR, 0 )
            CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS(1) ),
     $                             NR*NR )
            PLHS(2) = mxCreateDoubleMatrix( NR, NR, 0 )
            CALL mxCopyReal8ToPtr( %VAL( E ), mxGetPr( PLHS(2) ),
     $                             NR*NR )
            PLHS(3) = mxCreateDoubleMatrix( NR, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( B ), mxGetPr( PLHS(3) ),
     $                             NR*M )
            PLHS(4) = mxCreateDoubleMatrix( P, NR, 0 )
            CALL mxCopyReal8ToPtr( %VAL( C ), mxGetPr( PLHS(4) ),
     $                             P*NR )
            IF ( NLHS.GE.IP ) THEN
C
               DO 20 I = 1, 7
                  TMP(I) = INFRED(I)
   20          CONTINUE
C
               PLHS(IP) = mxCreateDoubleMatrix( 7, 1, 0 )
               CALL mxCopyReal8ToPtr( TMP, mxGetPr( PLHS(IP) ), 7 )
               IP = IP + 1
            END IF
            IF ( NLHS.GE.IP ) THEN
               NM = INFRED(7)
               ALLOCATE( DW( NM ) )
C
               DO 30 I = 1, NM
                  DW(I) = IWORK(I)
   30          CONTINUE
C
               PLHS(IP) = mxCreateDoubleMatrix( NM, MIN( NM, 1 ), 0 )
               CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(IP) ),
     $                                NM*MIN( NM, 1 ) )
               DEALLOCATE( DW )
            END IF
         ELSE
            PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS(1) ),
     $                             N*N )
            PLHS(2) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( E ), mxGetPr( PLHS(2) ),
     $                             N*N )
            PLHS(3) = mxCreateDoubleMatrix( N, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( B ), mxGetPr( PLHS(3) ),
     $                             N*M )
            PLHS(4) = mxCreateDoubleMatrix( P, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( C ), mxGetPr( PLHS(4) ),
     $                             P*N )
            IF ( NLHS.GE.IP .AND. WITHQ ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
               CALL mxCopyReal8ToPtr( %VAL( Q ),
     $                                mxGetPr( PLHS(IP) ), N*N )
               IP = IP + 1
            END IF
            IF ( NLHS.GE.IP .AND. WITHZ ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
               CALL mxCopyReal8ToPtr( %VAL( Z ),
     $                                mxGetPr( PLHS(IP) ), N*N )
               IP = IP + 1
            END IF
            IF ( NLHS.GE.IP ) THEN
               IF ( ITASK.EQ.1 ) THEN
                  FLAGR(1) = NCONT
                  FLAGR(2) = NIUCON
                  FLAGR(3) = NRBLCK
                  NM       = NRBLCK
               ELSE
                  FLAGR(1) = NOBSV
                  FLAGR(2) = NIUOBS
                  FLAGR(3) = NLBLCK
                  NM       = NLBLCK
               END IF
               PLHS(IP) = mxCreateDoubleMatrix( 3, 1, 0 )
               CALL mxCopyReal8ToPtr( FLAGR, mxGetPr( PLHS(IP) ), 3 )
               IP = IP + 1
            END IF
            IF ( NLHS.GE.IP ) THEN
               ALLOCATE( DW( NM ) )
               IF ( ITASK.EQ.1 ) THEN
C
                  DO 40 I = 1, NM
                     DW(I) = RTAU(I)
   40             CONTINUE
C
               ELSE
C
                  DO 50 I = 1, NM
                     DW(I) = CTAU(I)
   50             CONTINUE
C
               END IF
               PLHS(IP) = mxCreateDoubleMatrix( NM, MIN( NM, 1 ), 0 )
               CALL mxCopyReal8ToPtr( DW, mxGetPr( PLHS(IP) ),
     $                                NM*MIN( NM, 1 ) )
               DEALLOCATE( DW )
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( IWORK )
      CALL mxFree( A     )
      CALL mxFree( B     )
      CALL mxFree( C     )
      CALL mxFree( DWORK )
      CALL mxFree( E     )
      IF ( WITHQ .OR. WITHQ1 )
     $   CALL mxFree( Q )
      IF ( WITHZ .OR. WITHZ1 )
     $   CALL mxFree( Z )
      IF ( ITASK.EQ.1 ) THEN
         DEALLOCATE ( RTAU )
      ELSE IF ( ITASK.EQ.2 ) THEN
         DEALLOCATE ( CTAU )
      ELSE
         DEALLOCATE ( INFRED, TMP )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( ITASK.EQ.1 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM TG01HD'')' )
     $         INFO
         ELSE IF ( ITASK.EQ.2 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM TG01ID'')' )
     $         INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM TG01JD'')' )
     $         INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of GSYSCOM ***
      END
