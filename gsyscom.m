% GSYSCOM.F  - MEX-function to transform a descriptor system, by
%              equivalence transformations, to a controllable or
%              observable staircase form, or to a reduced (controllable,
%              observable, or irreducible) form, using SLICOT routines
%              TG01HD, TG01ID, and TG01JD.
%
%   [Ao,Eo,Bo,Co(,Q,Z)(,orders,sizes)] =
%                       GSYSCOM(task,A,E,B,C(,flag)(,Q1,Z1))
%
%   [Ao,Eo,Bo,Co(,Q,Z)(,orders,sizes)] =
%                       GSYSCOM(1,A,E,B,C(,flag)(,Q1,Z1))
%   [Ao,Eo,Bo,Co(,Q,Z)(,orders,sizes)] =
%                       GSYSCOM(2,A,E,B,C(,flag)(,Q1,Z1))
%   [Ao,Eo,Bo,Co(,infred,sizes)] =
%                       GSYSCOM(3,A,E,B,C(,flag))
%
%   GSYSCOM performs one of the equivalence transformations specified by
%   the value of the parameter task (task = 1, 2, or 3), for a
%   descriptor triple (A-lambda E,B,C), to reduce it to a staircase form,
%   or to a controllable, observable, or irreducible form:
%
%   1) To compute orthogonal transformation matrices Q and Z which
%   reduce the n-th order descriptor system (A-lambda*E,B,C) to the form
%
%               ( Ac  *  )             ( Ec  *  )           ( Bc )
%      Q'*A*Z = (        ) ,  Q'*E*Z = (        ) ,  Q'*B = (    ) ,
%               ( 0  Anc )             ( 0  Enc )           ( 0  )
%
%         C*Z = ( Cc Cnc ) ,                                         (1)
%
%   where the c-th order descriptor system (Ac-lambda*Ec,Bc,Cc) is
%   finite and/or infinite controllable. The pencil Anc - lambda*Enc
%   is regular of order n-c and contains the uncontrollable finite
%   and/or infinite eigenvalues of the pencil A-lambda*E. The reduced
%   order descriptor system (Ac-lambda*Ec,Bc,Cc) has the same
%   transfer-function matrix as the original system (A-lambda*E,B,C).
%   The left and/or right orthogonal transformations Q and Z performed
%   to reduce the system matrices can be optionally accumulated.
%
%   For jobcon = 0 or 2 (see parameter flag below), the pencil
%   ( Bc Ec-lambda*Ac ) has full row rank c for all finite lambda and
%   is in a staircase form with
%                    _      _          _        _
%                  ( E1,0   E1,1  ...  E1,k-1   E1,k  )
%                  (        _          _        _     )
%      ( Bc Ec ) = (  0     E2,1  ...  E2,k-1   E2,k  ) ,            (2)
%                  (              ...  _        _     )
%                  (  0       0   ...  Ek,k-1   Ek,k  )
%
%                    _          _        _
%                  ( A1,1  ...  A1,k-1   A1,k  )
%                  (            _        _     )
%        Ac      = (   0   ...  A2,k-1   A2,k  ) ,                   (3)
%                  (       ...           _     )
%                  (   0   ...    0      Ak,k  )
%         _
%   where Ei,i-1 is an rtau(i)-by-rtau(i-1) full row rank matrix
%                          _
%   (with rtau(0) = m) and Ai,i is an rtau(i)-by-rtau(i) upper
%   triangular matrix.
%
%   For jobcon = 1, the pencil ( Bc Ac-lambda*Ec ) has full row
%   rank c for all finite lambda and is in a staircase form with
%                    _     _          _        _
%                  ( A1,0  A1,1  ...  A1,k-1   A1,k  )
%                  (       _          _        _     )
%      ( Bc Ac ) = (  0    A2,1  ...  A2,k-1   A2,k  ) ,             (4)
%                  (             ...  _        _     )
%                  (  0      0   ...  Ak,k-1   Ak,k  )
%
%                    _          _        _
%                  ( E1,1  ...  E1,k-1   E1,k  )
%                  (            _        _     )
%        Ec      = (   0   ...  E2,k-1   E2,k  ) ,                   (5)
%                  (       ...           _     )
%                  (   0   ...    0      Ek,k  )
%         _
%   where Ai,i-1 is an rtau(i)-by-rtau(i-1) full row rank matrix
%                          _
%   (with rtau(0) = m) and Ei,i is an rtau(i)-by-rtau(i) upper
%   triangular matrix.
%
%   For jobcon = 0, the (n-c)-by-(n-c) regular pencil Anc - lambda*Enc
%   has the form
%
%                         ( Ainc - lambda*Einc         *          )
%      Anc - lambda*Enc = (                                       ) ,
%                         (        0           Afnc - lambda*Efnc )
%
%   where:
%     1) the inc-by-inc regular pencil Ainc - lambda*Einc, with Ainc
%        upper triangular and nonsingular, contains the uncontrollable
%        infinite eigenvalues of A - lambda*E;
%     2) the (n-c-inc)-by-(n-c-inc) regular pencil Afnc - lambda*Efnc,
%        with Efnc upper triangular and nonsingular, contains the
%        uncontrollable finite eigenvalues of A - lambda*E.
%
%   Note: The significance of the two diagonal blocks can be
%         interchanged by calling the gateway with the arguments A and E
%         interchanged. In this case, Ainc - lambda*Einc contains the
%         uncontrollable zero eigenvalues of A - lambda*E, while
%         Afnc - lambda*Efnc contains the uncontrollable nonzero finite
%         and infinite eigenvalues of A - lambda*E.
%
%   For jobcon = 1, the pencil Anc - lambda*Enc has the form
%
%      Anc - lambda*Enc = Afnc - lambda*Efnc ,
%
%   where the regular pencil Afnc - lambda*Efnc, with Efnc upper
%   triangular and nonsingular, contains the uncontrollable finite
%   eigenvalues of A - lambda*E.
%
%   For jobcon = 2, the pencil Anc - lambda*Enc has the form
%
%      Anc - lambda*Enc = Ainc - lambda*Einc ,
%
%   where the regular pencil Ainc - lambda*Einc, with Ainc upper
%   triangular and nonsingular, contains the uncontrollable nonzero
%   finite and infinite eigenvalues of A - lambda*E.
%
%   2) To compute orthogonal transformation matrices Q and Z which
%   reduce the n-th order descriptor system (A-lambda*E,B,C) to the form
%
%               ( Ano  * )             ( Eno  * )           ( Bno )
%      Q'*A*Z = (        ) ,  Q'*E*Z = (        ) ,  Q'*B = (     ) ,
%               ( 0   Ao )             ( 0   Eo )           ( Bo  )
%
%         C*Z = ( 0   Co ) ,                                         (6)
%
%   where the o-th order descriptor system (Ao-lambda*Eo,Bo,Co) is a
%   finite and/or infinite observable. The pencil Ano - lambda*Eno is
%   regular of order n-o and contains the unobservable finite and/or
%   infinite eigenvalues of the pencil A-lambda*E. The reduced order
%   descriptor system (Ao-lambda*Eo,Bo,Co) has the same
%   transfer-function matrix as the original system (A-lambda*E,B,C).
%   The left and/or right orthogonal transformations Q and Z performed
%   to reduce the system matrices can be optionally accumulated.
%
%   For jobobs = 0 or 2, the pencil ( Eo-lambda*Ao ) has full column
%                                   (      Co      )
%   rank o for all finite lambda and is in a staircase form with
%                    _      _            _      _
%                  ( Ek,k   Ek,k-1   ... Ek,2   Ek,1   )
%                  ( _      _            _      _      )
%      ( Eo )  =   ( Ek-1,k Ek-1,k-1 ... Ek-1,2 Ek-1,1 ) ,           (7)
%      ( Co )      (     ...         ... _      _      )
%                  (  0       0      ... E1,2   E1,1   )
%                  (                            _      )
%                  (  0       0      ... 0      E0,1   )
%                    _          _      _
%                  ( Ak,k  ...  Ak,2   Ak,1 )
%                  (       ...  _      _    )
%        Ao    =   (   0   ...  A2,2   A2,1 ) ,                      (8)
%                  (                   _    )
%                  (   0   ...    0    A1,1 )
%         _
%   where Ei-1,i is a ctau(i-1)-by-ctau(i) full column rank matrix
%                          _
%   (with ctau(0) = p) and Ai,i is a ctau(i)-by-ctau(i) upper
%   triangular matrix.
%
%   For jobobs = 1, the pencil ( Ao-lambda*Eo ) has full column
%                              (      Co      )
%   rank o for all finite lambda and is in a staircase form with
%                    _      _            _      _
%                  ( Ak,k   Ak,k-1   ... Ak,2   Ak,1   )
%                  ( _      _            _      _      )
%      ( Ao )  =   ( Ak-1,k Ak-1,k-1 ... Ak-1,2 Ak-1,1 ) ,           (9)
%      ( Co )      (     ...         ... _      _      )
%                  (  0       0      ... A1,2   A1,1   )
%                  (                            _      )
%                  (  0       0      ... 0      A0,1   )
%                    _          _      _
%                  ( Ek,k  ...  Ek,2   Ek,1 )
%                  (       ...  _      _    )
%        Eo    =   (   0   ...  E2,2   E2,1 ) ,                     (10)
%                  (                   _    )
%                  (   0   ...    0    E1,1 )
%         _
%   where Ai-1,i is a ctau(i-1)-by-ctau(i) full column rank matrix
%                          _
%   (with ctau(0) = p) and Ei,i is a ctau(i)-by-ctau(i) upper
%   triangular matrix.
%
%   For jobobs = 0, the (n-o)-by-(n-o) regular pencil Ano - lambda*Eno
%   has the form
%
%                         ( Afno - lambda*Efno         *          )
%      Ano - lambda*Eno = (                                       ) ,
%                         (        0           Aino - lambda*Eino )
%
%   where:
%     1) the ino-by-ino regular pencil Aino - lambda*Eino, with Aino
%        upper triangular and nonsingular, contains the unobservable
%        infinite eigenvalues of A - lambda*E;
%     2) the (n-o-ino)-by-(n-o-ino) regular pencil Afno - lambda*Efno,
%        with Efno upper triangular and nonsingular, contains the
%        unobservable finite eigenvalues of A - lambda*E.
%
%   Note: The significance of the two diagonal blocks can be
%         interchanged by calling the gateway with the
%         arguments A and E interchanged. In this case,
%         Aino - lambda*Eino contains the unobservable zero
%         eigenvalues of A - lambda*E, while Afno - lambda*Efno
%         contains the unobservable nonzero finite and infinite
%         eigenvalues of A - lambda*E.
%
%   For jobobs = 1, the pencil Ano - lambda*Eno has the form
%
%      Ano - lambda*Eno = Afno - lambda*Efno ,
%
%   where the regular pencil Afno - lambda*Efno, with Efno upper
%   triangular and nonsingular, contains the unobservable finite
%   eigenvalues of A - lambda*E.
%
%   For jobobs = 2, the pencil Ano - lambda*Eno has the form
%
%      Ano - lambda*Eno = Aino - lambda*Eino ,
%
%   where the regular pencil Aino - lambda*Eino, with Aino upper
%   triangular and nonsingular, contains the unobservable nonzero
%   finite and infinite eigenvalues of A - lambda*E.
%
%   3) To find a reduced (controllable, observable, or irreducible)
%   descriptor representation (Ar-lambda*Er,Br,Cr) for an original
%   descriptor representation (A-lambda*E,B,C). The pencil Ar-lambda*Er
%   is in an upper block Hessenberg form, with either Ar or Er upper
%   triangular.
%
%   Description of input parameters:
%   task   - integer specifying the computations to be performed.
%            task = 1 :  compute controllability staircase form (1);
%            task = 2 :  compute observability staircase form (6);
%            task = 3 :  compute controllable, observable, or
%                        irreducible form.
%   A      - the n-by-n state dynamics matrix A.
%   E      - the n-by-n descriptor matrix E.
%   B      - the n-by-m input/state matrix B.
%   C      - the p-by-n state/output matrix C.
%   flag   - (optional) real vector of length 4 specifying various
%            options, depending on task.
%
%            For task = 1, flag contains:
%            flag(1) = jobcon : indicates what to do, as follows:
%               jobcon = 0 :  separate both finite and infinite
%                             uncontrollable eigenvalues;
%               jobcon = 1 :  separate only finite uncontrollable
%                             eigenvalues;
%               jobcon = 2 :  separate only nonzero finite and infinite
%                             uncontrollable eigenvalues.
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
%            flag(4) = tol : see below.
%            Default :  flag = [ 0; 0; 0; 0 ].
%
%            For task = 2, flag contains:
%            flag(1) = jobobs : indicates what to do, as follows:
%               jobobs = 0 :  separate both finite and infinite
%                             unobservable eigenvalues;
%               jobobs = 1 :  separate only finite unobservable
%                             eigenvalues;
%               jobobs = 2 :  separate only nonzero finite and infinite
%                             unobservable eigenvalues.
%            flag(2) = compq : see task = 1 above.
%            flag(3) = compz : see above.
%            flag(4) = tol : see below.
%            Default :  flag = [ 0; 0; 0; 0 ].
%
%            For task = 3, flag contains:
%            flag(1) = job : indicates what to do, as follows:
%               job = 0 :  remove both the uncontrollable and
%                          unobservable parts to get an irreducible
%                          descriptor representation;
%               job = 1 :  remove the uncontrollable part only to get a
%                          controllable descriptor representation;
%               job = 2 :  remove the unobservable part only to get an
%                          observable descriptor representation.
%            flag(2) = systyp : indicates the type of descriptor system
%            algorithm to be applied according to the assumed
%            transfer-function matrix as follows:
%               systyp = 0 :  rational transfer-function matrix;
%               systyp = 1 :  proper (standard) transfer-function
%                             matrix;
%               systyp = 2 :  polynomial transfer-function matrix.
%            flag(3) = equil : specifies whether the user wishes to
%            preliminarily scale the system (A-lambda*E,B,C) as follows:
%               equil = 0 :  perform scaling;
%               equil = 1 :  do not perform scaling.
%            flag(4) = tol : see below.
%            Default :  flag = [ 0; 0; 0; 0 ].
%
%            tol is a real scalar indicating the tolerance to be used
%            in rank determinations when transforming (A-lambda*E,B,C),
%            or (A-lambda*E, B) if task = 1, or (A'-lambda*E',C')', if
%            task = 2. If tol > 0, then the given value of tol is used
%            as a lower bound for reciprocal condition numbers in rank
%            determinations; a (sub)matrix whose estimated condition
%            number is less than 1/tol is considered to be of full rank.
%            If tol <= 0, the default tolerance toldef = eps*n*n is
%            used instead, where eps is the machine precision. tol < 1.
%   Q1     - (optional) if compq = 2 the n-by-n given orthogonal
%            matrix Q1.
%   Z1     - (optional) if compz = 2 the n-by-n given orthogonal
%            matrix Z1.
%
%   Description of output parameters:
%   Ao     - If task = 1, the n-by-n transformed state matrix Q'*A*Z,
%
%                         ( Ac   *  )
%                Q'*A*Z = (         ) ,
%                         ( 0   Anc )
%
%            where Ac is c-by-c and Anc is (n-c)-by-(n-c).
%            If jobcon = 1, the matrix ( Bc Ac ) is in the
%            controllability staircase form (4).
%            If jobcon = 0 or 2, the submatrix Ac is upper triangular.
%            If jobcon = 0, the Anc matrix has the form
%
%                      ( Ainc   *  )
%                Anc = (           ) ,
%                      (  0   Afnc )
%
%            where the inc-by-inc matrix Ainc is nonsingular and upper
%            triangular.
%            If jobcon = 2, Anc is nonsingular and upper triangular.
%
%            If task = 2, the n-by-n transformed state matrix Q'*A*Z,
%
%                         ( Ano  *  )
%                Q'*A*Z = (         ) ,
%                         ( 0    Ao )
%
%            where Ao is o-by-o and Ano is (n-o)-by-(n-o).
%            If jobobs = 1, the matrix ( Ao ) is in the observability
%                                      ( Co )
%            staircase form (9).
%            If jobobs = 0 or 2, the submatrix Ao is upper triangular.
%            If jobobs = 0, the submatrix Ano has the form
%
%                      ( Afno   *  )
%                Ano = (           ) ,
%                      (  0   Aino )
%
%            where the ino-by-ino matrix Aino is nonsingular and upper
%            triangular.
%            If jobobs = 2, Ano is nonsingular and upper triangular.
%
%            If task = 3, the nr-by-nr reduced order state matrix Ar of
%            an irreducible, controllable, or observable realization for
%            the original system, depending on the value of job,
%            job = 0, job = 1, or job = 2, respectively.
%            The matrix Ar is upper triangular if systyp = 0 or 2.
%            If systyp = 1 and job = 1, the matrix [Br Ar] is in a
%            controllable staircase form.
%            If systyp = 1 and job = 0 or 2, the matrix ( Ar ) is in an
%                                                       ( Cr )
%            observable staircase form.
%            The block structure of staircase forms is contained
%            in the leading infred(7) elements of the vector sizes.
%   Eo     - If task = 1, the n-by-n transformed descriptor matrix
%            Q'*E*Z,
%
%                        ( Ec   *  )
%               Q'*E*Z = (         ) ,
%                        ( 0   Enc )
%
%            where Ec is c-by-c and Enc is (n-c)-by-(n-c).
%            If jobcon = 0 or 2, the matrix ( Bc Ec ) is in the
%            controllability staircase form (2).
%            If jobcon = 1, the submatrix Ec is upper triangular.
%            If jobcon = 0, the Enc matrix has the form
%
%                      ( Einc   *  )
%                Enc = (           ) ,
%                      (  0   Efnc )
%
%            where the inc-by-inc matrix Einc is nilpotent and the
%            (n-c-inc)-by-(n-c-inc) matrix Efnc is nonsingular and
%            upper triangular.
%            If jobcon = 1, Enc is nonsingular and upper triangular.
%
%            If task = 2, the n-by-n transformed descriptor matrix
%            Q'*E*Z,
%
%                         ( Eno  *  )
%                Q'*E*Z = (         ) ,
%                         ( 0    Eo )
%
%            where Eo is o-by-o and Eno is (n-o)-by-(n-o).
%            If jobobs = 0 or 2, the matrix ( Eo ) is in the
%                                           ( Co )
%            observability staircase form (7).
%            If jobobs = 1, the submatrix Eo is upper triangular.
%            If jobobs = 0, the Eno matrix has the form
%
%                      ( Efno   *  )
%                Eno = (           ) ,
%                      (  0   Eino )
%
%            where the ino-by-ino matrix Eino is nilpotent and the
%            (n-o-ino)-by-(n-o-ino) matrix Efno is nonsingular and
%            upper triangular.
%            If jobobs = 1, Eno is nonsingular and upper triangular.
%
%            If task = 3, the nr-by-nr reduced order descriptor matrix
%            Er of an irreducible, controllable, or observable
%            realization for the original system, depending on the value
%            of job, job = 0, job = 1, or job = 2, respectively.
%            The resulting Er has infred(6) nonzero sub-diagonals.
%            If at least for one k = 1,...,4, infred(k) >= 0, then the
%            resulting Er is structured being either upper triangular
%            or block Hessenberg, in accordance to the last
%            performed order reduction phase (see Method).
%            The block structure of staircase forms is contained
%            in the leading infred(7) elements of the vector sizes.
%   Bo     - If task = 1 or task = 2, the n-by-m transformed state/input
%            matrix Q'*B. If task = 1, this matrix has the form
%
%                       ( Bc )
%                Q'*B = (    ) ,
%                       ( 0  )
%
%            where Bc is c-by-m.
%            For jobcon = 0 or 2, the matrix ( Bc Ec ) is in the
%            controllability staircase form (2).
%            For jobcon = 1, the matrix ( Bc Ac ) is in the
%            controllability staircase form (4).
%            If task = 3, the nr-by-m reduced input matrix Br of an
%            irreducible, controllable, or observable realization for
%            the original system, depending on the value of job,
%            job = 0, job = 1, or job = 2, respectively.
%            If job = 1, only the first sizes(1) rows of B are nonzero.
%   Co     - If task = 1, the p-by-n transformed matrix C*Z.
%            If task = 2, the p-by-n transformed matrix
%
%                C*Z = (  0   Co ) ,
%
%            where Co is p-by-o.
%            If jobobs = 0 or 2, the matrix ( Eo ) is in the
%                                           ( Co )
%            observability staircase form (7).
%            If jobobs = 1, the matrix ( Ao ) is in the observability
%                                      ( Co )
%            staircase form (9).
%            If task = 3, the p-by-nr transformed state/output matrix Cr
%            of an irreducible, controllable, or observable realization
%            for the original system, depending on the value of job,
%            job = 0, job = 1, or job = 2, respectively.
%            If job = 0, or job = 2, only the last sizes(1) columns
%            (in the first nr columns) of C are nonzero.
%   Q      - the n-by-n orthogonal matrix Q.
%            If compq = 1, Q' is the product of the transformations
%            which are applied to A, E, and B on the left.
%            If compq = 2, Q is the orthogonal matrix product Q1*Q.
%   Z      - the n-by-n orthogonal matrix Z.
%            If compz = 1, Z is the product of the transformations
%            which are applied to A, E, and C on the right.
%            If compz = 2, Z is the orthogonal matrix product Z1*Z.
%   orders - (optional) if task < 3, integer vector of length 3
%            indicating the orders of the subsystems and the number of
%            full rank blocks.
%            If task = 1, orders contains c, inc, and nrb:
%               c is the order of the reduced matrices Ac and Ec, and
%            the number of rows of reduced matrix Bc; also the order
%            of the controllable part of the pair (A-lambda*E,B).
%               For jobcon = 0, inc is the order of the reduced matrices
%            Ainc and Einc, and also the number of uncontrollable
%            infinite eigenvalues of the pencil A - lambda*E.
%            For jobcon <> 0, inc has no significance and it is set to
%            zero.
%                                                            _
%               nrb is the number k, of full row rank blocks Ei,i in the
%            staircase form of the pencil (Bc Ec-lambda*Ac) in (2) and
%                                                                  _
%            (3) (for jobcon = 0 or 2), or of full row rank blocks Ai,i
%            in the staircase form of the pencil (Bc Ac-lambda*Ec) in
%            (4) and (5) (for jobcon = 1).
%
%            If task = 2, orders contains o, ino, and nrb:
%               o is the order of the reduced matrices Ao and Eo, and
%            the number of columns of reduced matrix Co; also the order
%            of the observable part of the pair (C, A-lambda*E).
%               For jobobs = 0, ino is the order of the reduced matrices
%            Aino and Eino, and also the number of unobservable
%            infinite eigenvalues of the pencil A - lambda*E.
%            For jobobs <> 0, ino has no significance and it is set to
%            zero.
%                                                               _
%               nrb is the number k, of full column rank blocks Ei-1,i
%            in the staircase form of the pencil (Eo-lambda*Ao) in (7)
%                                                (    Co      )
%            and (8) (for jobobs = 0 or 2), or of full column rank
%                   _
%            blocks Ai-1,i in the staircase form of the pencil
%            (Ao-lambda*Eo) in (9) and (10) (for jobobs = 1).
%            (     Co     )
%   infred - (optional) if task = 3, integer array of dimension 7,
%            containing information on performed reduction and on
%            structure of resulting system matrices as follows:
%            infred(k) >= 0 (k = 1, 2, 3, or 4) if Phase k of reduction
%                           (see Method) has been performed. In this
%                           case, infred(k) is the achieved order
%                           reduction in Phase k.
%            infred(k) < 0  (k = 1, 2, 3, or 4) if Phase k was not
%                           performed.
%            infred(5)  -   the number of nonzero sub-diagonals of A.
%            infred(6)  -   the number of nonzero sub-diagonals of E.
%            infred(7)  -   the number of blocks in the resulting
%                           staircase form at last performed reduction
%                           phase. The block dimensions are contained
%                           in the first infred(7) elements of sizes.
%   sizes -  (optional) if task = 1 or 2, integer nrb-vector containing
%            rtau, if task = 1, or ctau, if task = 2.
%               rtau(i), for i = 1, ..., nrb, is the row dimension of
%                                    _         _
%            the full row rank block Ei,i-1 or Ai,i-1 in the staircase
%            form (2) or (4) for jobcon = 0 or 2, or for jobcon = 1,
%            respectively.
%               ctau(i), for i = 1, ..., nrb, is the column dimension
%                                          _         _
%            of the full column rank block Ei-1,i or Ai-1,i in the
%            staircase form (7) or (9) for jobobs = 0 or 2, or
%            for jobobs = 1, respectively.
%            If task = 3, integer vector of dimension infred(7), whose
%            elements contain the orders of the diagonal blocks of
%            Ar-lambda*Er.
%
%   Method
%   If task = 3, the order reduction is performed in 4 phases:
%   Phase 1: Eliminate all finite uncontrolable eigenvalues.
%            The resulting matrix ( Br Ar ) is in a controllable
%            staircase form, and Er is upper triangular.
%            This phase is performed if job = 0 or 1 and systyp = 0
%            or 1.
%   Phase 2: Eliminate all infinite and finite nonzero uncontrollable
%            eigenvalues. The resulting matrix ( Br Er ) is in a
%            controllable staircase form, and Ar is upper triangular.
%            This phase is performed if job = 0 or 1 and systyp = 0
%            or 2.
%   Phase 3: Eliminate all finite unobservable eigenvalues.
%            The resulting matrix ( Ar ) is in an observable
%                                 ( Cr )
%            staircase form, and Er is upper triangular.
%            This phase is performed if job = 0 or 2 and systyp = 0
%            or 1.
%   Phase 4: Eliminate all infinite and finite nonzero unobservable
%            eigenvalues. The resulting matrix ( Er ) is in an
%                                              ( Cr )
%            observable staircase form, and Ar is upper triangular.
%            This phase is performed if job = 0 or 2 and systyp = 0
%            or 2.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003.
%
% Revisions:
%   V. Sima, March 2004.
%
