#include "fintrf.h"
C KFILTUPD.F - Gateway function for computing a combined measurement and
C              time update of various Kalman filters for discrete-time
C              systems, using SLICOT routines FB01QD, FB01RD, FB01SD,
C              FB01TD, FB01VD, and FD01AD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2003-2020 NICONET e.V.
C
C Matlab call:
C   [(S1,K,Rin)(Sinv1,X1,Qin,E)(P1,K,Rin)(,rcnd)
C    (Ef1,Xf1,Eb1,Cs1,Sn1,YQ1,Err,Salph)] =
C          Kfiltupd(task(,opt)(,S,A,B,C(,Q),R)
C                             (,Sinv,Ainv,B,C,Qinv(,Rinv),X,Rinvy,Z)
C                             (,Sinv,Ainv,AinvB,C,Qinv(,Rinv),X,Rinvy,Z)
C                             (,P,A,B,C,Q,R)
C                             (,L,Lambda,Si,Ef,Xf,Eb,Cs,Sn(,YQ))
C                             (,tol))
C
C   [S1(,K,Rin,rcnd)] = Kfiltupd(-1/1,[WithK,BxQ],S,A,B,C(,Q),R(,tol))
C   [Sinv1(,X1,Qin,E,rcnd)] =
C                       Kfiltupd(-2,[WithX,AixB,RixC],Sinv,Ainv,B,C,Qinv
C                                (,Rinv),X,Rinvy,Z(,tol))
C   [Sinv1(,X1,Qin,E,rcnd)] =
C                       Kfiltupd(2,[WithX,RixC],Sinv,Ainv,AinvB,C,Qinv
C                                (,Rinv),X,Rinvy,Z(,tol))
C   [P1(,K,Rin,rcnd)] = Kfiltupd(0,P,A,B,C,Q,R(,tol))
C   [Ef1,Xf1,Eb1,Cs1,Sn1(,YQ1),Err,Salph] =
C                       Kfiltupd(3,[WithF,printw],L,Lambda,Si,Ef,Xf,Eb,
C                                Cs,Sn(,YQ))
C
C Purpose:
C   To compute a combined measurement and time update of one iteration
C   of the time-varying or time-invariant Kalman filter in the square
C   root covariance or information form, or of the conventional Kalman
C   filter. One time update of the recursive least-squares filtering
C   problem can also be computed, using a fast QR-decomposition based
C   approach.
C
C Input parameters:
C   task   - scalar indicating the calculation to be performed:
C            =-2 :  calculate a combined measurement and time update of
C                   one iteration of the time-varying Kalman filter in
C                   the square root information form;
C            =-1 :  calculate a combined measurement and time update of
C                   one iteration of the time-varying Kalman filter in
C                   the square root covariance form;
C            = 0 :  calculate one recursion of the conventional Kalman
C                   filter equations;
C            = 1 :  calculate a combined measurement and time update of
C                   one iteration of the time-invariant Kalman filter in
C                   the square root covariance form;
C            = 2 :  calculate a combined measurement and time update of
C                   one iteration of the time-invariant Kalman filter in
C                   the square root information form;
C            = 3 :  calculate one time update of the recursive least-
C                   squares filtering problem.
C   opt    - (optional) integer scalar containing option values.
C          - If |task| = 1, opt has length 2, opt = [WithK,BxQ].
C            WithK indicates whether the Kalman filter gain matrix K
C            should be computed, as follows:
C            = 0 :  K is not required;
C            = 1 :  K is computed and stored in array K.
C            Default: WithK = 0.
C            BxQ indicates how matrices B and Q are to be passed to the
C            gateway, as follows:
C            = 0 :  arrays B and Q must contain the matrices as
C                   described below;
C            = 1 :  B must contain the product B*sqrt(Q) and Q is not
C                   used.
C            Default: BxQ = 0.
C          - If task = -2, opt has length 3, opt = [WithX,AixB,RixC];
C            if task =  2, opt has length 2, opt = [WithX,RixC].
C            WithX indicates whether X    is to be computed, as follows:
C                                     i+1
C            = 0 :  X is not required;
C            = 1 :  X is computed and stored in array X.
C            Default: WithX = 0.          -1
C            AixB indicates how matrices A   and B  are to be passed to
C                                         i       i
C            the gateway, as follows:
C            = 0 :  arrays Ainv and B must contain the matrices as
C                   described below;                    -1
C            = 1 :  array Ainv must contain the matrix A   and the array
C                                               -1      i
C                   B must contain the product A  B .
C                                               i  i
C            Default: AixB = 0.           -1/2
C            RixC indicates how matrices R     and C    are to be passed
C                                         i+1       i+1
C            to the gateway, as follows:
C            = 0 :  arrays Rinv and C must contain the matrices as
C                   described below;                         -1/2
C            = 1 :  array C must contain the matrix product R    C
C                                                            i+1  i+1
C
C            and array Rinv is not used.
C            Default: RixC = 0.
C          - If task = 3, opt has length 2, opt = [WithF,printw].
C            WithF indicates whether both prediction and filtering parts
C            are to be applied, as follows:
C            = 0 :  only the prediction section is to be applied;
C            = 1 :  both prediction and filtering parts are to be
C                   applied.
C            Default: WithF = 0.
C            printw is a switch for printing the warning messages.
C            = 0:  do not print warning messages;
C            = 1:  print warning messages.
C            Default: printw = 0.
C          - If task = 0, opt must not be specified as an input
C            parameter.
C   S      - the n-by-n lower triangular part contains S   , the square
C                                                       i-1
C            root (left Cholesky factor) of the state covariance matrix
C            at instant (i-1). The strict upper triangular part is not
C            used.
C   A      - the n-by-n state transition matrix.
C          - If task = -1 or task = 0, it contains A , the state
C                                                   i
C            transition matrix at instant i.
C          - If task = 1, it contains A, the state transition matrix of
C            the system in lower observer Hessenberg form (e.g., as
C            produced by SLICOT Library Routine TB01ND).
C   B      - the n-by-m input weight matrix B  at instant i.
C                                            i                1/2
C          - If |task| = 1, it contains B  (or the product B Q   ,  if
C                                        i                  i i
C            BxQ = opt(2) = 1).                            -1
C          - If task = -2, it contains B  (or the product A  B ,  if
C                                       i                  i  i
C            AixB = opt(2) = 1).
C   C      - if task = -1 or task = 0, the p-by-n output weight
C            matrix C , at instant i.
C                    i
C          - If task = 1, the p-by-n output weight matrix C of the
C            system in lower observer Hessenberg form (e.g., as produced
C            by SLICOT Library routine TB01ND).
C          - If |task| = 2, the p-by-n output weight matrix C   , (or
C                        -1/2                                i+1
C            the product R    C    if RixC = 1) at instant i+1.
C                         i+1  i+1
C   Q      - if |task| = 1 and BxQ = 0, the m-by-m lower triangular
C                           1/2
C            part contains Q   , the square root (left Cholesky factor)
C                           i
C            of the input (process) noise covariance matrix at
C            instant i.
C            The strict upper triangular part is not used.
C          - If |task| = 1 and BxQ = 1, Q must not be specified as an
C            input parameter.
C          - If task = 0, the m-by-m matrix Q , the process noise
C                                            i
C            covariance matrix at instant i.
C   R      - if |task| = 1, the p-by-p lower triangular part contains
C             1/2
C            R   , the square root (left Cholesky factor) of the output
C             i
C            (measurement) noise covariance matrix at instant i.
C            The strict upper triangular part of this array is not used.
C          - If task = 0, the p-by-p output noise covariance matrix
C            R , at instant i.
C             i
C   tol    - (optional) real scalar containing a tolerance value used to
C                                                           1/2
C            test for near singularity of a matrix: (RINOV )   , if
C                                                -1       i
C            |task| = 1 and WithK = opt(1) = 1; S   , if |task| = 2 and
C                                                i+1
C            WithX = opt(1) = 1; and (RINOV) , if task = 0.
C                                           i
C            If tol > 0, then the given value of tol is used as a
C            lower bound for the reciprocal condition number of that
C            matrix; a matrix whose estimated condition number is less
C            than 1/tol is considered to be nonsingular. If tol <= 0,
C            then ps*eps, is used instead, where ps is the product of
C            the matrix dimensions, and eps is the machine precision.
C          - If task = 3, tol must not be specified as an input
C            parameter.
C            Default: tol = 0.
C                                                       -1
C   Sinv   - the n-by-n upper triangular part contains S  , the inverse
C                                                       i
C            of the square root (right Cholesky factor) of the state
C            covariance matrix P    (hence the information square root)
C                               i|i
C            at instant i.
C            The strict lower triangular part of this array is not used.
C                                             -1
C   Ainv   - if task = -2, the n-by-n matrix A  , the inverse of the
C                                             i
C            state transition matrix of the system at instant i.
C            The matrix A is time-varying if task = -2.
C                                            -1
C          - If task = 2, the n-by-n matrix A  , the inverse of the
C            state transition matrix of the system in controller
C            Hessenberg form (e.g., as produced by SLICOT Library
C            Routine TB01MD).
C                               -1                   -1
C   AinvB  - the n-by-m matrix A  B, the product of A   and the input
C            weight matrix B of the system, in upper controller
C            Hessenberg form (e.g., as produced by SLICOT Library
C            Routine TB01MD).
C                                                       -1/2
C   Qinv   - the m-by-m upper triangular part contains Q    , the
C                                                       i
C            inverse of the covariance square root (right Cholesky
C            factor) of the process noise (hence the information square
C            root) at instant i.
C            The strict lower triangular part of this array is not used.
C   Rinv   - if RixC = 0, the p-by-p upper triangular part contains
C             -1/2
C            R    , the inverse of the covariance square root (right
C             i+1
C            Cholesky factor) of the output noise (hence the information
C            square root) at instant i+1.
C            The strict lower triangular part of this array is not used.
C            Otherwise, Rinv must not be specified as an input
C            parameter.
C   X      - the n-vector X , the estimated filtered state at instant i.
C                          i
C                          -1/2
C   Rinvy  - the p-vector R    Y   , the product of the upper triangular
C                          i+1  i+1
C                    -1/2
C            matrix R     and the measured output vector Y    at instant
C                    i+1                                  i+1
C            i+1.
C   Z      - the m-vector Z , the mean value of the state process noise
C                          i
C            at instant i.
C   P      - the n-by-n state covariance matrix P      at instant (i-1).
C                                                i|i-1
C            The upper triangular part only is needed.
C   L      - integer scalar specifying the length of the impulse
C            response of the equivalent transversal filter model.
C            L >= 1.
C   Lambda - real scalar specifying the square root of the forgetting
C            factor. For tracking capabilities and exponentially stable
C            error propagation, Lambda < 1.0 (strict inequality) should
C            be used.  0.0 < Lambda <= 1.0.
C   Si     - Si(1) = Ui contains the input sample at instant i.
C            If WithF = 1, Si(2) = Yi contains the reference sample at
C            instant i. Otherwise, Yi is not used.
C            (The situation just before and just after the call of the
C            gateway are denoted by instant (i-1) and instant i,
C            respectively.)
C   Ef     - the square root of exponentially weighted forward
C            prediction error energy at instant (i-1).  Ef >= 0.0.
C   Xf     - L-vector containing the transformed forward prediction
C            variables at instant (i-1).
C   Eb     - (L+1)-vector; the leading L elements must contain the
C            normalized a posteriori backward prediction error residuals
C            of orders zero through L-1, respectively, at instant (i-1),
C            and Eb(L+1) must contain the square-root of the so-called
C            "conversion factor" at instant (i-1).
C   Cs     - L-vector containing the cosines of the rotation angles used
C            in time updates, at instant (i-1).
C   Sn     - L-vector containing the sines of the rotation angles used
C            in time updates, at instant (i-1).
C   YQ     - if WithF = 1, L-vector containing the orthogonally
C            transformed reference vector at instant (i-1). These
C            elements are also the tap multipliers of an equivalent
C            normalized lattice least-squares filter.
C            Otherwise, this parameter is not used.
C
C Output parameters:
C   S1     - the n-by-n lower triangular part contains S , the square
C                                                       i
C            root (left Cholesky factor) of the state covariance matrix
C            at instant i.
C   K      - if |task| = 1 and WithK = 1, or if task = 0, the n-by-p
C            Kalman filter gain matrix K  at instant i. If |task| = 1
C                                       i
C            and WithK = 0, the n-by-p matrix AK , the product of the
C                                               i
C            input matrix A and the Kalman filter gain matrix at
C            instant i.
C   Rin    - if |task| = 1, the p-by-p lower triangular part contains
C                    1/2
C            (RINOV )   , the square root (left Cholesky factor) of the
C                  i
C            covariance matrix of the innovations at instant i.
C          - If task = 0, the p-by-p upper triangular part contains
C                    1/2
C            (RINOV )   , the square root (left Cholesky factor) of the
C                  i
C            covariance matrix of the innovations at instant i.
C   rcnd   - (optional) if |task| = 1 and WithK = 1, an estimate of the
C            reciprocal of the condition number (in the 1-norm) of
C                    1/2
C            (RINOV )   .
C                  i
C          - If task = 0, an estimate of the reciprocal condition number
C            (in the 1-norm) of RINOV .
C                                    i
C          - If |task| = 2 and WithX = 1, an estimate of the reciprocal
C                                                        -1
C            of the condition number (in the 1-norm) of S   .
C                                                        i+1
C            If |task| = 1 and WithK = 0, or |task| = 2 and WithX = 0,
C            rcnd is not computed, but set to 1.
C                                                       -1
C   Sinv1  - the n-by-n upper triangular part contains S   , the inverse
C                                                       i+1
C            of the square root (right Cholesky factor) of the state
C            covariance matrix P        (hence the information square
C                               i+1|i+1
C            root) at instant i+1.
C   X1     - if WithX = 1, the n-vector X   , the estimated filtered
C                                        i+1
C            state at instant i+1.
C                                        -1
C          - If WithX = 0, the n-vector S   X   .
C                                        i+1 i+1
C                                                              -1/2
C   Qin    - the m-by-m upper triangular part contains (QINOV )    , the
C                                                            i
C            inverse of the covariance square root (right Cholesky
C            factor) of the process noise (hence the information square
C            root) at instant i.
C   E      - the p-vector E   , the estimated error at instant i+1.
C                          i+1
C   P1     - the n-by-n upper triangular part contains the upper
C            triangular part of the state covariance matrix P     ,
C                                                            i+1|i
C            at instant i.
C   Ef1    - the square root of the exponentially weighted forward
C            prediction error energy at instant i.
C   Xf1    - L-vector containing the transformed forward prediction
C            variables at instant i.
C   Eb1    - (L+1)-vector containing the normalized a posteriori
C            backward prediction error residuals, plus the square root
C            of the conversion factor at instant i.
C   Cs1    - L-vector containing the cosines of the rotation angles at
C            instant i.
C   Sn1    - L-vector containing the sines of the rotation angles at
C            instant i.
C   YQ1    - if WithF = 1, L-vector containing the orthogonally
C            transformed reference vector at instant i.
C            If WithF = 0, YQ1 must not be specified as an output
C            parameter.
C   Err    - Err(1) = Ep contains the a posteriori forward prediction
C            error residual.
C            If WithF = 1, Err(2) = Eo contains the a posteriori output
C            error residual from the least-squares filter at instant i.
C   Salph  - L-vector: the element Salph(i), i=1,...,L, contains the
C            opposite of the i-(th) reflection coefficient for the
C            least-squares normalized lattice predictor (whose value is
C            -Salph(i)).
C
C Method:
C   If task = -2, the gateway performs one recursion of the square root
C   information filter algorithm, summarized as follows:
C
C     |    -1/2             -1/2    |     |         -1/2             |
C     |   Q         0      Q    Z   |     | (QINOV )     *     *     |
C     |    i                i    i  |     |       i                  |
C     |                             |     |                          |
C     |  -1 -1     -1 -1    -1      |     |             -1    -1     |
C   T | S  A  B   S  A     S  X     |  =  |    0       S     S   X   |
C     |  i  i  i   i  i     i  i    |     |             i+1   i+1 i+1|
C     |                             |     |                          |
C     |           -1/2      -1/2    |     |                          |
C     |    0     R    C    R    Y   |     |    0         0     E     |
C     |           i+1  i+1  i+1  i+1|     |                     i+1  |
C
C                (Pre-array)                      (Post-array)
C
C   where T is an orthogonal transformation triangularizing the
C                      -1/2
C   pre-array, (QINOV )     is the inverse of the covariance square
C                    i
C   root (right Cholesky factor) of the process noise (hence the
C   information square root) at instant i, and E    is the estimated
C                                               i+1
C   error at instant i+1.
C
C   The inverse of the corresponding state covariance matrix P
C                                                             i+1|i+1
C   (hence the information matrix I) is then factorized as
C
C                 -1         -1     -1
C      I       = P       = (S   )' S
C       i+1|i+1   i+1|i+1    i+1   i+1
C
C   and one combined time and measurement update for the state is
C   given by X   .
C             i+1
C
C   The triangularization is done entirely via Householder
C   transformations exploiting the zero pattern of the pre-array.
C
C   If task = 2, the gateway performs one recursion of the square root
C   information filter algorithm, summarized as follows:
C
C     |    -1/2             -1/2    |     |         -1/2             |
C     |   Q         0      Q    Z   |     | (QINOV )     *     *     |
C     |    i                i    i  |     |       i                  |
C     |                             |     |                          |
C     |           -1/2      -1/2    |     |             -1    -1     |
C   T |    0     R    C    R    Y   |  =  |    0       S     S   X   |
C     |           i+1  i+1  i+1  i+1|     |             i+1   i+1 i+1|
C     |                             |     |                          |
C     |  -1 -1     -1 -1    -1      |     |                          |
C     | S  A  B   S  A     S  X     |     |    0         0     E     |
C     |  i         i        i  i    |     |                     i+1  |
C
C                 (Pre-array)                      (Post-array)
C
C   where T is an orthogonal transformation triangularizing the
C                      -1/2
C   pre-array, (QINOV )     is the inverse of the covariance square
C                    i
C   root (right Cholesky factor) of the process noise (hence the
C                                               -1  -1
C   information square root) at instant i and (A  ,A  B) is in upper
C   controller Hessenberg form.
C
C   An example of the pre-array is given below (where n = 6, m = 2,
C   and p = 3):
C
C       |x x |             | x|
C       |  x |             | x|
C       _______________________
C       |    | x x x x x x | x|
C       |    | x x x x x x | x|
C       |    | x x x x x x | x|
C       _______________________ .
C       |x x | x x x x x x | x|
C       |  x | x x x x x x | x|
C       |    | x x x x x x | x|
C       |    |   x x x x x | x|
C       |    |     x x x x | x|
C       |    |       x x x | x|
C
C   The inverse of the corresponding state covariance matrix P
C                                                             i+1|i+1
C   (hence the information matrix I) is then factorized as
C
C                  -1         -1    -1
C       I       = P       = (S   )' S
C        i+1|i+1   i+1|i+1    i+1   i+1
C
C   and one combined time and measurement update for the state is
C   given by X   .
C             i+1
C
C   The triangularization is done entirely via Householder
C   transformations exploiting the zero pattern of the pre-array.
C
C   If task = -1, the gateway performs one recursion of the square root
C   covariance filter algorithm, summarized as follows:
C
C    |  1/2                      |     |         1/2          |
C    | R      C x S      0       |     | (RINOV )     0     0 |
C    |  i      i   i-1           |     |       i              |
C    |                      1/2  | T = |                      |
C    | 0      A x S    B x Q     |     |     AK       S     0 |
C    |         i   i-1  i   i    |     |       i       i      |
C
C        (Pre-array)                      (Post-array)
C
C   where T is an orthogonal transformation triangularizing the
C   pre-array.
C
C   The state covariance matrix P    is factorized as
C                                i|i-1
C      P     = S  S'
C       i|i-1   i  i
C
C   and one combined time and measurement update for the state X
C                                                               i|i-1
C   is given by
C
C      X     = A X      + K (Y - C X     ),
C       i+1|i   i i|i-1    i  i   i i|i-1
C
C                        -1/2
C   where K = AK (RINOV )     is the Kalman filter gain matrix and Y
C          i    i      i                                            i
C   is the observed output of the system.
C
C   The triangularization is done entirely via Householder
C   transformations exploiting the zero pattern of the pre-array.
C
C   If task = 1, the gateway performs one recursion of the square root
C   covariance filter algorithm, summarized as follows:
C
C    |  1/2                     |     |         1/2          |
C    | R      0        C x S    |     | (RINOV )     0     0 |
C    |  i                   i-1 |     |       i              |
C    |             1/2          | T = |                      |
C    | 0      B x Q    A x S    |     |     AK       S     0 |
C    |         i   i        i-1 |     |       i       i      |
C
C         (Pre-array)                      (Post-array)
C
C   where T is unitary and (A,C) is in lower observer Hessenberg form.
C
C   An example of the pre-array is given below (where n = 6, p = 2
C   and m = 3):
C
C        |x   |      | x          |
C        |x x |      | x x        |
C        |____|______|____________|
C        |    | x x x| x x x      |
C        |    | x x x| x x x x    | .
C        |    | x x x| x x x x x  |
C        |    | x x x| x x x x x x|
C        |    | x x x| x x x x x x|
C        |    | x x x| x x x x x x|
C
C   The corresponding state covariance matrix P      is then
C                                              i|i-1
C   factorized as
C
C       P     = S  S'
C        i|i-1   i  i
C
C   and one combined time and measurement update for the state X
C                                                               i|i-1
C   is given by
C
C       X     = A X      + K (Y - C X     )
C        i+1|i     i|i-1    i  i     i|i-1
C
C                        -1/2
C   where K = AK (RINOV )     is the Kalman filter gain matrix and Y
C          i    i      i                                            i
C   is the observed output of the system.
C
C   The triangularization is done entirely via Householder
C   transformations exploiting the zero pattern of the pre-array.
C
C   If task = 0, the gateway performs one recursion of the conventional
C   Kalman filter. The Kalman filter gain used at the i-th recursion
C   step is of the form
C
C                          -1
C      K  = P     C'  RINOV  ,
C       i    i|i-1 i       i
C
C   where RINOV  = C P     C' + R , and the state covariance matrix
C              i    i i|i-1 i    i
C
C   P      is updated by the discrete-time Riccati equation
C    i|i-1
C
C      P      = A  (P      - K C P     ) A'  + B Q B'.
C       i+1|i    i   i|i-1    i i i|i-1   i     i i i
C
C   Using these two updates, the combined time and measurement update
C   of the state X      is given by
C                 i|i-1
C
C      X      = A X      + A K (Y  - C X     ),
C       i+1|i    i i|i-1    i i  i    i i|i-1
C
C   where Y  is the new observation at step i.
C          i
C
C   If task = 3, the gateway performs one recursion of the recursive
C   least-squares filter. The output error Eo at instant i, denoted by
C   Eo(i), is the reference sample minus a linear combination of L
C   successive input samples:
C
C                      L-1
C      Eo(i) = Yi(i) - SUM h_j * Ui(j-i),
C                      j=0
C
C   where Yi(i) and Ui(i) are the scalar samples at instant i.
C   A least-squares filter uses those h_0,...,h_{L-1} which minimize
C   an exponentially weighted sum of successive output errors squared:
C
C       i
C      SUM [Lambda**(2(i-k)) * Eo(k)**2].
C      k=1
C
C   Each gateway call performs a time update of the least-squares
C   filter using a fast least-squares algorithm derived from a
C   QR decomposition. The algorithm does not compute the parameters
C   h_0,...,h_{L-1} from the above formula, but instead delivers the
C   parameters of an equivalent normalized least-squares lattice filter,
C   which are available from the arrays Salph (reflection coefficients)
C   and YQ (tap multipliers), as well as the exponentially weighted
C   input signal energy
C
C       i                                           L
C      SUM [Lambda**(2(i-k)) * Ui(k)**2] = Ef**2 + SUM Xf(j)**2.
C      k=1                                         j=1
C
C Comments for task = 3:
C   1.  For tracking capabilities and exponentially stable error
C       propagation, Lambda < 1.0 should be used.  Lambda is typically
C       chosen slightly less than 1.0 so that "past" data are
C       exponentially forgotten.
C   2.  Prior to the first gateway call, the variables must be
C       initialized. The following initial values are recommended:
C
C       Xf(i)   = 0.0,        i=1,...,L
C       Eb(i)   = 0.0,        i=1,...,L
C       Eb(L+1) = 1.0
C       Cs(i)   = 1.0,        i=1,...,L
C       Sn(i)   = 0.0,        i=1,...,L
C       YQ(i)   = 0.0,        i=1,...,L
C
C       Ef      = 0.0         (exact start)
C       Ef      = "small positive constant" (soft start).
C
C       Soft starts are numerically more reliable, but result in a
C       biased least-squares solution during the first few iterations.
C       This bias decays exponentially fast provided Lambda < 1.0.
C       If sigma is the standard deviation of the input sequence
C       Xi, then initializing Ef = sigma*1.0E-02 usually works well.
C   3.  The algorithm is backward consistent for all input sequences Xi,
C       and backward stable for persistently exciting input sequences,
C       assuming Lambda < 1.0.
C       If the condition of the signal is very poor (IWARN = 1), then
C       the results are not guaranteed to be reliable.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2003.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003,
C   Apr. 2009, Dec. 2010, July 2012, Aug. 2012.
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
      CHARACTER         JOBK, JOBX, JP, MULTAB, MULTBQ, MULTRC
      INTEGER           INFO, IWARN, L, LDA, LDAINB, LDAINV, LDB, LDC,
     $                  LDK, LDQ, LDQINV, LDR, LDRINV, LDS, LDSINV, LDW,
     $                  LDWORK, M, N, P
      DOUBLE PRECISION  EFOR, EOUT, EPOS, LAMBDA, TOL, XIN, YIN
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER, ALLOCATABLE ::          IWORK(:)
C
      mwPointer         A, AINV, AINVB, B, C, CTETA, DWORK, E, EPSBCK,
     $                  K, Q, QINV, R, RINV, RINVY, S, SALPH, SINV,
     $                  STETA, W, X, XF, YQ, Z
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           AIXB, BXQ, PRINTW, RIXC, WITHF, WITHK, WITHX
      CHARACTER*5       NAME
      CHARACTER*120     TEXT
      INTEGER           IOPT(3), IATASK, IP, IPC, IRCND, ISIZE, ITASK,
     $                  LOPT, LIWORK, LSI, MC, PC
      DOUBLE PRECISION  OP(3), RCND, TEMP
C
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C
C .. External subroutines ..
      EXTERNAL          DCOPY, FB01QD, FB01RD, FB01SD, FB01TD, FB01VD,
     $                  FD01AD, MA02KV
C
C .. Intrinsic functions ..
      INTRINSIC         ABS, MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.7 ) THEN
         CALL mexErrMsgTxt
     $        ( 'KFILTUPD requires at least 7 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'KFILTUPD requires at least 1 output argument' )
      END IF
C
C   task, (opt,) S(nxn), A(nxn), B(nxm), C(pxn), (Q(nxn)), (R(pxp),)
C   (X(n), Rinvy(p), Z(m)) (, tol)
C
C   task
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
      IF ( ITASK.LT.-2 .OR. ITASK.GT.3 ) THEN
         CALL mexErrMsgTxt
     $    ( 'ITASK has -2, -1, 0, 1, 2 or 3 the only admissible values')
      END IF
      IATASK = ABS( ITASK )
C
      IF ( ITASK.EQ.-2 ) THEN
         LOPT = 3
      ELSE
         LOPT = 2
      END IF
      IP = 2
C
C   opt
C
      IOPT(1) = 0
      IOPT(2) = 0
      IOPT(3) = 0
      IF ( ITASK.NE.0 ) THEN
         ISIZE = mxGetM( PRHS(2) )*mxGetN( PRHS(2) )
         IF ( ISIZE.GT.LOPT ) THEN
            WRITE( TEXT, '('' OPT must have at most '',I7,
     $                     '' entries'')' ) LOPT
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(2) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'OPT must be an integer vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), OP, ISIZE )
         IF( ISIZE.GE.1 ) IOPT(1) = OP(1)
         IF( ISIZE.GE.2 ) IOPT(2) = OP(2)
         IF( ISIZE.GE.3 ) IOPT(3) = OP(3)
         IF ( MIN( IOPT(1), IOPT(2), IOPT(3) ).LT.0 .OR.
     $        MAX( IOPT(1), IOPT(2), IOPT(3) ).GT.1 ) THEN
            CALL mexErrMsgTxt
     $   ( 'The entries of OPT have 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
      END IF
C
      IF ( IATASK.EQ.1 ) THEN
         WITHK = IOPT(1).EQ.1
         BXQ   = IOPT(2).EQ.1
      ELSE
         WITHK = .FALSE.
         BXQ   = .FALSE.
      END IF
C
      IF ( IATASK.EQ.2 ) THEN
         WITHX = IOPT(1).EQ.1
         IF ( ITASK.EQ.-2 ) THEN
            AIXB = IOPT(2).EQ.1
         ELSE
            AIXB = .TRUE.
         END IF
         RIXC  = IOPT(LOPT).EQ.1
      ELSE
         WITHX = .FALSE.
         AIXB  = .FALSE.
         RIXC  = .FALSE.
      END IF
C
      IF ( ITASK.EQ.0 )
     $   WITHK = .TRUE.
C
      IF ( ITASK.EQ.3 ) THEN
         WITHF  = IOPT(1).EQ.1
         PRINTW = IOPT(2).EQ.1
      ELSE
         WITHF  = .FALSE.
         PRINTW = .FALSE.
      END IF
C
      IF ( IATASK.GE.2 .AND. NRHS.LT.10 )
     $   CALL mexErrMsgTxt
     $      ( '(KFILTUPD requires at least 10 input arguments'')' )
C
      IF ( ITASK.LT.3 ) THEN
C
C   task < 3
C
         N = mxGetM( PRHS(IP) )
         M = mxGetN( PRHS(IP+2) )
         P = mxGetM( PRHS(IP+3) )
C
C   S/Sinv/P
C
         IF ( IATASK.EQ.1 ) THEN
            NAME = 'S'
         ELSE IF ( IATASK.EQ.2 ) THEN
            NAME = 'Sinv'
         ELSE
            NAME = 'P'
         END IF
C
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            WRITE( TEXT, '(A, '' must be a real matrix'')' ) NAME
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxGetN( PRHS(IP) ).NE.N ) THEN
            WRITE( TEXT, '(A, '' must be a square matrix'')' ) NAME
            CALL mexErrMsgTxt( TEXT )
         END IF
         IP = IP + 1
C
C   A/Ainv
C
         IF ( IATASK.EQ.2 ) THEN
            NAME = 'Ainv'
         ELSE
            NAME = 'A'
         END IF
C
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            WRITE( TEXT, '(A, '' must be a real matrix'')' ) NAME
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxGetM( PRHS(IP) ).NE.N ) THEN
            WRITE( TEXT, '(A, '' must have '', I7, '' rows'')' ) NAME, N
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxGetN( PRHS(IP) ).NE.N ) THEN
            WRITE( TEXT, '(A, '' must be a square matrix'')' ) NAME
            CALL mexErrMsgTxt( TEXT )
         END IF
         IP = IP + 1
C
C   B/AinvB
C
         IF ( ITASK.EQ.2 .OR. AIXB ) THEN
            NAME = 'AinvB'
         ELSE
            NAME = 'B'
         END IF
C
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            WRITE( TEXT, '(A, '' must be a real matrix'')' ) NAME
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxGetM( PRHS(IP) ).NE.N ) THEN
            WRITE( TEXT, '(A, '' must have '', I7, '' rows'')' ) NAME, N
            CALL mexErrMsgTxt( TEXT )
         END IF
         IP = IP + 1
C
C   C
C
         NAME = 'C'
C
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            WRITE( TEXT, '(A, '' must be a real matrix'')' ) NAME
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxGetN( PRHS(IP) ).NE.N ) THEN
            WRITE( TEXT, '(A, '' must have '', I7, '' columns'')' )
     $             NAME, N
            CALL mexErrMsgTxt( TEXT )
         END IF
         IP = IP + 1
C
C   Q/Qinv
C
         IF ( IATASK.EQ.2 ) THEN
            NAME = 'Qinv'
         ELSE IF ( .NOT.BXQ ) THEN
            NAME = 'Q'
         END IF
C
         IF ( .NOT.BXQ ) THEN
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               WRITE( TEXT, '(A, '' must be a real matrix'')' ) NAME
               CALL mexErrMsgTxt( TEXT )
            END IF
            IF ( mxGetM( PRHS(IP) ).NE.M .OR.
     $           mxGetN( PRHS(IP) ).NE.M ) THEN
               WRITE( TEXT, '(A, '' must have '', I7,
     $                           '' rows and columns'')' ) NAME, M
               CALL mexErrMsgTxt( TEXT )
            END IF
            IP = IP + 1
         END IF
C
C   R/Rinv
C
         IF ( IATASK.NE.2 ) THEN
            NAME = 'R'
         ELSE IF ( .NOT.RIXC ) THEN
            NAME = 'Rinv'
         END IF
C
         IF ( .NOT.RIXC ) THEN
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               WRITE( TEXT, '(A, '' must be a real matrix'')' ) NAME
               CALL mexErrMsgTxt( TEXT )
            END IF
            IF ( mxGetM( PRHS(IP) ).NE.P .OR.
     $           mxGetN( PRHS(IP) ).NE.P ) THEN
               WRITE( TEXT, '(A, '' must have '', I7,
     $                           '' rows and columns'')' ) NAME, P
               CALL mexErrMsgTxt( TEXT )
            END IF
            IP = IP + 1
         END IF
C
         IF ( IATASK.EQ.2 ) THEN
C
C   X
C
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'X must be a real vector' )
            END IF
            IF ( mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) ).NE.N ) THEN
               WRITE( TEXT, '('' X must have '', I7, '' entries'')' ) N
               CALL mexErrMsgTxt( TEXT )
            END IF
            IP = IP + 1
C
C   Rinvy
C
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'Rinvy must be a real vector' )
            END IF
            IF ( mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) ).NE.P ) THEN
               WRITE( TEXT, '('' Rinvy must have '', I7,
     $                        '' entries'')' ) P
               CALL mexErrMsgTxt( TEXT )
            END IF
            IP = IP + 1
C
C   Z
C
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'Z must be a real vector' )
            END IF
            IF ( mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) ).NE.M ) THEN
               WRITE( TEXT, '('' Z must have '', I7, '' entries'')' ) M
               CALL mexErrMsgTxt( TEXT )
            END IF
            IP = IP + 1
         END IF
C
         IF ( NRHS.GE.IP ) THEN
C
C   tol
C
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'TOL must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'TOL must be a real scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TOL, 1 )
         ELSE
            TOL = ZERO
         END IF
C
      ELSE
C
C   task = 3
C
C   L
C
         IF ( mxGetM( PRHS(3) ).NE.1 .OR.
     $        mxGetN( PRHS(3) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'L must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(3) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'L must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TEMP, 1 )
         L = TEMP
         IF ( L.LT.1 ) THEN
            CALL mexErrMsgTxt( 'L must satisfy L >= 1' )
         END IF
C
C   Lambda
C
         IF ( mxGetM( PRHS(4) ).NE.1 .OR.
     $        mxGetN( PRHS(4) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'LAMBDA must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(4) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'LAMBDA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), LAMBDA, 1 )
         IF ( LAMBDA.LE.ZERO .OR. LAMBDA.GT.ONE ) THEN
            CALL mexErrMsgTxt
     $         ( 'LAMBDA must satisfy 0.0 < LAMBDA <= 1.0' )
         END IF
C
C   Si
C
         ISIZE = mxGetM( PRHS(5) )*mxGetN( PRHS(5) )
         IF ( WITHF ) THEN
            LSI = 2
         ELSE
            LSI = 1
         END IF
         IF ( ISIZE.GT.LSI ) THEN
            WRITE( TEXT, '('' Si must have at most '',I1,
     $                     '' entries'')' ) LSI
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(5) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'Si must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), OP, ISIZE )
         IF( ISIZE.GE.1 ) XIN = OP(1)
         IF( ISIZE.GE.2 ) YIN = OP(2)
C
C   Ef
C
         IF ( mxGetM( PRHS(6) ).NE.1 .OR.
     $        mxGetN( PRHS(6) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'Ef must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'Ef must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), EFOR, 1 )
         IF ( EFOR.LT.ZERO ) THEN
            CALL mexErrMsgTxt
     $         ( 'Ef must satisfy Ef >= 0.0' )
         END IF
C
C   Xf
C
         IF ( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(7) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'Xf must be a real vector' )
         END IF
         IF ( mxGetM( PRHS(7) )*mxGetN( PRHS(7) ).NE.L ) THEN
            WRITE( TEXT, '('' Xf must have '', I7, '' entries'')' ) L
            CALL mexErrMsgTxt( TEXT )
         END IF
C
C   Eb
C
         IF ( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(8) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'Eb must be a real vector' )
         END IF
         IF ( mxGetM( PRHS(8) )*mxGetN( PRHS(8) ).NE.L+1 ) THEN
            WRITE( TEXT, '('' Eb must have '', I7, '' entries'')' ) L+1
            CALL mexErrMsgTxt( TEXT )
         END IF
C
C   Cs
C
         IF ( mxIsNumeric( PRHS(9) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(9) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'Cs must be a real vector' )
         END IF
         IF ( mxGetM( PRHS(9) )*mxGetN( PRHS(9) ).NE.L ) THEN
            WRITE( TEXT, '('' Cs must have '', I7, '' entries'')' ) L
            CALL mexErrMsgTxt( TEXT )
         END IF
C
C   Sn
C
         IF ( mxIsNumeric( PRHS(10) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(10) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'Sn must be a real vector' )
         END IF
         IF ( mxGetM( PRHS(10) )*mxGetN( PRHS(10) ).NE.L ) THEN
            WRITE( TEXT, '('' Sn must have '', I7, '' entries'')' ) L
            CALL mexErrMsgTxt( TEXT )
         END IF
C
         IF ( WITHF ) THEN
C
C   YQ
C
            IF ( mxIsNumeric( PRHS(11) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(11) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'YQ must be a real vector' )
            END IF
            IF ( mxGetM( PRHS(11) )*mxGetN( PRHS(11) ).NE.L ) THEN
               WRITE( TEXT, '('' YQ must have '', I7, '' entries'')' ) L
               CALL mexErrMsgTxt( TEXT )
            END IF
         END IF
      END IF
C
C   Set options.
C
      IF ( WITHK ) THEN
         JOBK = 'K'
      ELSE
         JOBK = 'N'
      END IF
C
      IF ( BXQ ) THEN
         MULTBQ = 'P'
      ELSE
         MULTBQ = 'N'
      END IF
C
      IF ( WITHX ) THEN
         JOBX = 'X'
      ELSE
         JOBX = 'N'
      END IF
C
      IF ( AIXB ) THEN
         MULTAB = 'P'
      ELSE
         MULTAB = 'N'
      END IF
C
      IF ( RIXC ) THEN
         MULTRC = 'P'
      ELSE
         MULTRC = 'N'
      END IF
C
      IF ( WITHF ) THEN
         JP = 'B'
      ELSE
         JP = 'P'
      END IF
C
C Determine the lenghts of working arrays.
C For task <> 0 or task <> 3, use a larger value for LDWORK for enabling
C calls of block algorithms in FB01QD, FB01RD, FB01SD, and FB01TD.
C
      IF ( IATASK.LE.2 ) THEN
         LDA = MAX( 1, N )
         LDB = LDA
         LDC = MAX( 1, P )
         LDS = LDA
         LDQ = MAX( 1, M )
         LDR = LDC
         LDK = LDA
         IF ( IATASK.EQ.2 ) THEN
            LDAINV = LDA
            IF ( ITASK.EQ.2 )
     $         LDAINB = LDA
            LDSINV = LDA
            LDQINV = LDQ
            IF ( .NOT.RIXC ) THEN
               LDRINV = LDR
               PC     = P
            ELSE
               LDRINV = 1
               PC     = 0
            END IF
         ELSE IF ( IATASK.EQ.1 ) THEN
            IF ( .NOT.BXQ ) THEN
               MC = M
            ELSE
               LDQ = 1
               MC  = 0
            END IF
         ELSE
            LDW = LDA
            MC  = M
         END IF
      END IF
C
      IF ( IATASK.EQ.1 ) THEN
         LDWORK = MAX( N*( P + N ) + 2*P, N*( N + M + 2 ) )
         IF ( WITHK ) THEN
            LDWORK = MAX( 2, 3*P, LDWORK )
            LIWORK = P
         ELSE
            LDWORK = MAX( 1, LDWORK )
            LIWORK = 1
         END IF
         IF ( ITASK.EQ.1 )
     $      LDWORK = MAX( N*( P + N + 1 ), LDWORK )
      ELSE IF ( IATASK.EQ.2 ) THEN
         LDWORK = N*( N + 2*M ) + 3*M
         IF ( ITASK.EQ.-2 ) THEN
            LDWORK = MAX( ( N + P )*( N + 1 ) + 2*N, LDWORK )
         ELSE
            LDWORK = MAX( ( N + P )*( N + 1 ) + N + MAX( N - 1, M + 1 ),
     $                    LDWORK )
         END IF
         IF ( WITHX ) THEN
            LDWORK = MAX( 2, 3*N, LDWORK )
            LIWORK = N
         ELSE
            LDWORK = MAX( 1, LDWORK )
            LIWORK = 1
         END IF
      ELSE IF ( ITASK.EQ.0 ) THEN
         LIWORK = P
         LDWORK = MAX( 1, P*N + 3*P, N*N, N*M )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      IF ( ITASK.LE.2 ) THEN
         ALLOCATE ( IWORK( LIWORK ) )
         C     = mxCalloc(  LDC*N, 8 )
         DWORK = mxCalloc( LDWORK, 8 )
         IF ( ITASK.LE.1 ) THEN
            B = mxCalloc( LDB*M, 8 )
            IF ( IATASK.LE.1  ) THEN
               A = mxCalloc( LDA*N,  8 )
               K = mxCalloc( LDK*P,  8 )
               Q = mxCalloc( LDQ*MC, 8 )
               R = mxCalloc( LDR*P,  8 )
               IF ( IATASK.EQ.1  ) THEN
                  S = mxCalloc( LDS*N, 8 )
               ELSE
                  W = mxCalloc( LDW*N, 8 )
               END IF
            END IF
         END IF
         IF ( IATASK.EQ.2 ) THEN
            AINV  = mxCalloc( LDAINV*N,  8 )
            E     = mxCalloc(        P,  8 )
            QINV  = mxCalloc( LDQINV*M,  8 )
            RINV  = mxCalloc( LDRINV*PC, 8 )
            RINVY = mxCalloc(        P,  8 )
            SINV  = mxCalloc( LDSINV*N,  8 )
            X     = mxCalloc(        N,  8 )
            Z     = mxCalloc(        M,  8 )
            IF ( ITASK.EQ.2 )
     $         AINVB = mxCalloc( LDAINB*M, 8 )
         END IF
      ELSE
         CTETA  = mxCalloc(   L, 8 )
         EPSBCK = mxCalloc( L+1, 8 )
         SALPH  = mxCalloc(   L, 8 )
         STETA  = mxCalloc(   L, 8 )
         XF     = mxCalloc(   L, 8 )
         IF ( WITHF )
     $      YQ = mxCalloc( L, 8 )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      IF ( ITASK.LE.2 ) THEN
         IF ( IATASK.EQ.2 ) THEN
            IP = 8
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ),
     $                             %VAL( SINV ), N*N )
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ),
     $                             %VAL( AINV ), N*N )
            IF ( ITASK.EQ.2 ) THEN
               CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ),
     $                                %VAL( AINVB ), N*M )
            ELSE
               CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ),
     $                                %VAL( B ), N*M )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ),
     $                             %VAL( QINV ), M*M )
            IF ( .NOT.RIXC ) THEN
               CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ),
     $                                %VAL( RINV ), P*P )
               IP = IP + 1
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ),
     $                             %VAL( X ), N )
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP+1) ),
     $                             %VAL( RINVY ), P )
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP+2) ),
     $                             %VAL( Z ), M )
            IPC = 6
         ELSE
            IF ( IATASK.EQ.1 ) THEN
               IP = 3
               CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ),
     $                                %VAL( S ), N*N )
            ELSE
               IP = 2
               CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ),
     $                                %VAL( W ), N*N )
            END IF
            IP  = IP + 1
            IPC = IP + 2
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ),
     $                             %VAL( A ), N*N )
            IP = IP + 1
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ),
     $                             %VAL( B ), N*M )
            IP = IP + 2
            IF ( .NOT.BXQ ) THEN
               CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ),
     $                                %VAL( Q ), M*M )
               IP = IP + 1
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ),
     $                             %VAL( R ), P*P )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IPC) ), %VAL( C ),
     $                          P*N )
C
      ELSE
C
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(7)  ), %VAL( XF ),
     $                          L )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(8)  ),
     $                          %VAL( EPSBCK ), L+1 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(9)  ),
     $                          %VAL( CTETA ), L )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(10) ),
     $                          %VAL( STETA ), L )
         IF ( WITHF )
     $      CALL mxCopyPtrToReal8( mxGetPr( PRHS(11) ),
     $                             %VAL( YQ ), L )
      END IF
C
C Do the actual computations.
C
      IRCND = 2
      IWARN = 0
      IF ( ITASK.EQ.-2 ) THEN
         CALL FB01SD( JOBX, MULTAB, MULTRC, N, M, P, %VAL( SINV ),
     $                LDSINV, %VAL( AINV ), LDAINV, %VAL( B ), LDB,
     $                %VAL( RINV ), LDRINV, %VAL( C ), LDC,
     $                %VAL( QINV ), LDQINV, %VAL( X ), %VAL( RINVY ),
     $                %VAL( Z ), %VAL( E ), TOL, IWORK, %VAL( DWORK ),
     $                LDWORK, INFO )
      ELSE IF ( ITASK.EQ.-1 ) THEN
         CALL FB01QD( JOBK, MULTBQ, N, M, P, %VAL( S ), LDS, %VAL( A ),
     $                LDA, %VAL( B ), LDB, %VAL( Q ), LDQ, %VAL( C ),
     $                LDC, %VAL( R ), LDR, %VAL( K ), LDK, TOL, IWORK,
     $                %VAL( DWORK ), LDWORK, INFO )
      ELSE IF ( ITASK.EQ.0 ) THEN
         CALL FB01VD( N, M, P, %VAL( W ), LDW, %VAL( A ), LDA,
     $                %VAL( B ), LDB, %VAL( C ), LDC, %VAL( Q ), LDQ,
     $                %VAL( R ), LDR, %VAL( K ), LDK, TOL, IWORK,
     $                %VAL( DWORK ), LDWORK, INFO )
         IRCND = 1
      ELSE IF ( ITASK.EQ.1 ) THEN
         CALL FB01RD( JOBK, MULTBQ, N, M, P, %VAL( S ), LDS, %VAL( A ),
     $                LDA, %VAL( B ), LDB, %VAL( Q ), LDQ, %VAL( C ),
     $                LDC, %VAL( R ), LDR, %VAL( K ), LDK, TOL, IWORK,
     $                %VAL( DWORK ), LDWORK, INFO )
      ELSE IF ( ITASK.EQ.2 ) THEN
         CALL FB01TD( JOBX, MULTRC, N, M, P, %VAL( SINV ), LDSINV,
     $                %VAL( AINV ), LDAINV, %VAL( AINVB ), LDAINB,
     $                %VAL( RINV ), LDRINV, %VAL( C ), LDC,
     $                %VAL( QINV ), LDQINV, %VAL( X ), %VAL( RINVY ),
     $                %VAL( Z ), %VAL( E ), TOL, IWORK, %VAL( DWORK ),
     $                LDWORK, INFO )
      ELSE
         CALL FD01AD( JP, L, LAMBDA, XIN, YIN, EFOR, %VAL( XF ),
     $                %VAL( EPSBCK ), %VAL( CTETA ), %VAL( STETA ),
     $                %VAL( YQ ), EPOS, EOUT, %VAL( SALPH ), IWARN,
     $                INFO )
      END IF
C
      IF ( ITASK.NE.3 ) THEN
         IF ( ITASK.EQ.0 .OR. ( IATASK.EQ.1 .AND. WITHK ) .OR.
     $                        ( IATASK.EQ.2 .AND. WITHX ) ) THEN
            CALL MA02KV( 1, IRCND, 1, %VAL( DWORK ), OP )
            RCND = OP(1)
         ELSE
            RCND = ONE
         END IF
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         IF ( ITASK.LE.2 ) THEN
            IF ( IATASK.EQ.2 ) THEN
               PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
               CALL mxCopyReal8ToPtr( %VAL( SINV ),
     $                                mxGetPr( PLHS(1) ), N*N )
               IP = 2
               IF ( NLHS.GE.IP ) THEN
                  PLHS(IP) = mxCreateDoubleMatrix( N, 1, 0 )
                  CALL mxCopyReal8ToPtr( %VAL( X ),
     $                                   mxGetPr( PLHS(IP) ), N )
                  IP = IP + 1
               END IF
               IF ( NLHS.GE.IP ) THEN
                  PLHS(IP) = mxCreateDoubleMatrix( M, M, 0 )
                  CALL mxCopyReal8ToPtr( %VAL( QINV ),
     $                                   mxGetPr( PLHS(IP) ), M*M )
                  IP = IP + 1
                  IF ( NLHS.GE.IP ) THEN
                     PLHS(IP) = mxCreateDoubleMatrix( P, 1, 0 )
                     CALL mxCopyReal8ToPtr( %VAL( E ),
     $                                      mxGetPr( PLHS(IP) ), P )
                     IP = IP + 1
                  END IF
               END IF
            ELSE
               IF ( IATASK.EQ.1 ) THEN
                  PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
                  CALL mxCopyReal8ToPtr( %VAL( S ),
     $                                   mxGetPr( PLHS(1) ), N*N )
               ELSE
                  PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
                  CALL mxCopyReal8ToPtr( %VAL( W ),
     $                                   mxGetPr( PLHS(1) ), N*N )
               END IF
               IP = 2
               IF ( NLHS.GE.IP ) THEN
                  PLHS(IP) = mxCreateDoubleMatrix( N, P, 0 )
                  CALL mxCopyReal8ToPtr( %VAL( K ),
     $                                   mxGetPr( PLHS(IP) ), N*P )
                  IP = IP + 1
               END IF
               IF ( NLHS.GE.IP ) THEN
                  PLHS(IP) = mxCreateDoubleMatrix( P, P, 0 )
                  CALL mxCopyReal8ToPtr( %VAL( R ),
     $                                   mxGetPr( PLHS(IP) ), P*P )
                  IP = IP + 1
               END IF
            END IF
            IF ( NLHS.GE.IP ) THEN
               PLHS(IP) = mxCreateDoubleMatrix( 1, 1, 0 )
               CALL mxCopyReal8ToPtr( RCND, mxGetPr( PLHS(IP) ),
     $                                1 )
            END IF
C
         ELSE
C
            PLHS(1) = mxCreateDoubleMatrix( 1, 1, 0 )
            CALL mxCopyReal8ToPtr( EFOR, mxGetPr( PLHS(1) ), 1 )
            IF ( NLHS.GE.2 ) THEN
               PLHS(2) = mxCreateDoubleMatrix( L, 1, 0 )
               CALL mxCopyReal8ToPtr( %VAL( XF ),
     $                                mxGetPr( PLHS(2) ), L )
               IF ( NLHS.GE.3 ) THEN
                  PLHS(3) = mxCreateDoubleMatrix( L+1, 1, 0 )
                  CALL mxCopyReal8ToPtr( %VAL( EPSBCK ),
     $                                   mxGetPr( PLHS(3) ), L+1 )
                  IF ( NLHS.GE.4 ) THEN
                     PLHS(4) = mxCreateDoubleMatrix( L, 1, 0 )
                     CALL mxCopyReal8ToPtr( %VAL( CTETA ),
     $                                      mxGetPr( PLHS(4) ), L )
                     IF ( NLHS.GE.5 ) THEN
                        PLHS(5) = mxCreateDoubleMatrix( L, 1, 0 )
                        CALL mxCopyReal8ToPtr( %VAL( STETA ),
     $                                         mxGetPr( PLHS(5) ),
     $                                         L )
                        IP = 6
                        OP(1) = EPOS
                        IF ( WITHF ) THEN
                           OP(2) = EOUT
                           ISIZE = 2
                           IF ( NLHS.GE.6 ) THEN
                              PLHS(6) = mxCreateDoubleMatrix( L, 1,
     $                                                        0 )
                              CALL mxCopyReal8ToPtr( %VAL( YQ ),
     $                                           mxGetPr( PLHS(6) ),
     $                                               L )
                              IP = IP + 1
                           END IF
                        ELSE
                           ISIZE = 1
                        END IF
                        IF ( NLHS.GE.IP ) THEN
                           PLHS(IP) = mxCreateDoubleMatrix( ISIZE,
     $                                                      1, 0 )
                           CALL mxCopyReal8ToPtr( OP,
     $                                          mxGetPr( PLHS(IP) ),
     $                                            ISIZE )
                           IP = IP + 1
                           IF ( NLHS.GE.IP ) THEN
                              PLHS(IP) = mxCreateDoubleMatrix( L, 1,
     $                                                         0 )
                              CALL mxCopyReal8ToPtr( %VAL( SALPH ),
     $                                           mxGetPr( PLHS(IP)),
     $                                               L )
                           END IF
                        END IF
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      IF ( ITASK.LE.2  ) THEN
         DEALLOCATE ( IWORK )
         CALL mxFree( C     )
         CALL mxFree( DWORK )
         IF ( ITASK.LE.1  ) THEN
            CALL mxFree( B )
            IF ( IATASK.LE.1  ) THEN
               CALL mxFree( A )
               CALL mxFree( K )
               CALL mxFree( Q )
               CALL mxFree( R )
               IF ( IATASK.EQ.1  ) THEN
                  CALL mxFree( S )
               ELSE
                  CALL mxFree( W )
               END IF
            END IF
         END IF
         IF ( IATASK.EQ.2 ) THEN
            CALL mxFree( AINV )
            CALL mxFree( E )
            CALL mxFree( QINV )
            CALL mxFree( RINV )
            CALL mxFree( RINVY )
            CALL mxFree( SINV )
            CALL mxFree( X )
            CALL mxFree( Z )
            IF ( ITASK.EQ.2 )
     $         CALL mxFree( AINVB )
         END IF
      ELSE
         CALL mxFree( CTETA )
         CALL mxFree( EPSBCK )
         CALL mxFree( SALPH )
         CALL mxFree( STETA )
         CALL mxFree( XF )
         IF ( WITHF )
     $      CALL mxFree( YQ )
      END IF
C
C Error and warning handling.
C
      IF ( IWARN.NE.0 .AND. PRINTW ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM FD01AD'')'
     $        ) IWARN
         CALL mexPrintf( TEXT )
      END IF
C
      IF ( INFO.NE.0 ) THEN
         IF ( ITASK.EQ.-2 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM FB01SD'')'
     $           ) INFO
         ELSE IF ( ITASK.EQ.-1 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM FB01QD'')'
     $           ) INFO
         ELSE IF ( ITASK.EQ.0 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM FB01VD'')'
     $           ) INFO
         ELSE IF ( ITASK.EQ.1 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM FB01RD'')'
     $           ) INFO
         ELSE IF ( ITASK.EQ.2 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM FB01TD'')'
     $           ) INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of KFILTUPD ***
      END
