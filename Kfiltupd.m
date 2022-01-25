% KFILTUPD.F - MEX-function for computing a combined measurement and
%              time update of various Kalman filters for discrete-time
%              systems, using SLICOT routines FB01QD, FB01RD, FB01SD,
%              FB01TD, FB01VD, and FD01AD.
%
%   [(S1,K,Rin)(Sinv1,X1,Qin,E)(P1,K,Rin)(,rcnd)
%    (Ef1,Xf1,Eb1,Cs1,Sn1,YQ1,Err,Salph)] =
%          Kfiltupd(task(,opt)(,S,A,B,C(,Q),R)
%                             (,Sinv,Ainv,B,C,Qinv(,Rinv),X,Rinvy,Z)
%                             (,Sinv,Ainv,AinvB,C,Qinv(,Rinv),X,Rinvy,Z)
%                             (,P,A,B,C,Q,R)
%                             (,L,Lambda,Si,Ef,Xf,Eb,Cs,Sn(,YQ))
%                             (,tol))
%
%   [S1(,K,Rin,rcnd)] = Kfiltupd(-1/1,[WithK,BxQ],S,A,B,C(,Q),R(,tol))
%   [Sinv1(,X1,Qin,E,rcnd)] =
%                       Kfiltupd(-2,[WithX,AixB,RixC],Sinv,Ainv,B,C,Qinv
%                                (,Rinv),X,Rinvy,Z(,tol))
%   [Sinv1(,X1,Qin,E,rcnd)] =
%                       Kfiltupd(2,[WithX,RixC],Sinv,Ainv,AinvB,C,Qinv
%                                (,Rinv),X,Rinvy,Z(,tol))
%   [P1(,K,Rin,rcnd)] = Kfiltupd(0,P,A,B,C,Q,R(,tol))
%   [Ef1,Xf1,Eb1,Cs1,Sn1(,YQ1),Err,Salph] =
%                       Kfiltupd(3,[WithF,printw],L,Lambda,Si,Ef,Xf,Eb,
%                                Cs,Sn(,YQ))
%
%   KFILTUPD computes a combined measurement and time update of one
%   iteration of the time-varying or time-invariant Kalman filter in
%   the square root covariance or information form, or of the
%   conventional Kalman filter. One time update of the recursive
%   least-squares filtering problem can also be computed, using a fast
%   QR-decomposition based approach.
%
%   Description of input parameters:
%   task   - scalar indicating the calculation to be performed:
%            =-2 :  calculate a combined measurement and time update of
%                   one iteration of the time-varying Kalman filter in
%                   the square root information form;
%            =-1 :  calculate a combined measurement and time update of
%                   one iteration of the time-varying Kalman filter in
%                   the square root covariance form;
%            = 0 :  calculate one recursion of the conventional Kalman
%                   filter equations;
%            = 1 :  calculate a combined measurement and time update of
%                   one iteration of the time-invariant Kalman filter in
%                   the square root covariance form;
%            = 2 :  calculate a combined measurement and time update of
%                   one iteration of the time-invariant Kalman filter in
%                   the square root information form;
%            = 3 :  calculate one time update of the recursive least-
%                   squares filtering problem.
%   opt    - (optional) integer scalar containing option values.
%          - If |task| = 1, opt has length 2, opt = [WithK,BxQ].
%            WithK indicates whether the Kalman filter gain matrix K
%            should be computed, as follows:
%            = 0 :  K is not required;
%            = 1 :  K is computed and stored in array K.
%            Default: WithK = 0.
%            BxQ indicates how matrices B and Q are to be passed to the
%            gateway, as follows:
%            = 0 :  arrays B and Q must contain the matrices as
%                   described below;
%            = 1 :  B must contain the product B*sqrt(Q) and Q is not
%                   used.
%            Default: BxQ = 0.
%          - If task = -2, opt has length 3, opt = [WithX,AixB,RixC];
%            if task =  2, opt has length 2, opt = [WithX,RixC].
%            WithX indicates whether X    is to be computed, as follows:
%                                     i+1
%            = 0 :  X is not required;
%            = 1 :  X is computed and stored in array X.
%            Default: WithX = 0.          -1
%            AixB indicates how matrices A   and B  are to be passed to
%                                         i       i
%            the gateway, as follows:
%            = 0 :  arrays Ainv and B must contain the matrices as
%                   described below;                    -1
%            = 1 :  array Ainv must contain the matrix A   and the array
%                                               -1      i
%                   B must contain the product A  B .
%                                               i  i
%            Default: AixB = 0.           -1/2
%            RixC indicates how matrices R     and C    are to be passed
%                                         i+1       i+1
%            to the gateway, as follows:
%            = 0 :  arrays Rinv and C must contain the matrices as
%                   described below;                         -1/2
%            = 1 :  array C must contain the matrix product R    C
%                                                            i+1  i+1
%
%            and array Rinv is not used.
%            Default: RixC = 0.
%          - If task = 3, opt has length 2, opt = [WithF,printw].
%            WithF indicates whether both prediction and filtering parts
%            are to be applied, as follows:
%            = 0 :  only the prediction section is to be applied;
%            = 1 :  both prediction and filtering parts are to be
%                   applied.
%            Default: WithF = 0.
%            printw is a switch for printing the warning messages.
%            = 0:  do not print warning messages;
%            = 1:  print warning messages.
%            Default: printw = 0.
%          - If task = 0, opt must not be specified as an input
%            parameter.
%   S      - the n-by-n lower triangular part contains S   , the square
%                                                       i-1
%            root (left Cholesky factor) of the state covariance matrix
%            at instant (i-1). The strict upper triangular part is not
%            used.
%   A      - the n-by-n state transition matrix.
%          - If task = -1 or task = 0, it contains A , the state
%                                                   i
%            transition matrix at instant i.
%          - If task = 1, it contains A, the state transition matrix of
%            the system in lower observer Hessenberg form (e.g., as
%            produced by SLICOT Library Routine TB01ND).
%   B      - the n-by-m input weight matrix B  at instant i.
%                                            i                1/2
%          - If |task| = 1, it contains B  (or the product B Q   ,  if
%                                        i                  i i
%            BxQ = opt(2) = 1).                            -1
%          - If task = -2, it contains B  (or the product A  B ,  if
%                                       i                  i  i
%            AixB = opt(2) = 1).
%   C      - if task = -1 or task = 0, the p-by-n output weight
%            matrix C , at instant i.
%                    i
%          - If task = 1, the p-by-n output weight matrix C of the
%            system in lower observer Hessenberg form (e.g., as produced
%            by SLICOT Library routine TB01ND).
%          - If |task| = 2, the p-by-n output weight matrix C   , (or
%                        -1/2                                i+1
%            the product R    C    if RixC = 1) at instant i+1.
%                         i+1  i+1
%   Q      - if |task| = 1 and BxQ = 0, the m-by-m lower triangular
%                           1/2
%            part contains Q   , the square root (left Cholesky factor)
%                           i
%            of the input (process) noise covariance matrix at
%            instant i.
%            The strict upper triangular part is not used.
%          - If |task| = 1 and BxQ = 1, Q must not be specified as an
%            input parameter.
%          - If task = 0, the m-by-m matrix Q , the process noise
%                                            i
%            covariance matrix at instant i.
%   R      - if |task| = 1, the p-by-p lower triangular part contains
%             1/2
%            R   , the square root (left Cholesky factor) of the output
%             i
%            (measurement) noise covariance matrix at instant i.
%            The strict upper triangular part of this array is not used.
%          - If task = 0, the p-by-p output noise covariance matrix
%            R , at instant i.
%             i
%   tol    - (optional) real scalar containing a tolerance value used to
%                                                           1/2
%            test for near singularity of a matrix: (RINOV )   , if
%                                                -1       i
%            |task| = 1 and WithK = opt(1) = 1; S   , if |task| = 2 and
%                                                i+1
%            WithX = opt(1) = 1; and (RINOV) , if task = 0.
%                                           i
%            If tol > 0, then the given value of tol is used as a
%            lower bound for the reciprocal condition number of that
%            matrix; a matrix whose estimated condition number is less
%            than 1/tol is considered to be nonsingular. If tol <= 0,
%            then ps*eps, is used instead, where ps is the product of
%            the matrix dimensions, and eps is the machine precision.
%          - If task = 3, tol must not be specified as an input
%            parameter.
%            Default: tol = 0.
%                                                       -1
%   Sinv   - the n-by-n upper triangular part contains S  , the inverse
%                                                       i
%            of the square root (right Cholesky factor) of the state
%            covariance matrix P    (hence the information square root)
%                               i|i
%            at instant i.
%            The strict lower triangular part of this array is not used.
%                                             -1
%   Ainv   - if task = -2, the n-by-n matrix A  , the inverse of the
%                                             i
%            state transition matrix of the system at instant i.
%            The matrix A is time-varying if task = -2.
%                                            -1
%          - If task = 2, the n-by-n matrix A  , the inverse of the
%            state transition matrix of the system in controller
%            Hessenberg form (e.g., as produced by SLICOT Library
%            Routine TB01MD).
%                               -1                   -1
%   AinvB  - the n-by-m matrix A  B, the product of A   and the input
%            weight matrix B of the system, in upper controller
%            Hessenberg form (e.g., as produced by SLICOT Library
%            Routine TB01MD).
%                                                       -1/2
%   Qinv   - the m-by-m upper triangular part contains Q    , the
%                                                       i
%            inverse of the covariance square root (right Cholesky
%            factor) of the process noise (hence the information square
%            root) at instant i.
%            The strict lower triangular part of this array is not used.
%   Rinv   - if RixC = 0, the p-by-p upper triangular part contains
%             -1/2
%            R    , the inverse of the covariance square root (right
%             i+1
%            Cholesky factor) of the output noise (hence the information
%            square root) at instant i+1.
%            The strict lower triangular part of this array is not used.
%            Otherwise, Rinv must not be specified as an input
%            parameter.
%   X      - the n-vector X , the estimated filtered state at instant i.
%                          i
%                          -1/2
%   Rinvy  - the p-vector R    Y   , the product of the upper triangular
%                          i+1  i+1
%                    -1/2
%            matrix R     and the measured output vector Y    at instant
%                    i+1                                  i+1
%            i+1.
%   Z      - the m-vector Z , the mean value of the state process noise
%                          i
%            at instant i.
%   P      - the n-by-n state covariance matrix P      at instant (i-1).
%                                                i|i-1
%            The upper triangular part only is needed.
%   L      - integer scalar specifying the length of the impulse
%            response of the equivalent transversal filter model.
%            L >= 1.
%   Lambda - real scalar specifying the square root of the forgetting
%            factor. For tracking capabilities and exponentially stable
%            error propagation, Lambda < 1.0 (strict inequality) should
%            be used.  0.0 < Lambda <= 1.0.
%   Si     - Si(1) = Ui contains the input sample at instant i.
%            If WithF = 1, Si(2) = Yi contains the reference sample at
%            instant i. Otherwise, Yi is not used.
%            (The situation just before and just after the call of the
%            gateway are denoted by instant (i-1) and instant i,
%            respectively.)
%   Ef     - the square root of exponentially weighted forward
%            prediction error energy at instant (i-1).  Ef >= 0.0.
%   Xf     - L-vector containing the transformed forward prediction
%            variables at instant (i-1).
%   Eb     - (L+1)-vector; the leading L elements must contain the
%            normalized a posteriori backward prediction error residuals
%            of orders zero through L-1, respectively, at instant (i-1),
%            and Eb(L+1) must contain the square-root of the so-called
%            "conversion factor" at instant (i-1).
%   Cs     - L-vector containing the cosines of the rotation angles used
%            in time updates, at instant (i-1).
%   Sn     - L-vector containing the sines of the rotation angles used
%            in time updates, at instant (i-1).
%   YQ     - if WithF = 1, L-vector containing the orthogonally
%            transformed reference vector at instant (i-1). These
%            elements are also the tap multipliers of an equivalent
%            normalized lattice least-squares filter.
%            Otherwise, this parameter is not used.
%
%   Description of output parameters:
%   S1     - the n-by-n lower triangular part contains S , the square
%                                                       i
%            root (left Cholesky factor) of the state covariance matrix
%            at instant i.
%   K      - if |task| = 1 and WithK = 1, or if task = 0, the n-by-p
%            Kalman filter gain matrix K  at instant i. If |task| = 1
%                                       i
%            and WithK = 0, the n-by-p matrix AK , the product of the
%                                               i
%            input matrix A and the Kalman filter gain matrix at
%            instant i.
%   Rin    - if |task| = 1, the p-by-p lower triangular part contains
%                    1/2
%            (RINOV )   , the square root (left Cholesky factor) of the
%                  i
%            covariance matrix of the innovations at instant i.
%          - If task = 0, the p-by-p upper triangular part contains
%                    1/2
%            (RINOV )   , the square root (left Cholesky factor) of the
%                  i
%            covariance matrix of the innovations at instant i.
%   rcnd   - (optional) if |task| = 1 and WithK = 1, an estimate of the
%            reciprocal of the condition number (in the 1-norm) of
%                    1/2
%            (RINOV )   .
%                  i
%          - If task = 0, an estimate of the reciprocal condition number
%            (in the 1-norm) of RINOV .
%                                    i
%          - If |task| = 2 and WithX = 1, an estimate of the reciprocal
%                                                        -1
%            of the condition number (in the 1-norm) of S   .
%                                                        i+1
%            If |task| = 1 and WithK = 0, or |task| = 2 and WithX = 0, 
%            rcnd is not computed, but set to 1.
%                                                       -1
%   Sinv1  - the n-by-n upper triangular part contains S   , the inverse
%                                                       i+1
%            of the square root (right Cholesky factor) of the state
%            covariance matrix P        (hence the information square
%                               i+1|i+1
%            root) at instant i+1.
%   X1     - if WithX = 1, the n-vector X   , the estimated filtered
%                                        i+1
%            state at instant i+1.
%                                        -1
%          - If WithX = 0, the n-vector S   X   .
%                                        i+1 i+1
%                                                              -1/2
%   Qin    - the m-by-m upper triangular part contains (QINOV )    , the
%                                                            i
%            inverse of the covariance square root (right Cholesky
%            factor) of the process noise (hence the information square
%            root) at instant i.
%   E      - the p-vector E   , the estimated error at instant i+1.
%                          i+1
%   P1     - the n-by-n upper triangular part contains the upper
%            triangular part of the state covariance matrix P     ,
%                                                            i+1|i
%            at instant i.
%   Ef1    - the square root of the exponentially weighted forward
%            prediction error energy at instant i.
%   Xf1    - L-vector containing the transformed forward prediction
%            variables at instant i.
%   Eb1    - (L+1)-vector containing the normalized a posteriori
%            backward prediction error residuals, plus the square root
%            of the conversion factor at instant i.
%   Cs1    - L-vector containing the cosines of the rotation angles at
%            instant i.
%   Sn1    - L-vector containing the sines of the rotation angles at
%            instant i.
%   YQ1    - if WithF = 1, L-vector containing the orthogonally
%            transformed reference vector at instant i.
%            If WithF = 0, YQ1 must not be specified as an output
%            parameter.
%   Err    - Err(1) = Ep contains the a posteriori forward prediction
%            error residual.
%            If WithF = 1, Err(2) = Eo contains the a posteriori output
%            error residual from the least-squares filter at instant i.
%   Salph  - L-vector: the element Salph(i), i=1,...,L, contains the
%            opposite of the i-(th) reflection coefficient for the
%            least-squares normalized lattice predictor (whose value is
%            -Salph(i)).
%
%   Method:
%   If task = -2, the gateway performs one recursion of the square root
%   information filter algorithm, summarized as follows:
%
%     |    -1/2             -1/2    |     |         -1/2             |
%     |   Q         0      Q    Z   |     | (QINOV )     *     *     |
%     |    i                i    i  |     |       i                  |
%     |                             |     |                          |
%     |  -1 -1     -1 -1    -1      |     |             -1    -1     |
%   T | S  A  B   S  A     S  X     |  =  |    0       S     S   X   |
%     |  i  i  i   i  i     i  i    |     |             i+1   i+1 i+1|
%     |                             |     |                          |
%     |           -1/2      -1/2    |     |                          |
%     |    0     R    C    R    Y   |     |    0         0     E     |
%     |           i+1  i+1  i+1  i+1|     |                     i+1  |
%
%                (Pre-array)                      (Post-array)
%
%   where T is an orthogonal transformation triangularizing the
%                      -1/2
%   pre-array, (QINOV )     is the inverse of the covariance square
%                    i
%   root (right Cholesky factor) of the process noise (hence the
%   information square root) at instant i, and E    is the estimated
%                                               i+1
%   error at instant i+1.
%
%   The inverse of the corresponding state covariance matrix P
%                                                             i+1|i+1
%   (hence the information matrix I) is then factorized as
%
%                 -1         -1     -1
%      I       = P       = (S   )' S
%       i+1|i+1   i+1|i+1    i+1   i+1
%
%   and one combined time and measurement update for the state is
%   given by X   .
%             i+1
%
%   The triangularization is done entirely via Householder
%   transformations exploiting the zero pattern of the pre-array.
%
%   If task = 2, the gateway performs one recursion of the square root
%   information filter algorithm, summarized as follows:
%
%     |    -1/2             -1/2    |     |         -1/2             |
%     |   Q         0      Q    Z   |     | (QINOV )     *     *     |
%     |    i                i    i  |     |       i                  |
%     |                             |     |                          |
%     |           -1/2      -1/2    |     |             -1    -1     |
%   T |    0     R    C    R    Y   |  =  |    0       S     S   X   |
%     |           i+1  i+1  i+1  i+1|     |             i+1   i+1 i+1|
%     |                             |     |                          |
%     |  -1 -1     -1 -1    -1      |     |                          |
%     | S  A  B   S  A     S  X     |     |    0         0     E     |
%     |  i         i        i  i    |     |                     i+1  |
%
%                 (Pre-array)                      (Post-array)
%
%   where T is an orthogonal transformation triangularizing the
%                      -1/2
%   pre-array, (QINOV )     is the inverse of the covariance square
%                    i
%   root (right Cholesky factor) of the process noise (hence the
%                                               -1  -1
%   information square root) at instant i and (A  ,A  B) is in upper
%   controller Hessenberg form.
%
%   An example of the pre-array is given below (where n = 6, m = 2,
%   and p = 3):
%
%       |x x |             | x|
%       |  x |             | x|
%       _______________________
%       |    | x x x x x x | x|
%       |    | x x x x x x | x|
%       |    | x x x x x x | x|
%       _______________________ .
%       |x x | x x x x x x | x|
%       |  x | x x x x x x | x|
%       |    | x x x x x x | x|
%       |    |   x x x x x | x|
%       |    |     x x x x | x|
%       |    |       x x x | x|
%
%   The inverse of the corresponding state covariance matrix P
%                                                             i+1|i+1
%   (hence the information matrix I) is then factorized as
%
%                  -1         -1    -1
%       I       = P       = (S   )' S
%        i+1|i+1   i+1|i+1    i+1   i+1
%
%   and one combined time and measurement update for the state is
%   given by X   .
%             i+1
%
%   The triangularization is done entirely via Householder
%   transformations exploiting the zero pattern of the pre-array.
%
%   If task = -1, the gateway performs one recursion of the square root
%   covariance filter algorithm, summarized as follows:
%
%    |  1/2                      |     |         1/2          |
%    | R      C x S      0       |     | (RINOV )     0     0 |
%    |  i      i   i-1           |     |       i              |
%    |                      1/2  | T = |                      |
%    | 0      A x S    B x Q     |     |     AK       S     0 |
%    |         i   i-1  i   i    |     |       i       i      |
%
%        (Pre-array)                      (Post-array)
%
%   where T is an orthogonal transformation triangularizing the
%   pre-array.
%
%   The state covariance matrix P    is factorized as
%                                i|i-1
%      P     = S  S'
%       i|i-1   i  i
%
%   and one combined time and measurement update for the state X
%                                                               i|i-1
%   is given by
%
%      X     = A X      + K (Y - C X     ),
%       i+1|i   i i|i-1    i  i   i i|i-1
%
%                        -1/2
%   where K = AK (RINOV )     is the Kalman filter gain matrix and Y
%          i    i      i                                            i
%   is the observed output of the system.
%
%   The triangularization is done entirely via Householder
%   transformations exploiting the zero pattern of the pre-array.
%
%   If task = 1, the gateway performs one recursion of the square root
%   covariance filter algorithm, summarized as follows:
%
%    |  1/2                     |     |         1/2          |
%    | R      0        C x S    |     | (RINOV )     0     0 |
%    |  i                   i-1 |     |       i              |
%    |             1/2          | T = |                      |
%    | 0      B x Q    A x S    |     |     AK       S     0 |
%    |         i   i        i-1 |     |       i       i      |
%
%         (Pre-array)                      (Post-array)
%
%   where T is unitary and (A,C) is in lower observer Hessenberg form.
%
%   An example of the pre-array is given below (where n = 6, p = 2
%   and m = 3):
%
%        |x   |      | x          |
%        |x x |      | x x        |
%        |____|______|____________|
%        |    | x x x| x x x      |
%        |    | x x x| x x x x    | .
%        |    | x x x| x x x x x  |
%        |    | x x x| x x x x x x|
%        |    | x x x| x x x x x x|
%        |    | x x x| x x x x x x|
%
%   The corresponding state covariance matrix P      is then
%                                              i|i-1
%   factorized as
%
%       P     = S  S'
%        i|i-1   i  i
%
%   and one combined time and measurement update for the state X
%                                                               i|i-1
%   is given by
%
%       X     = A X      + K (Y - C X     )
%        i+1|i     i|i-1    i  i     i|i-1
%
%                        -1/2
%   where K = AK (RINOV )     is the Kalman filter gain matrix and Y
%          i    i      i                                            i
%   is the observed output of the system.
%
%   The triangularization is done entirely via Householder
%   transformations exploiting the zero pattern of the pre-array.
%
%   If task = 0, the gateway performs one recursion of the conventional
%   Kalman filter. The Kalman filter gain used at the i-th recursion
%   step is of the form
%
%                          -1
%      K  = P     C'  RINOV  ,
%       i    i|i-1 i       i
%
%   where RINOV  = C P     C' + R , and the state covariance matrix
%              i    i i|i-1 i    i
%
%   P      is updated by the discrete-time Riccati equation
%    i|i-1
%
%      P      = A  (P      - K C P     ) A'  + B Q B'.
%       i+1|i    i   i|i-1    i i i|i-1   i     i i i
%
%   Using these two updates, the combined time and measurement update
%   of the state X      is given by
%                 i|i-1
%
%      X      = A X      + A K (Y  - C X     ),
%       i+1|i    i i|i-1    i i  i    i i|i-1
%
%   where Y  is the new observation at step i.
%          i
%
%   If task = 3, the gateway performs one recursion of the recursive
%   least-squares filter. The output error Eo at instant i, denoted by
%   Eo(i), is the reference sample minus a linear combination of L
%   successive input samples:
%
%                      L-1
%      Eo(i) = Yi(i) - SUM h_j * Ui(j-i),
%                      j=0
%
%   where Yi(i) and Ui(i) are the scalar samples at instant i.
%   A least-squares filter uses those h_0,...,h_{L-1} which minimize
%   an exponentially weighted sum of successive output errors squared:
%
%       i
%      SUM [Lambda**(2(i-k)) * Eo(k)**2].
%      k=1
%
%   Each gateway call performs a time update of the least-squares
%   filter using a fast least-squares algorithm derived from a
%   QR decomposition. The algorithm does not compute the parameters
%   h_0,...,h_{L-1} from the above formula, but instead delivers the
%   parameters of an equivalent normalized least-squares lattice filter,
%   which are available from the arrays Salph (reflection coefficients)
%   and YQ (tap multipliers), as well as the exponentially weighted
%   input signal energy
%
%       i                                           L
%      SUM [Lambda**(2(i-k)) * Ui(k)**2] = Ef**2 + SUM Xf(j)**2.
%      k=1                                         j=1
%
%   Comments for task = 3:
%   1.  For tracking capabilities and exponentially stable error
%       propagation, Lambda < 1.0 should be used.  Lambda is typically
%       chosen slightly less than 1.0 so that "past" data are 
%       exponentially forgotten.
%   2.  Prior to the first gateway call, the variables must be
%       initialized. The following initial values are recommended:
%
%       Xf(i)   = 0.0,        i=1,...,L
%       Eb(i)   = 0.0,        i=1,...,L
%       Eb(L+1) = 1.0
%       Cs(i)   = 1.0,        i=1,...,L
%       Sn(i)   = 0.0,        i=1,...,L
%       YQ(i)   = 0.0,        i=1,...,L
%
%       Ef      = 0.0         (exact start)
%       Ef      = "small positive constant" (soft start).
%
%       Soft starts are numerically more reliable, but result in a
%       biased least-squares solution during the first few iterations.
%       This bias decays exponentially fast provided Lambda < 1.0.
%       If sigma is the standard deviation of the input sequence
%       Xi, then initializing Ef = sigma*1.0E-02 usually works well.
%   3.  The algorithm is backward consistent for all input sequences Xi,
%       and backward stable for persistently exciting input sequences,
%       assuming Lambda < 1.0.
%       If the condition of the signal is very poor (IWARN = 1), then
%       the results are not guaranteed to be reliable.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2003.
%
% Revisions:
%   -
%
