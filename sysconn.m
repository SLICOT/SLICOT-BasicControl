% SYSCONN.F  - MEX-function to compute a state-space model
%              (A,B,C,D) for various inter-connections of two systems,
%              using SLICOT routines AB05MD, AB05ND, AB05OD, AB05PD,
%              and AB05QD.
%
%   [A,B,C,D] = SYSCONN(task,A1,B1,C1,D1,A2,B2,C2,D2(,uplo)(,alpha))
%
%   [A,B,C,D] = SYSCONN(1,A1,B1,C1,D1,A2,B2,C2,D2(,uplo))
%   [A,B,C,D] = SYSCONN(2,A1,B1,C1,D1,A2,B2,C2,D2(,alpha))
%   [A,B,C,D] = SYSCONN(3,A1,B1,C1,D1,A2,B2,C2,D2(,alpha))
%   [A,B,C,D] = SYSCONN(4,A1,B1,C1,D1,A2,B2,C2,D2(,alpha))
%   [A,B,C,D] = SYSCONN(5,A1,B1,C1,D1,A2,B2,C2,D2)
%
%   SYSCONN computes a state-space model (A,B,C,D) for various inter-
%   connections (depending on the value of task) of two systems,
%   given in state-space form, (A1,B1,C1,D1), and (A2,B2,C2,D2).
%   Specifically,
%
%   task = 1: compute cascaded inter-connection;
%   task = 2: compute feedback inter-connection;
%   task = 3: compute rowwise concatenation (parallel inter-connection
%             on outputs, with separate inputs);
%   task = 4: compute the state-space model (A,B,C,D) corresponding to
%             the sum G = G1 + alpha*G2, where G, G1, and G2 are the 
%             transfer-function matrices of the corresponding state-
%             space models (A,B,C,D), (A1,B1,C1,D1), and (A2,B2,C2,D2),
%             respectively;
%   task = 5: append two systems in state-space form, (A1,B1,C1,D1) and
%             (A2,B2,C2,D2), with the transfer-function matrices G1
%             and G2, respectively, and obtain the state-space model
%             (A,B,C,D) corresponding to the transfer-function matrix
%
%                           ( G1 0  )
%                       G = (       ) .                              (1)
%                           ( 0  G2 )
%
%   Description of input parameters: 
%   task   - integer specifying the computations to be performed.
%            task = 1 :  compute the cascaded inter-connection;
%            task = 2 :  compute the feedback inter-connection;
%            task = 3 :  compute the rowwise  inter-connection;
%            task = 4 :  compute the parallel inter-connection;
%            task = 5 :  compute the compound system in (1).
%   A1     - the n1-by-n1 state dynamics matrix A1.
%   B1     - the n1-by-m1 input/state matrix B1.
%   C1     - the p1-by-n1 state/output matrix C1.
%   D1     - the p1-by-m1 input/output matrix D1.
%   A2     - the n2-by-n2 state dynamics matrix A2.
%   B2     - the n2-by-m2 input/state matrix B2.  If task = 1, m2 = p1.
%                                                 If task = 2, m2 = p1.
%                                                 If task = 4, m2 = m1.
%   C2     - the p2-by-n2 state/output matrix C2. If task = 2, p2 = m1.
%                                                 If task = 3, p2 = p1.
%                                                 If task = 4, p2 = p1.
%   D2     - the p2-by-m2 input/output matrix D2.
%   uplo   - (optional) integer indicating whether the matrix A should
%            be obtained in the upper or lower block diagonal form:
%            uplo = 1 :  Obtain A in the lower block diagonal form;
%            uplo = 2 :  Obtain A in the upper block diagonal form.
%            Default:  uplo = 1.
%   alpha  - (optional) real coefficient multiplying the transfer-
%            function matrix (or the output equation) of the second
%            system. For task = 2,
%            alpha = +1 corresponds to positive feedback, and
%            alpha = -1 corresponds to negative feedback.
%            For task = 3 or task = 4, alpha is not constrained, but
%            alpha = 0 is not dealt with as a special case.
%            Default:  alpha = -1, for task = 2;
%                      alpha =  1, for task = 3 or task = 4.
%
%   Description of output parameters: 
%   A      - the n-by-n state dynamics matrix A of the obtained system,
%            where n = n1+n2.
%   B      - the n-by-m input/state matrix B, where m = m1, if task <= 2
%            or task = 4, and m = m1 + m2, if task = 3 or task = 5.
%   C      - the p-by-n state/output matrix C, where p = p2 if task = 1,
%            p = p1 if task = 2, 3, or 4, and p = p1 + p2, if task = 5.
%   D      - the p-by-m input/output matrix D.
%
%   Method:
%   task = 1: After cascaded inter-connection of the two systems
%
%   X1'     = A1*X1 + B1*U,
%   V       = C1*X1 + D1*U,  of order n1,
%
%   X2'     = A2*X2 + B2*V,
%   Y       = C2*X2 + D2*V,  of order n2,
%
%   where  '  denotes differentiation with respect to time,
%   the following state-space model is obtained:
%
%   X'      = A*X + B*U,
%   Y       = C*X + D*U,                                             (2)
%
%   where
%
%   A = ( A1     0  ),   B = (  B1   ),
%       ( B2*C1  A2 )        ( B2*D1 )
%
%   C = ( D2*C1  C2 ),   D = ( D2*D1 ).
%
%   This form is returned when uplo = 1. When A1 and A2 are block lower
%   triangular, the resulting state matrix is also block lower
%   triangular.
%
%   By applying a similarity transformation to the system above,
%   using the matrix  ( 0  I ),  where I is the identity matrix of
%                     ( J  0 )
%   order  n2,  and  J  is the identity matrix of order n1, the
%   system matrices become
%
%   A = ( A2  B2*C1 ),   B = ( B2*D1 ),
%       ( 0     A1  )        (  B1   )
%
%   C = ( C2  D2*C1 ),   D = ( D2*D1 ).
%
%   This form is returned when uplo = 2. When A1 and A2 are block upper
%   triangular (for instance, in the real Schur form), the resulting
%   state matrix is also block upper triangular.
%
%   task = 2: After feedback inter-connection of the two systems,
%
%   X1'     = A1*X1 + B1*U1,
%   Y1      = C1*X1 + D1*U1,  of order n1,
%
%   X2'     = A2*X2 + B2*U2,
%   Y2      = C2*X2 + D2*U2,  of order n2,
%
%   the state-space model (2) will be obtained, where
%
%   U = U1 + alpha*Y2,    X  =  ( X1 ),
%   Y = Y1 = U2,                ( X2 )
%
%   A = ( A1  -  alpha*B1*E12*D2*C1       -  alpha*B1*E12*C2    ),
%       (        B2*E21*C1            A2  -  alpha*B2*E21*D1*C2 )
%
%   B = (  B1*E12    ),
%       (  B2*E21*D1 )
%
%   C = (  E21*C1     -  alpha*E21*D1*C2 ),
%
%   D = (  E21*D1 ),
%
%   E21  =  inv( I + alpha*D1*D2 ) and
%   E12  =  inv( I + alpha*D2*D1 ) = I - alpha*D2*E21*D1.
%
%   Taking n1 = 0 and/or n2 = 0 on the function call will solve the
%   constant plant and/or constant feedback cases.
%
%   task = 3: After rowwise concatenation (parallel inter-connection
%   with separate inputs) of the two systems,
%
%   X1'     = A1*X1 + B1*U
%   Y1      = C1*X1 + D1*U
%
%   X2'     = A2*X2 + B2*V
%   Y2      = C2*X2 + D2*V
%
%   (where  '  denotes differentiation with respect to time),
%
%   with the output equation for the second system multiplied by a
%   scalar alpha, the following state-space model will be obtained:
%
%   X'      = A*X + B*(U),
%                     (V)
%
%   Y       = C*X + D*(U),
%                     (V)
%
%   where
%
%   A = ( A1   0  ),         B = ( B1   0  ),
%       ( 0    A2 )              ( 0    B2 )
%
%   C = ( C1   alpha*C2 ),   D = ( D1   alpha*D2 ).
%
%   task = 4: The matrices of the resulting systems are determined as:
%
%       ( A1   0  )             ( B1 )
%   A = (         ) ,       B = (    ) ,
%       ( 0    A2 )             ( B2 )
%
%   C = ( C1  alpha*C2 ) ,  D = D1 + alpha*D2 .
%
%   task = 5: The matrices of the resulting systems are determined as:
%
%       ( A1   0  )         ( B1  0  )
%   A = (         ) ,   B = (        ) ,
%       ( 0    A2 )         ( 0   B2 )
%
%       ( C1   0  )         ( D1  0  )
%   C = (         ) ,   D = (        ) .
%       ( 0    C2 )         ( 0   D2 )
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, July 2003.
%
% Revisions:
%   -
%
