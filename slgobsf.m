function [sysob,no,s,Q,Z] = slgobsf(sys,job,tol,compq,compz,Q1,Z1)
%SLGOBSF Transform a descriptor system (A-lambda*E,B,C) to the
%        observability staircase form.
%
%        [SYSOB,NO,S] = SLGOBSF(SYS)  computes orthogonal transformation
%        matrices Q and Z which reduce the n-th order descriptor system
%        (A-lambda*E,B,C) to the form
%     
%                    ( Ano  * )             ( Eno  * )           ( Bno )
%           Q'*A*Z = (        ) ,  Q'*E*Z = (        ) ,  Q'*B = (     ) ,
%                    ( 0   Ao )             ( 0   Eo )           ( Bo  )
%  
%           C*Z = ( 0   Co ) ,                                        (1)
%     
%        where the o-th order descriptor system (Ao-lambda*Eo,Bo,Co) is
%        finite and/or infinite observable, hence o is the order of the
%        observable part of the pair (A'-lambda*E',C')'. The pencil
%        Ano - lambda*Eno is regular of order n-o and contains the
%        unobservable finite and/or infinite eigenvalues of the pencil
%        A-lambda*E. The reduced order descriptor system (Ao-lambda*Eo,Bo,Co)
%        has the same transfer-function matrix as the original system
%        (A-lambda*E,B,C). The matrices Ao and Eo have a staircase form 
%        (see (2) and (3) below). The transformed system is returned in the
%        LTI object SYSOB. The size of the observable subsystem, and the
%        block sizes in the staircase form of the matrices Ao and Eo, are
%        also returned in the vectors NO and S, respectively, but the
%        matrices Q and Z are not computed.
%
%        [SYSOB,NO,S,Q,Z] = SLGOBSF(SYS,JOB,TOL,COMPQ,COMPZ,Q1,Z1)  has
%        additional input and output arguments.
%
%        JOB is an integer scalar specifying what to do, as follows: 
%        JOB = 0 :  separate both finite and infinite unobservable
%                   eigenvalues;
%        JOB = 1 :  separate only finite unobservable eigenvalues;
%        JOB = 2 :  separate only nonzero finite and infinite
%                   unobservable eigenvalues;
%        Default:  JOB = 0.
%
%        TOL is the tolerance to be used in rank determinations when
%        transforming (A'-lambda*E',C')'. If TOL > 0, then the given value
%        of TOL is used as a lower bound for reciprocal condition numbers;
%        a (sub)matrix whose estimated condition number is less than
%        1/tol is considered to be of full rank.
%        If TOL <= 0, the default tolerance toldef = eps*prod(size(A))
%        is used instead, where eps is the machine precision. TOL < 1.
%        Default :  TOL = 0.
%
%        COMPQ is an integer scalar indicating what should be done
%        with matrix Q, as follows:
%        COMPQ = 0 :  do not compute Q;
%        COMPQ = 1 :  Q is initialized to the unit matrix, and
%                     the orthogonal matrix Q is returned;
%        COMPQ = 2 :  Q is initialized to an orthogonal matrix Q1
%                     and the product Q1*Q is returned.
%        Default:  COMPQ = 0.
%        If COMPQ = 0, Q must not be an output argument.
%
%        COMPZ : is an integer scalar indicating what should be done
%        with matrix Z, as follows:
%        COMPZ = 0 :  do not compute Z;
%        COMPZ = 1 :  Z is initialized to the unit matrix, and
%                     the orthogonal matrix Z is returned;
%        COMPZ = 2 :  Z is initialized to an orthogonal matrix Z1
%                     and the product Z1*Z is returned.
%        Default:  COMPZ = 0.
%        If COMPZ = 0, Z must not be an output argument.
%
%        S is an integer vector containing ctau, where ctau(i), 
%        for i = 1, ..., k, is the column dimension of the full column
%                   _         _
%        rank block Ei-1,i or Ai-1,i in the staircase form (2)
%        or (4) for JOB = 0 or 2, or for JOB = 1, respectively.
%
%        For JOB = 0 or 2, the pencil ( Eo-lambda*Ao ) has full column
%                                     (      Co      )
%        rank o for all finite lambda and is in a staircase form with
%                         _      _            _      _
%                       ( Ek,k   Ek,k-1   ... Ek,2   Ek,1   )
%                       ( _      _            _      _      )
%           ( Eo )  =   ( Ek-1,k Ek-1,k-1 ... Ek-1,2 Ek-1,1 ) ,      (2)
%           ( Co )      (     ...         ... _      _      )
%                       (  0       0      ... E1,2   E1,1   )
%                       (                            _      )
%                       (  0       0      ... 0      E0,1   )
%                         _          _      _
%                       ( Ak,k  ...  Ak,2   Ak,1 )
%                       (       ...  _      _    )
%             Ao    =   (   0   ...  A2,2   A2,1 ) ,                 (3)
%                       (                   _    )
%                       (   0   ...    0    A1,1 )
%              _
%        where Ei-1,i is a ctau(i-1)-by-ctau(i) full column rank matrix
%                               _
%        (with ctau(0) = p) and Ai,i is a ctau(i)-by-ctau(i) upper
%        triangular matrix, hence Ao is also upper triangular. The number o
%        is stored in NO(1), and the vector ctau with elements ctau(i),
%        for i = 1 : k, is stored in S.
%     
%        For JOB = 1, the pencil ( Ao-lambda*Eo ) has full column
%                                (      Co      )
%        rank o for all finite lambda and is in a staircase form with
%                         _      _            _      _
%                       ( Ak,k   Ak,k-1   ... Ak,2   Ak,1   )
%                       ( _      _            _      _      )
%           ( Ao )  =   ( Ak-1,k Ak-1,k-1 ... Ak-1,2 Ak-1,1 ) ,      (4)
%           ( Co )      (     ...         ... _      _      )
%                       (  0       0      ... A1,2   A1,1   )
%                       (                            _      )
%                       (  0       0      ... 0      A0,1   )
%                         _          _      _
%                       ( Ek,k  ...  Ek,2   Ek,1 )
%                       (       ...  _      _    )
%             Eo    =   (   0   ...  E2,2   E2,1 ) ,                 (5)
%                       (                   _    )
%                       (   0   ...    0    E1,1 )
%              _
%        where Ai-1,i is a ctau(i-1)-by-ctau(i) full column rank matrix
%                               _
%        (with ctau(0) = p) and Ei,i is a ctau(i)-by-ctau(i) upper
%        triangular and nonsingular matrix, hence Eo is also upper 
%        triangular and nonsingular.
%     
%        The pairs ( Eo ) in (2) and ( Ao ) in (4) are in the 
%                  ( Co )            ( Co )
%        observability staircase form.
%     
%        For JOB = 0, the (n-o)-by-(n-o) regular pencil Ano - lambda*Eno
%        has the form
%
%                              ( Afno - lambda*Efno         *          )
%           Ano - lambda*Eno = (                                       ) ,
%                              (        0           Aino - lambda*Eino )
%
%        where:
%          1) the ino-by-ino regular pencil Aino - lambda*Eino, with Aino
%             upper triangular and nonsingular, and Eino nilpotent, contains
%             the unobservable infinite eigenvalues of A - lambda*E; the
%             value ino is the number of unobservable infinite eigenvalues
%             of the pencil A - lambda*E, and it is stored in NO(2).
%          2) the (n-o-ino)-by-(n-o-ino) regular pencil Afno - lambda*Efno,
%             with Efno upper triangular and nonsingular, contains the
%             unobservable finite eigenvalues of A - lambda*E.
%
%        Note: The significance of the two diagonal blocks can be
%              interchanged by calling the function with the arguments A
%              and E interchanged. In this case, Aino - lambda*Eino contains
%              the unobservable zero eigenvalues of A - lambda*E, while
%              Afno - lambda*Efno contains the unobservable nonzero finite
%              and infinite eigenvalues of A - lambda*E.
%
%        For JOB = 1, the pencil Ano - lambda*Eno has the form
%
%           Ano - lambda*Eno = Afno - lambda*Efno ,
%
%        where the regular pencil Afno - lambda*Efno, with Efno upper
%        triangular and nonsingular, contains the unobservable finite
%        eigenvalues of A - lambda*E.
%
%        For JOB = 2, the pencil Ano - lambda*Eno has the form
%
%           Ano - lambda*Eno = Aino - lambda*Eino ,
%
%        where the regular pencil Aino - lambda*Eino, with Aino (hence Ano) 
%        upper triangular and nonsingular, contains the unobservable
%        nonzero finite and infinite eigenvalues of A - lambda*E.
%
%        Note: This function cannot work with a singular matrix E.
%              Use GSYSCOM instead.
%
%        See also SLGCONF, SLGMINR, GSYSCOM
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 15-05-2003.
%
%        Revisions: 03-03-2009.

ni  = nargin;
nou = nargout;
%
if ni == 1,
    job   = 0;
    tol   = 0;
    compq = 0;
    compz = 0;
elseif ni == 2,
    tol   = 0;
    compq = 0;
    compz = 0;
elseif ni == 3,
    compq = 0;
    compz = 0;
elseif ni == 4,
    compz = 0;
end
%
if compq < 0 || compq > 2,  error('Wrong value for COMPQ');  end
if compz < 0 || compz > 2,  error('Wrong value for COMPZ');  end
if compq + compz == 4 && ni < 7,  
   error('Q1 and Z1 must be specified as input arguments');  
elseif max( compq, compz ) == 2 && ni < 6,  
   if compq == 2,  
      error('Q1 must be specified as input argument');  
   else
      error('Z1 must be specified as input argument');  
   end
end
if ( ni == 7 &&     compq + compz   ~= 4 ) || ...
   ( ni == 6 && max( compq, compz ) ~= 2 ),
   warning('SLICOT:slgobsf', 'Too many input arguments. Q1 and/or Z1 are not used.');  
end
%
[A,B,C,D,E] = dssdata( sys );
if norm( E - eye( size( E ) ), 1 ) == 0,
   % Standard system.
   [sysob,no,Q,s] = slobsf(sys,tol);
   if compq == 2,  Q = Q1*Q;  end
   if nou == 5,  Z = Q;  end
   if job == 0,  no(2) = 0;  end
   %
else
   % Descriptor system.
   flag = [ job; compq; compz; tol ];
   %
   if nou <= 3 || max( compq, compz ) == 0,
      [A,E,B,C,no,s] = gsyscom( 2,A,E,B,C,flag );
   else  
      if ni <= 5,
         if min( compq, compz ) == 0,
            [A,E,B,C,Q,no,s] = gsyscom( 2,A,E,B,C,flag );    % Q means Z, if COMPZ == 1.
         else
            [A,E,B,C,Q,Z,no,s] = gsyscom( 2,A,E,B,C,flag );
         end
      elseif ni == 6,
         if compq*compz == 0,
            [A,E,B,C,Q,no,s] = gsyscom( 2,A,E,B,C,flag,Q1 );  % Q means Z, Q1 means Z1, if COMPZ == 2.
         else
            [A,E,B,C,Q,Z,no,s] = gsyscom( 2,A,E,B,C,flag,Q1 );  % Q1 means Z1, if COMPZ == 2.
         end
      else
         [A,E,B,C,Q,Z,no,s] = gsyscom( 2,A,E,B,C,flag,Q1,Z1 );
      end
   end
   sysob = dss( A,B,C,D,E,sys );
   if job == 0,  no = no(1:2);  else  no = no(1);  end
end
%
% end slgobsf
