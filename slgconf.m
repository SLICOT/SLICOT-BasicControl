function [syscf,nc,s,Q,Z] = slgconf(sys,job,tol,compq,compz,Q1,Z1)
%SLGCONF Transform a descriptor system (A-lambda*E,B,C) to the
%        controllability staircase form.
%
%        [SYSCF,NC,S] = SLGCONF(SYS)  computes orthogonal transformation
%        matrices Q and Z which reduce the n-th order descriptor system
%        (A-lambda*E,B,C) to the form
%     
%                    ( Ac  *  )             ( Ec  *  )           ( Bc )
%           Q'*A*Z = (        ) ,  Q'*E*Z = (        ) ,  Q'*B = (    ) ,
%                    ( 0  Anc )             ( 0  Enc )           ( 0  )
%     
%              C*Z = ( Cc Cnc ) ,                                     (1)
%     
%        where the c-th order descriptor system (Ac-lambda*Ec,Bc,Cc) is
%        finite and/or infinite controllable, hence c is the order of the
%        controllable part of the pair (A-lambda*E,B). The pencil
%        Anc - lambda*Enc is regular of order n-c and contains the
%        uncontrollable finite and/or infinite eigenvalues of the pencil
%        A-lambda*E. The reduced order descriptor system (Ac-lambda*Ec,Bc,Cc)
%        has the same transfer-function matrix as the original system
%        (A-lambda*E,B,C). The matrices Ac and Ec have a staircase form 
%        (see (2) and (3) below). The transformed system is returned in the
%        LTI object SYSCF. The size of the controllable subsystem, and the
%        block sizes in the staircase form of the matrices Ac and Ec, are
%        also returned in the vectors NC and S, respectively, but the
%        matrices Q and Z are not computed.
%
%        [SYSCF,NC,S,Q,Z] = SLGCONF(SYS,JOB,TOL,COMPQ,COMPZ,Q1,Z1)  has
%        additional input and output arguments.
%
%        JOB is an integer scalar specifying what to do, as follows: 
%        JOB = 0 :  separate both finite and infinite uncontrollable
%                   eigenvalues;
%        JOB = 1 :  separate only finite uncontrollable eigenvalues;
%        JOB = 2 :  separate only nonzero finite and infinite
%                   uncontrollable eigenvalues;
%        Default:  JOB = 0.
%
%        TOL is the tolerance to be used in rank determinations when
%        transforming (A-lambda*E,B). If TOL > 0, then the given value
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
%        S is an integer vector containing rtau, where rtau(i), 
%        for i = 1, ..., k, is the row dimension of the full row
%                   _         _
%        rank block Ei,i-1 or Ai,i-1 in the staircase form (2)
%        or (4) for JOB = 0 or 2, or for JOB = 1, respectively.
%
%        For JOB = 0 or 2, the pencil ( Bc Ec-lambda*Ac ) has full row
%        rank c for all finite lambda and is in a staircase form with
%                         _      _          _        _
%                       ( E1,0   E1,1  ...  E1,k-1   E1,k  )
%                       (        _          _        _     )
%           ( Bc Ec ) = (  0     E2,1  ...  E2,k-1   E2,k  ) ,       (2)
%                       (              ...  _        _     )
%                       (  0       0   ...  Ek,k-1   Ek,k  )
%     
%                         _          _        _
%                       ( A1,1  ...  A1,k-1   A1,k  )
%                       (            _        _     )
%             Ac      = (   0   ...  A2,k-1   A2,k  ) ,              (3)
%                       (       ...           _     )
%                       (   0   ...    0      Ak,k  )
%              _
%        where Ei,i-1 is an rtau(i)-by-rtau(i-1) full row rank matrix
%                               _
%        (with rtau(0) = m) and Ai,i is an rtau(i)-by-rtau(i) upper
%        triangular matrix, hence Ac is also upper triangular. The number c
%        is stored in NC(1), and the vector rtau with elements rtau(i),
%        for i = 1 : k, is stored in S.
%     
%        For JOB = 1, the pencil ( Bc Ac-lambda*Ec ) has full row
%        rank c for all finite lambda and is in a staircase form with
%                         _     _          _        _
%                       ( A1,0  A1,1  ...  A1,k-1   A1,k  )
%                       (       _          _        _     )
%           ( Bc Ac ) = (  0    A2,1  ...  A2,k-1   A2,k  ) ,        (4)
%                       (             ...  _        _     )
%                       (  0      0   ...  Ak,k-1   Ak,k  )
%     
%                         _          _        _
%                       ( E1,1  ...  E1,k-1   E1,k  )
%                       (            _        _     )
%             Ec      = (   0   ...  E2,k-1   E2,k  ) ,              (5)
%                       (       ...           _     )
%                       (   0   ...    0      Ek,k  )
%              _
%        where Ai,i-1 is an rtau(i)-by-rtau(i-1) full row rank matrix
%                               _
%        (with rtau(0) = m) and Ei,i is an rtau(i)-by-rtau(i) upper
%        triangular and nonsingular matrix, hence Ec is also upper 
%        triangular and nonsingular.
%     
%        The pairs (Bc,Ec) in (2) and (Bc,Ac) in (4) are in the 
%        controllability staircase form.
%     
%        For JOB = 0, the (n-c)-by-(n-c) regular pencil Anc - lambda*Enc
%        has the form
%     
%                              ( Ainc - lambda*Einc         *          )
%           Anc - lambda*Enc = (                                       ) ,
%                              (        0           Afnc - lambda*Efnc )
%     
%        where:
%          1) the inc-by-inc regular pencil Ainc - lambda*Einc, with Ainc
%             upper triangular and nonsingular, and Einc nilpotent, contains
%             the uncontrollable infinite eigenvalues of A - lambda*E; the
%             value inc is the number of uncontrollable infinite eigenvalues
%             of the pencil A - lambda*E, and it is stored in NC(2).
%          2) the (n-c-inc)-by-(n-c-inc) regular pencil Afnc - lambda*Efnc,
%             with Efnc upper triangular and nonsingular, contains the
%             uncontrollable finite eigenvalues of A - lambda*E.
%     
%        Note: The significance of the two diagonal blocks can be
%              interchanged by calling the function with the arguments A
%              and E interchanged. In this case, Ainc - lambda*Einc contains
%              the uncontrollable zero eigenvalues of A - lambda*E, while
%              Afnc - lambda*Efnc contains the uncontrollable nonzero finite
%              and infinite eigenvalues of A - lambda*E.
%     
%        For JOB = 1, the pencil Anc - lambda*Enc has the form
%     
%           Anc - lambda*Enc = Afnc - lambda*Efnc ,
%     
%        where the regular pencil Afnc - lambda*Efnc, with Efnc upper
%        triangular and nonsingular, contains the uncontrollable finite
%        eigenvalues of A - lambda*E.
%     
%        For JOB = 2, the pencil Anc - lambda*Enc has the form
%     
%           Anc - lambda*Enc = Ainc - lambda*Einc ,
%     
%        where the regular pencil Ainc - lambda*Einc, with Ainc (hence Anc)
%        upper triangular and nonsingular, contains the uncontrollable
%        nonzero finite and infinite eigenvalues of A - lambda*E.
%
%        Note: This function cannot work with a singular matrix E.
%              Use GSYSCOM instead.
%
%        See also SLGMINR, SLGOBSF, GSYSCOM
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 15-05-2003, 03-03-2009.
%

ni = nargin;
no = nargout;
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
   warning('SLICOT:slgconf', 'Too many input arguments. Q1 and/or Z1 are not used.');  
end
%
[A,B,C,D,E] = dssdata( sys );
if norm( E - eye( size( E ) ), 1 ) == 0,
   % Standard system.
   [syscf,nc,Q,s] = slconf(sys,tol);
   if compq == 2,  Q = Q1*Q;  end
   if no  == 5,  Z = Q;  end
   if job == 0,  nc(2) = 0;  end
   %
else
   % Descriptor system.
   flag = [ job; compq; compz; tol ];
   %
   if no <= 3 || max( compq, compz ) == 0,
      [A,E,B,C,nc,s] = gsyscom( 1,A,E,B,C,flag );
   else  
      if ni <= 5,
         if min( compq, compz ) == 0,
            [A,E,B,C,Q,nc,s] = gsyscom( 1,A,E,B,C,flag );    % Q means Z, if COMPZ == 1.
         else
            [A,E,B,C,Q,Z,nc,s] = gsyscom( 1,A,E,B,C,flag );
         end
      elseif ni == 6,
         if compq*compz == 0,
            [A,E,B,C,Q,nc,s] = gsyscom( 1,A,E,B,C,flag,Q1 );  % Q means Z, Q1 means Z1, if COMPZ == 2.
         else
            [A,E,B,C,Q,Z,nc,s] = gsyscom( 1,A,E,B,C,flag,Q1 );  % Q1 means Z1, if COMPZ == 2.
         end
      else
         [A,E,B,C,Q,Z,nc,s] = gsyscom( 1,A,E,B,C,flag,Q1,Z1 );
      end
   end
   syscf = dss( A,B,C,D,E,sys );
   if job == 0,  nc = nc(1:2);  else  nc = nc(1);  end
end
%
% end slgconf
