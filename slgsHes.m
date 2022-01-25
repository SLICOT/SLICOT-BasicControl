function [sysH,Q,Z] = slgsHes(sys,jobe,compq,compz,Q1,Z1)
%SLGSHES Transform the pair (A,E) of a descriptor system to a
%        generalized Hessenberg form.
%
%        SYSH = SLGSHES(SYS)  applies an orthogonal equivalence
%        transformation to the pair (A,E) of a descriptor system,
%        SYS = (A,E,B,C), so that the transformed pair (in SYSH),
%        (A,E) <- (Q'*A*Z,Q'*E*Z), is in a generalized Hessenberg form.
%        The transformation is also applied to the matrices B and C,
%        but Q and Z are not computed.
%
%        [SYSH,Q,Z] = SLGSHES(SYS,JOBE,COMPQ,COMPZ,Q1,Z1)  has 
%        additional input and output arguments.
%
%        JOBE is an integer scalar specifying whether E is a
%        general square or an upper triangular matrix, as follows:
%        JOBE = 0 :  E is a general square matrix;
%        JOBE = 1 :  E is an upper triangular matrix.
%        Default:  JOBE = 0.
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
%        Note: This function cannot work with a singular matrix E. 
%              Use GSYSTRA instead.
%
%        See also SLGSBAL, SLGSQRQ, SLGSRSF, SLGSSVD, GSYSTRA
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 30-04-2003.
%
%        Revisions: 03-03-2009.

ni = nargin;
no = nargout;
%
if ni <= 1,  jobe  = 0;  end
if ni <= 2 || no == 1,  compq = 0;  end
if ni <= 3 || no == 1,  compz = 0;  end
%
if compq < 0 || compq > 2,  error('Wrong value for COMPQ');  end
if compz < 0 || compz > 2,  error('Wrong value for COMPZ');  end
if compq + compz == 4 && ni < 6,  
   error('Q1 and Z1 must be specified as input arguments');  
elseif max( compq, compz ) == 2 && ni < 5,  
   if compq == 2,  
      error('Q1 must be specified as input argument');  
   else
      error('Z1 must be specified as input argument');  
   end
end
if ( ni == 6 &&     compq + compz   ~= 4 ) || ...
   ( ni == 5 && max( compq, compz ) ~= 2 ),
   warning('SLICOT:slgsHes', 'Too many input arguments. Q1 and/or Z1 are not used.');  
end
%
[A,B,C,D,E] = dssdata( sys );
if norm( E - eye( size( E ) ), 1 ) == 0,
   % Standard system: use just Hessenberg form.
   if max( compq, compz ) >= 1,
      [Q,A] = hess( A );  B = Q'*B;  C = C*Q;
      if compz >= 1 && no == 3,  Z = Q;     end
      if compq == 2 && no >= 2,  Q = Q1*Q;  end
      if no == 3 && compz == 2,
         if compq <= 1,  Z = Q1*Z;  else  Z = Z1*Z;  end
      end
   else
      [U,A] = hess( A );  B = U'*B;  C = C*U;
   end
   sysH = ss( A,B,C,D,sys );
%
else
   % Descriptor system: use real generalized Schur form.
   if jobe == 0, 
      if norm( tril( E, -1 ), 1 ) == 0,  jobe = 1;  end
   end
   flag = [ jobe; compq; compz ];
   %
   if no <= 1 || max( compq, compz ) == 0,
      [A,E,B,C] = gsystra( 1,A,E,B,C,flag );
   else  
      if ni <= 4,
         if min( compq, compz ) == 0,
            [A,E,B,C,Q] = gsystra( 1,A,E,B,C,flag );    % Q means Z, if COMPZ == 1.
         else
            [A,E,B,C,Q,Z] = gsystra( 1,A,E,B,C,flag );
         end
      elseif ni == 5,
         if compq*compz == 0,
            [A,E,B,C,Q] = gsystra( 1,A,E,B,C,flag,Q1 );  % Q means Z, Q1 means Z1, if COMPZ == 2.
         else
            [A,E,B,C,Q,Z] = gsystra( 1,A,E,B,C,flag,Q1 );  % Q1 means Z1, if COMPZ == 2.
         end
      else
         [A,E,B,C,Q,Z] = gsystra( 1,A,E,B,C,flag,Q1,Z1 );
      end
   end
   sysH = dss( A,B,C,D,E,sys );
end
%
% end slgsHes
