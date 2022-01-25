function [syssvd,Q,Z,ranks] = slgsSVD(sys,job,flags,Q1,Z1)
%SLGSSVD Transform the pair (A,E) of a descriptor system to a
%        singular value decomposition (SVD) or SVD-like coordinate
%        form.
%
%        SYSSVD = SLGSSVD(SYS)  applies an orthogonal equivalence
%        transformation to the pair (A,E) of a descriptor system,
%        SYS = (A,E,B,C), so that the transformed pair (in SYSSVD),
%        (A,E) <- (Q'*A*Z,Q'*E*Z), is in a SVD coordinate form, i.e.,
%     
%                    ( A11  A12 )             ( Er  0 )
%           Q'*A*Z = (          ) ,  Q'*E*Z = (       ) ,          (1)
%                    ( A21  A22 )             (  0  0 )
%     
%        where Er is an invertible diagonal matrix having on the diagonal
%        the decreasingly ordered nonzero singular values of E. The left
%        and/or right transformations performed to reduce E (and A22, if
%        JOBA > 0, see below) are also applied to the matrices B and C.
%
%        [SYSSVD,Q,Z,RANKS] = SLGSSVD(SYS,JOB,FLAGS,Q1,Z1)  has 
%        additional input and output arguments.
%
%        JOB is an integer scalar specifying whether the SVD or
%        SVD-like coordinate form is desired, as follows:
%        JOB = 0 :  SVD coordinate form in (1);
%        JOB = 1 :  SVD-like coordinate form (1), but with Er an
%        upper triangular invertible matrix.
%        Default:  JOB = 0.
%
%        For JOB = 1, FLAGS has length 4, and contains:
%        FLAGS(1) = JOBA : specifies whether or not A22 should be
%        further reduced, as follows:
%            JOBA = 0 :  do not reduce A22;
%            JOBA = 1 :  reduce A22 to an SVD-like form (2);
%            JOBA = 2 :  reduce A22 to an upper trapezoidal form.
%        If JOBA > 0, the A22 matrix is further reduced to the 
%        SVD-like form
%     
%                    ( Ar  X )
%                    (       ) ,                                   (2)
%                    (  0  0 )
%     
%        where Ar is an invertible upper triangular matrix.
%        If JOBA = 2, then X = 0. 
%        FLAGS(2) = COMPQ : indicates what should be done with
%        matrix Q, as follows:
%           COMPQ = 0 :  do not compute Q;
%           COMPQ = 1 :  Q is initialized to the unit matrix, and
%                        the orthogonal matrix Q is returned;
%           COMPQ = 2 :  Q is initialized to an orthogonal matrix Q1
%                        and the product Q1*Q is returned.
%        If COMPQ = 0, Q must not be an output argument.
%        FLAGS(3) = COMPZ : indicates what should be done with
%        matrix Z, as follows:
%           COMPZ = 0 :  do not compute Z;
%           COMPZ = 1 :  Z is initialized to the unit matrix, and
%                        the orthogonal matrix Z is returned;
%           COMPZ = 2 :  Z is initialized to an orthogonal matrix Z1
%                        and the product Z1*Z is returned.
%        If COMPZ = 0, Z must not be an output argument.
%        FLAGS(4) = TOL  : the tolerance to be used in determining
%        the rank of E and of A22, if JOBA > 0. If TOL > 0, then the
%        given value of TOL is used as a lower bound for the reciprocal
%        condition numbers of leading submatrices of R or R22 in
%        the QR decompositions E*P = Q*R of E or A22*P22 = Q22*R22
%        of A22. A submatrix whose estimated condition number is
%        less than 1/TOL is considered to be of full rank.
%        If TOL <= 0, the default tolerance toldef = eps*prod(size(A))
%        is used instead, where eps is the machine precision. TOL < 1.
%        Default :  FLAGS = [ 0; 0; 0; 0 ].
%
%        For JOB = 0, FLAGS has length 2, and contains JOBA and TOL.
%        In this case, JOBA cannot take the value 2. If JOBA = 1,
%        the A22 matrix is further reduced to the SVD form (2), where
%        Ar is an invertible diagonal matrix having on the diagonal the
%        decreasingly ordered nonzero singular values of A22 and X = 0.
%        If TOL > 0, then singular values less than TOL*svmax are
%        treated as zero, where svmax is the maximum singular value
%        of E or of its estimate for A.
%        If TOL <= 0, the default tolerance above is used instead.
%        Default :  FLAGS = [ 0; 0 ].
%
%        RANKS is a vector with dimension at most 2.
%        RANKS(1) contains RANKE, the effective (if JOB = 0), or
%        estimated (if JOB = 1) rank of matrix E, and thus also the
%        order of the invertible diagonal or upper triangular
%        submatrix Er in (1).
%        If JOBA = 1 or 2, RANKS(2) contains RNKA22, the effective
%        (if JOB = 0), or estimated (if JOB = 1) rank of matrix
%        A22, and thus also the order of the invertible diagonal
%        or upper triangular submatrix Ar in (2).
%        If JOBA = 0, then RNKA22 is not returned.
%
%        SYSSVD = SLGSSVD(SYS)  performs the same computations,
%        but Q, Z, and RANKS are not returned; Q and Z are not 
%        computed.
%
%        Note: This function cannot work with a singular matrix E. 
%              Use GSYSTRA instead.
%        Note that GSYSTRA can work with rectangular matrices A and E.
%
%        See also SLGSBAL, SLGSHES, SLGSQRQ, SLGSRSF, GSYSTRA
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 30-04-2003.
%
%        Revisions: 03-03-2009.
%

ni = nargin;
no = nargout;
%
if ni  <= 1,  job = 0;  end
if job == 1,  nf  = 4;  else  nf = 2;  end
if ni  <= 2,  flags = zeros( nf,1 );  end
if job == 1,  
   compq = flags(2);
   compz = flags(3);
   if compq < 0 || compq > 2,  error('Wrong value for COMPQ');  end
   if compz < 0 || compz > 2,  error('Wrong value for COMPZ');  end
   if compq + compz == 4 && ni < 5,  
      error('Q1 and Z1 must be specified as input arguments');  
   elseif max( compq, compz ) == 2 && ni < 4,  
      if compq == 2,  
         error('Q1 must be specified as input argument');  
      else
         error('Z1 must be specified as input argument');  
      end
   end
   if ( ni == 5 &&     compq + compz   ~= 4 ) || ...
      ( ni == 4 && max( compq, compz ) ~= 2 ),
      warning('SLICOT:slgsSVD', 'Too many input arguments. Q1 and/or Z1 are not used.');  
   end
else
   compq = 0;
   compz = 0;
end
%
task = job + 4;
[A,B,C,D,E] = dssdata( sys );
if no <= 1,
   [A,E,B,C] = gsystra( task,A,E,B,C,flags );
else
   if job == 0,
      [A,E,B,C,Q,Z,ranks] = gsystra( task,A,E,B,C,flags );
   elseif ni <= 3,
      if max( compq, compz ) == 0,
         [A,E,B,C,Q] = gsystra( task,A,E,B,C,flags );          % Q means ranks.
      elseif min( compq, compz ) == 0,
         [A,E,B,C,Q,Z] = gsystra( task,A,E,B,C,flags );    % Q means Z, if COMPZ == 1; Z means ranks.
      else
         [A,E,B,C,Q,Z,ranks] = gsystra( task,A,E,B,C,flags );
      end
   elseif ni == 4,
      if min( compq, compz ) == 0,                            % Q means Z, Q1 means Z1, if COMPZ == 2.
         [A,E,B,C,Q,Z] = gsystra( task,A,E,B,C,flags,Q1 );    % Z means ranks.
      else
         [A,E,B,C,Q,Z,ranks] = gsystra( task,A,E,B,C,flags,Q1 );  % Q1 means Z1, if COMPZ == 2.
      end
   else
      [A,E,B,C,Q,Z,ranks] = gsystra( task,A,E,B,C,flags,Q1,Z1 );
   end
end
syssvd = dss( A,B,C,D,E,sys );
%
% end slgsSVD
