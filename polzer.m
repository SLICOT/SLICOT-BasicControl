function [ p, z, rnk, infz, Kronl, Kronr, infe, niz, Af, Ef ] = polzer( sys, tol, bal )
%POLZER  Computes the normal rank, poles, zeros, and the Kronecker
%        structure of the system pencil for a standard or descriptor
%        system, SYS = (A,B,C,D,E), with possibly rectangular
%        matrices A and E.
%
%        [P,Z] = POLZER(SYS)  computes the poles P and zeros Z of the 
%        system SYS. 
%        
%        [P,Z,RNK,INFZ,KRONL,KRONR,INFE,NIZ,AF,EF] = POLZER(SYS,TOL,BAL)  
%        has additional input and output parameters.
%        
%        TOL is the tolerance to be used in rank decisions to determine
%        the effective rank, which is defined as the order of the
%        largest leading (or trailing) triangular submatrix in the
%        QR (or RQ) factorization with column (or row) pivoting
%        whose estimated condition number is less than 1/TOL.
%        If TOL <= 0, then default tolerances are used instead,
%        defined in terms of the size of the system matrix and eps,
%        where eps is the machine precision.
%        Default: TOL = 0.
%
%        BAL is an integer indicating whether the system should be
%        balanced (scaled).
%        BAL = 0 :  use balancing;
%        BAL = 1 :  do not use balancing.
%        Default: BAL = 0.
%        
%        RNK is the normal rank of the transfer-function matrix of
%        the system.
%
%        INFZ is an integer vector containing information on the
%        infinite elementary divisors as follows: the system has INFZ(i)
%        infinite elementary divisors of degree i (in the Smith form),
%        where i = 1,2,...,lenght(INFZ).
%
%        KRONL and KRONR are integer vectors containing the left
%        Kronecker (row) and right Kronecker (column) indices,
%        respectively.
%
%        INFE is an integer vector containing the multiplicities of
%        infinite zeros (for a descriptor system).
%
%        NIZ is the number of infinite zeros (for a descriptor system).
%
%        AF and BF are the matrices Af and Ef of the reduced pencil,
%        which define the finite zeros of the system.
%
%        Note: This function cannot work with a singular matrix E,
%              or with rectangular A and E. Use POLEZERO(Z) instead.
%
%        See also POLEZERO, POLEZEROZ
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Dec-22-2003.
%        Revised Jul-14-2004; Dec-14-2008; Mar-03-2009.
%

%
ni = nargin;
if ni < 1,  
    error( ['Usage: [P,Z] = POLZER(SYS)', sprintf('\n'), ...
            '       [P,Z,RNK,INFZ,KRONL,KRONR,INFE,NIZ,AF,EF] = POLZER(SYS,TOL,BAL)' ] );
end
%
if ni == 1,
   tol = 0;  
   bal = 0;  
elseif ni == 2,
   bal = 0;  
end
%
if isa( sys, 'ss' ),
   E = sys.e;
   iscmplx = any( any( imag( sys.a ) ) ) | any( any( imag( sys.b ) ) ) | ...
             any( any( imag( sys.c ) ) ) | any( any( imag( sys.d ) ) ) | ...
             any( any( imag( E ) ) );
   if iscmplx,
      if ~isempty( E ) && ~isequal( E, eye( size( sys.a ) ) ),
         [ poles, betap, zros, betaz, rnk, infz, Kronl, Kronr, infe, niz, Af, Ef ] ...
              = polezeroz( 4, sys.a, sys.b, sys.c, sys.d, E, tol, bal );
         p = poles ./ betap;
      else  
         [ p, zros, betaz, rnk, infz, Kronl, Kronr, infe, niz, Af, Ef ] ...
              = polezeroz( 4, sys.a, sys.b, sys.c, sys.d, [ ], tol, bal );
      end
      z = zros ./ betaz;
   else
      if ~isempty( E ) && ~isequal( E, eye( size( sys.a ) ) ),
         [ rpoles, ipoles, betap, rzeros, izeros, betaz, rnk, infz, ...
                Kronl, Kronr, infe, niz, Af, Ef ] ...
              = polezero( 4, sys.a, sys.b, sys.c, sys.d, E, tol, bal );
         p = complex( rpoles, ipoles ) ./ betap;
      else
         [ rpoles, ipoles, rzeros, izeros, betaz, rnk, infz, ...
                Kronl, Kronr, infe, niz, Af, Ef ] ...
              = polezero( 4, sys.a, sys.b, sys.c, sys.d, [ ], tol, bal );
         p = complex( rpoles, ipoles );
      end
      z = complex( rzeros, izeros ) ./ betaz;
   end
else
    error( 'SYS must be an ss/dss system object' );
end
%
% end polzer
