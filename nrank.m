function rnk = nrank( sys, tol, bal )
%NRANK   Computes the normal rank of the transfer-function matrix
%        of a standard system, SYS = (A,B,C,D).
%
%        RNK = NRANK(SYS)  computes the normal rank of the system SYS. 
%        
%        RNK = NRANK(SYS,TOL,BAL)  has  additional input parameters.
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
%        See also POLEZERO, POLEZEROZ
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Dec-22-2003.
%        Revised Dec-13-2008, Mar-03-2009.
%

%
ni = nargin;
if ni < 1,  
    error( ['Usage: RNK = NRANK(SYS)', sprintf('\n'), ...
            '       RNK = NRANK(A,TOL,BAL)' ] );
end
%
if ni == 1,
   tol = 0;  
   bal = 0;  
elseif ni == 2,
   bal = 0;  
end
%
if any( any( imag( sys.a ) ) ) || any( any( imag( sys.b ) ) ) || ...
   any( any( imag( sys.c ) ) ) || any( any( imag( sys.d ) ) ),
   rnk = polezeroz( 1, sys.a, sys.b, sys.c, sys.d, tol, bal );
else
   rnk = polezero( 1, sys.a, sys.b, sys.c, sys.d, tol, bal );
end
%
% end nrank
