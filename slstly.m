function X = slstly(A,C,flag,trans)
%SLSTLY  Solve stable continuous-time Lyapunov equations.
% 
%        X = SLSTLY(A,C,FLAG,TRANS)  computes a Cholesky factor X
%        of the unique symmetric positive semi-definite solution
%        op(X)'*op(X) of a stable (continuous-time) Lyapunov equation
%
%           op(A)'*op(X)'*op(X) + op(X)'*op(X)*op(A) = -op(C)'*op(C),
%
%        where op(M) = M or M', op(C) is a general matrix with the same
%        number of columns as A, and matrix X is upper triangular.
%
%        FLAG is a scalar characterizing the structure of A:
%        FLAG =  1 : A is quasi upper triangular;
%                    otherwise, A is in general form.
%        Default: FLAG = 0.
%
%        TRANS specifies if op(M) = M or M':
%        TRANS = 0 : op(M) = M;
%                    otherwise, op(M) = M'.
%        Default: TRANS = 0.
%
%        See also LINMEQ, SLDSYL, SLLYAP, SLSTEI, SLSTST, SLSYLV
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        H. Xu 7-12-1998. 
%        Revised V. Sima 30-04-1999, 29-04-2000, 03-03-2009.
%
%        Comments:
%        Matrix A must be stable, i.e., all eigenvalues must have
%        negative real parts.
%

ni = nargin;
fl = [0;0];
%
if ni >= 3
   [mf,nf] = size(flag);
   if mf > 1 || nf > 1
      error('flag must be a scalar.')
   end
   fl(2) = flag;
end            
%
if ni < 4
   trans = 0;
end
%
X = linmeq(3,A,C,fl,trans);
%
% end slstly
