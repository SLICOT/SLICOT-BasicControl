function [X,sep] = sllyap(A,C,flag,trans)
%SLLYAP  Solve continuous-time Lyapunov equations.
% 
%        X = SLLYAP(A,C,FLAG,TRANS)  computes the unique symmetric
%        solution X of a (continuous-time) Lyapunov equation
%
%                op(A)'*X + X*op(A) = C,
%
%        where op(M) = M or M', and C is a symmetric matrix.
%
%        FLAG is a scalar characterizing the structure of A:
%        FLAG = 1 : A is quasi upper triangular;
%                   otherwise, A is a general matrix.
%        Default: FLAG = 0.
%
%        TRANS specifies if op(M) = M or M':
%        TRANS = 0 : op(M) = M;
%                    otherwise, op(M) = M'.
%        Default: TRANS = 0.
%
%        [X,SEP] = SLLYAP(A,C,FLAG,TRANS)  computes the solution X,
%        as well as the separation, sep(A', -A).
%
%        See also LINMEQ, SLDSYL, SLSTEI, SLSTLY, SLSTST, SLSYLV
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        H. Xu 7-12-1998. 
%        Revised V. Sima 30-04-1999, 29-04-2000, Mar. 3, 2009.
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
if nargout < 2
   X = linmeq(2,A,C,fl,trans);
else
   [X,sep] = linmeq(2,A,C,fl,trans);
end
%
% end sllyap
