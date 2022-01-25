function [X,sep] = slgest(A,E,C,flag,trans)
%SLGEST  Solve generalized Stein equations.
%
%        X = SLGEST(A,E,C,FLAG,TRANS)  computes the unique symmetric
%        solution X of a generalized Stein (discrete-time Lyapunov)
%        equation
%
%                op(A)'*X*op(A) - op(E)'*X*op(E) = C,
%
%        where op(M) = M or M', and C is a symmetric matrix.
%
%        FLAG is a scalar characterizing the structures of A and E:
%        FLAG = 1 : (A,E) is in generalized Schur form;
%                   otherwise, (A,E) is in general form.
%        Default: FLAG = 0.
%
%        TRANS specifies if op(M) = M or M':
%        TRANS = 0 : op(M) = M;
%                    otherwise, op(M) = M'.
%        Default: TRANS = 0.
%
%        [X,SEP] = SLGEST(A,E,C,FLAG,TRANS)  computes the solution X,
%        as well as the separation, sep(A,E).
%
%        See also GENLEQ, SLGELY, SLGESG, SLGSLY, SLGSST
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 30-04-1999.
%
%        Revisions: 03-03-2009.

ni = nargin;
fl = [1;0];
%
if ni >= 4
   [mf,nf] = size(flag);
   if mf > 1 || nf > 1
      error('flag must be a scalar.')
   end
   fl(2) = flag;
end
%
if ni < 5
   trans = 0;
end
%
if nargout < 2
   X = genleq(3,A,E,C,fl,trans);
else
   [X,sep] = genleq(3,A,E,C,fl,trans);
end
%
% end slgest
