function [X] = slgsly(A,E,C,flag,trans)
%SLGSLY  Solve stable generalized continuous-time Lyapunov equations.
%
%        X = SLGSLY(A,E,C,FLAG,TRANS)  computes a Cholesky factor X
%        of the unique symmetric positive semi-definite solution
%        op(X)'*op(X) of a stable generalized continuous-time Lyapunov
%        equation
%
%           op(A)'*op(X)'*op(X)*op(E) + op(E)'*op(X)'*op(X)*op(A) =
%                                                   -op(C)'*op(C),
%
%        where op(M) = M or M', op(C) is a general matrix with the same
%        number of columns as A, and matrix X is upper triangular.
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
%        See also GENLEQ, SLGELY, SLGESG, SLGEST, SLGSST
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 30-04-1999.
%
%        Revisions: 03-03-2009.
%
%        Comments
%        The pencil (A,E) must be stable, i.e., all eigenvalues must
%        have negative real parts.
%

ni = nargin;
fl = [0;0];
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
X = genleq(4,A,E,C,fl,trans);
%
% end slgsly
