function [Q,R,X,Y] = fstqr(TC,TR,B,C)
%FSTQR   Compute the orthogonal-triangular decomposition of a (block)
%        Toeplitz matrix T and solve associated linear least-squares
%        problems, given the first block column TC and the first block
%        row TR of T. It is assumed that the first MIN(SIZE(T)) columns
%        of T have full rank.
%
%        R = FSTQR(TC,TR)  computes the upper Cholesky factor R of T'*T,
%        i.e., R is upper triangular and R'*R = T'*T.
%
%        [Q,R] = FSTQR(TC,TR)  also computes a factor Q such that
%        Q'*Q = I and T = Q*R.
%
%        [R,X]   = FSTQR(TC,TR,B)  also computes a matrix X such that
%        [Q,R,X] = FSTQR(TC,TR,B)  the Euclidean norms of the columns of
%                                  B-T*X are minimized.
%
%        [R,X,Y]   = FSTQR(TC,TR,B,C)  also computes a matrix Y such
%        [Q,R,X,Y] = FSTQR(TC,TR,B,C)  that T'*Y = C and the Euclidean
%        [R,Y]     = FSTQR(TC,TR,[],C) norms of the columns of Y are
%        [Q,R,Y]   = FSTQR(TC,TR,[],C) minimized.
%
%        When TC or TR is sparse and the factor Q is not required, a
%        banded factorization is used and a banded sparse Cholesky
%        factor R is returned.
%
%        See also FSTLSQ, FSTMUL.
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        D. Kressner 10-06-2002.
%        Revised: V. Sima, 03-03-2009. 
%

ni = nargin;  nout = nargout;
%
if ni < 2,
   error('Not enough input arguments.');
end
if ~(issparse(TC) || issparse(TR)) || (nout >= ni),
   if issparse(TC),
      TC = full(TC);
   end
   if issparse(TR),
      TR = full(TR);
   end
   if ni == 2,
      if nout <= 1,
         Q = fstoeq(1,TC,TR);
         Q = Q';
      elseif nout == 2,
         [Q,R] = fstoeq(1,TC,TR);
         R = R';
      else
         error('Improper number of output arguments.')
      end
   elseif ni == 3,
      if nout == 2,
         [Q,R] = fstoeq(2,TC,TR,B);
         Q = Q';
      elseif nout == 3,
         [Q,R,X] = fstoeq(2,TC,TR,B);
         R = R';
      else
         error('Improper number of output arguments.')
      end
   elseif ni == 4,
      if (isempty(B)) && nout == 2,
         [Q,R] = fstoeq(3,TC,TR,C);
         Q = Q';
      elseif (isempty(B)) && nout == 3,
         [Q,R,X] = fstoeq(3,TC,TR,C);
         R = R';
      elseif nout == 3,
         [Q,R,X] = fstoeq(4,TC,TR,B,C);
         Q = Q';
      elseif nout == 4,
         [Q,R,X,Y] = fstoeq(4,TC,TR,B,C);
         R = R';
      else
         error('Improper number of output arguments.')
      end    
   else
      error('Improper number of input arguments.')
   end
else
   if ~isnumeric(TC) || ~isreal(TC),
      error('TC must be a real matrix');
   end
   if ~isnumeric(TR) || ~isreal(TR),
      error('TR must be a real matrix');
   end
   l = size(TC,2);
   k = size(TR,1);
   m = size(TC,1);
   n = size(TR,2);
   if min([l k m n]) <= 0 || rem(m,k) || rem(n,l),
      error('Dimensions of TC and TR do not match');
   end
   m = m/k;
   n = n/l;
   [i,j] = find(TC);
   mb = ceil(max(i)/k);
   [i,j] = find(TR);
   nb = ceil(max(j)/l);
   TCB = full(TC(1:mb*k,:));
   TRB = full(TR(:,1:nb*l));
   if ni == 2,
      if nout <= 1,
         RB = fstoeq(6,m,n,TCB,TRB);
      else
         error('Improper number of output arguments.')
      end
    elseif ni == 3,
      if nout == 2,
         [RB,R] = fstoeq(7,m,n,TCB,TRB,B);
      else
         error('Improper number of output arguments.')
      end
   elseif ni == 4,
      if (isempty(B)) && nout == 2,
         [RB,R] = fstoeq(8,m,n,TCB,TRB,C);
      elseif nout == 3,
         [RB,R,X] = fstoeq(9,m,n,TCB,TRB,B,C);
      else
         error('Improper number of output arguments.')
      end    
   else
      error('Improper number of input arguments.')
   end
   % Create a sparse matrix from its diagonals
   nrr = min(m*k,n*l);
   ncr = n*l;
   RB = RB';
   ndiag = size(RB,2);
   d = 0:ndiag-1;
   a = zeros(ndiag,3);  i1 = 1;
   for j = 1:ndiag
      i = (1:min(nrr,ncr-d(j)))';  i2 = i1 + i(end) - 1;
      a(i1:i2,:) = [i i+d(j) RB(i,j)];
      i1 = i2 + 1;
   end
   Q = sparse(a(:,1),a(:,2),full(a(:,3)),nrr,ncr);
end

%
% end fstqr
