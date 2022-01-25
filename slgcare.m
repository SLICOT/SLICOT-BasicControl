function  [X,F,ev,rcond1] = slgcare(A,E,Q,R,B,L,flag)
%SLGCARE  Solve descriptor CARE with generalized Schur method.
%
%           [X,F,EV,RCOND1] = SLGCARE(A,E,Q,R,B,L,FLAG)  computes the
%         symmetric solution X of the continuous-time algebraic Riccati
%         equation
%                                                   -1
%           0 = Q + A'*X*E + E'*X*A - (L + E'*X*B)*R  *(L + E'*X*B)',  (1)
%         
%         the feedback gain matrix F, the associated eigenvalues EV,
%         as well as RCOND1, an estimate of the reciprocal of the
%         condition number of the system of algebraic equations from
%         which the solution X is obtained, by applying the generalized
%         Schur method to the related extended pencil.
%         
%           [X,EV,RCOND1] = SLGCARE(A,E,Q,R,B,L,FLAG)  does not compute
%         the feedback gain F.
%         
%           [X1,X2,EV,RCOND1] = SLGCARE(A,E,Q,R,B,L,FLAG)  computes a
%                                               [ X1 ]
%         basis of the deflating subspace range [    ], the associated
%                                               [ X2 ]
%         eigenvalues EV and RCOND1 = 1/kappa_1(X1). (X = X2/X1 when X1
%         is nonsingular.)
%         
%           [X,EV,RCOND1] = SLGCARE(A,E,Q,G,FLAG)  computes the solution
%         X of the following equation, equivalent to (1),
%
%                      -1
%           0 = Q - L*R  *L' + A1'*X*E + E'*X*A1 - E'*X*G*X*E,
%
%                      -1                 -1
%         where G = B*R  *B', A1 = A - B*R  *L'.
%         
%           [X1,X2,EV,RCOND1] = SLGCARE(A,E,Q,G,FLAG)  computes a
%         basis of the deflating subspace.
%         
%         FLAG - (optional) vector containing options.
%                FLAG(1) = 0 : compute the stabilizing solution;
%                              otherwise, compute the anti-stabilizing
%                              solution.
%                FLAG(2) = 0 : compute both the solution X and the
%                              feedback gain matrix F;
%                        > 0 : compute the solution X only;
%                        < 0 : compute the deflating subspace.
%                FLAG(3)     : tolerance to check the singularity of
%                              the matrix pencil. If tol <= 0 it will
%                              be replaced by EPS, the machine epsilon.
%         Default: FLAG = [0,1,0.0].
%         If only five input arguments are specified, and all are
%         1-by-1, the fifth parameter is taken as B, not FLAG.
%         Similarly, if the first and forth parameters have equal size,
%         but the fifth is not a vector, it is assumed to be B. In both
%         cases, L is then taken as a zero matrix of a suitable size.
%         
%         See also GARESOL, SLGDARE
%         
 
%         RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%         Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%         V. Sima 11-03-2003.
%         Revised 03-03-2009.
%         
ni = nargin;
fl = [0,0,1,0]';
%
if ni < 4
   error('Usage: X = SLGCARE(A,E,Q,G,FLAG) or X = SLGCARE(A,E,Q,R,B,L,FLAG)')
elseif ni == 4
   withG = 1;
elseif ni == 5
   if numel( A ) ~= numel( R ) || min( size( B ) ) ~= 1,
      withG = 0;
   else
      withG = 1;
      flag  = B;
   end
else
   withG = 0;
end
%
if ni >= 7 || ( ni == 5 && withG )
   [mf,nf] = size(flag);
   if mf > 1 && nf > 1
      error('flag must be a vector')
   end
%
   m = max(mf,nf);
%
   if m >= 1 && flag(1) ~= 0
      fl(2) = 1;
   end
   if m >= 2
      if flag(2) == 0
         fl(3) = 0;
      elseif flag(2) < 0
         fl(3) = -1;
      end
   end
   if m >= 3
      fl(4) = flag(3);
   end
end
%
if withG == 1
   if fl(3) == 0
      [X,F,ev,db,rcond1] = garesol(A,E,Q,R,fl);
      ev = ev./db;
   elseif fl(3) > 0
      [X,ev,db,rcond1] = garesol(A,E,Q,R,fl);
      F = ev./db;
      ev = rcond1;
   else
      [X,scale,ev,db,rcond1] = garesol(A,E,Q,R,fl);
      n = length(A);
      F = scale*X(n+1:n+n,:);
      X = X(1:n,:);
      ev = ev./db;
   end
   %
else
   %
   if ni == 5,  fl = fl(:);  fl = [ fl; 0; 0; 0; 0];  L = zeros( size( B ) );  end
   if fl(3) == 0
      [X,F,ev,db,rcond1] = garesol(A,E,Q,R,B,L,fl);
      ev = ev./db;
   elseif fl(3) > 0
      [X,ev,db,rcond1] = garesol(A,E,Q,R,B,L,fl);
      F = ev./db;
      ev = rcond1;
   else
      [X,scale,ev,db,rcond1] = garesol(A,E,Q,R,B,L,fl);
      n = length(A);
      F = scale*X(n+1:n+n,:);
      X = X(1:n,:);
      ev = ev./db;
   end
end
%
% end slgcare
