function  [X,F,ev,rcond1] = slcaregs(A,Q,R,B,L,flag)
%SLCAREGS Solve CARE with generalized Schur method on an extended matrix
%         pencil.
%
%           [X,F,EV,RCOND1] = SLCAREGS(A,Q,R,B,L,FLAG)  computes the
%         symmetric solution X of the continuous-time algebraic Riccati
%         equation (ARE)
%                                              -1
%              0 = Q + A'*X + X*A - (L + X*B)*R  *(L + X*B)' ,      (1)
%         
%         the feedback gain matrix F, the associated eigenvalues EV,
%         as well as RCOND1, an estimate of the reciprocal of the
%         condition number of the system of algebraic equations from
%         which the solution X is obtained, by applying the generalized
%         Schur method to the related extended pencil.
%         
%           [X,EV,RCOND1] = SLCAREGS(A,Q,R,B,L,FLAG)  does not compute
%         the feedback gain F.
%         
%           [X1,X2,EV,RCOND1] = SLCAREGS(A,Q,R,B,L,FLAG)  computes a
%                                               [ X1 ]
%         basis of the deflating subspace range [    ], the associated
%                                               [ X2 ]
%         eigenvalues EV and RCOND1 = 1/kappa_1(X1). (X = X2/X1 when X1
%         is nonsingular.)
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
%         
%         See also ARESOL, SLCARES, SLDAREGS, SLDAREGSV, SLDARES
%         

%         RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%         Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%         H. Xu 7-12-1998.
%         Revised V. Sima 30-04-1999, 08-06-1999, 11-01-2003, 03-03-2009.
%         
ni = nargin;
fl = [0,0,1,0]';
%
if ni >= 6
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
if fl(3) == 0
   [X,F,ev,db,rcond1] = aresol(3,A,Q,R,B,L,fl);
   ev = ev./db;
elseif fl(3) > 0
   [X,ev,db,rcond1] = aresol(3,A,Q,R,B,L,fl);
   F = ev./db;
   ev = rcond1;
else
   [X,scale,ev,db,rcond1] = aresol(3,A,Q,R,B,L,fl);
   n = length(A);
   F = scale*X(n+1:n+n,:);
   X = X(1:n,:);
   ev = ev./db;
end
%
% end slcaregs
