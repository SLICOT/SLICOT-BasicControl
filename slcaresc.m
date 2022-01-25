function  [X,F,ev,rcond1,acc] = slcaresc(A,Q,R,B,L,flag)
%SLCARESC Solve CARE with refined Schur method and estimate condition.
%
%         Case I.
%           [X,F,EV,RCOND1] = SLCARESC(A,Q,R,B,L,FLAG)  computes the
%         symmetric solution X of a continuous-time algebraic Riccati
%         equation (ARE)
%                                             -1
%             0 = Q + A'*X + X*A - (L + X*B)*R  *(L + X*B)' ,    (1)
%         
%         the feedback gain matrix F, the associated eigenvalues EV,
%         as well as RCOND1, an estimate of the reciprocal of the
%         condition number of the system of algebraic equations from
%         which the solution X is obtained, by using the Schur method
%         on the related Hamiltonian matrix.
%         
%           [X,EV,RCOND1,ACC] = SLCARESC(A,Q,R,B,L,FLAG)  does not compute
%         the feedback gain F, but computes the vector ACC, containing
%         estimates of the separation, reciprocal condition number and
%         forward error bound.
%         
%           [X1,X2,EV,RCOND1] = SLCARESC(A,Q,R,B,L,FLAG)  computes a
%                                               [ X1 ]
%         basis of the invariant subspace range [    ], the associated
%                                               [ X2 ]
%         eigenvalues EV and RCOND1 = 1/kappa_1(X1). (X = X2/X1 when
%         X1 is nonsingular.)
%         
%         Case II.
%           [X,EV,RCOND1,ACC] = SLCARESC(A,Q,G,FLAG)  computes the
%         symmetric solution X of the ARE
%         
%             0 = Q + A'*X + X*A - X*G*X,                        (2)
%         
%         the associated eigenvalues EV, as well as RCOND1 and ACC.
%         
%           [X1,X2,EV,RCOND1] = SLCARESC(A,Q,G,FLAG)  computes a basis
%                                         [ X1 ]
%         of the invariant subspace range [    ], the associated
%                                         [ X2 ]
%         eigenvalues EV and RCOND1.
%         
%         FLAG - (optional) vector containing options.
%                FLAG(1) = 0 : compute the stabilizing solution;
%                              otherwise, compute the anti-stabilizing
%                              solution.
%                FLAG(2) = 0 : use a scaling strategy;
%                              otherwise, do not use a scaling strategy.
%                FLAG(3) = 0 : compute both the solution X and the
%                              feedback gain matrix F;
%                        = 1 : compute the solution X only;
%                        = 2 : compute the solution X and estimates of
%                              the separation, reciprocal condition
%                              number and forward error bound.
%                        < 0 : compute the invariant subspace.
%         Default: FLAG = [0,0,1].
%         
%         See also ARESOLC, SLCAREGS, SLDAREGS, SLDAREGSV, SLDARESC
%         

%         RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%         Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%         V. Sima 11-08-1999.
%         Revised V. Sima 19-01-2003, 03-03-2009.
%         
%         Comments:
%         1. Matrix R must be nonsingular.
%         2. In case II the number of input parameters is at most 4.
%         
ni = nargin;
fl = [0,0,0,1]';
%
if ni == 4
   flag = B;
end
%
if ni >= 6 || ni == 4
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
   if m >= 2 && flag(2) ~= 0
      fl(3) = 1;
   end
   if m >= 3
      if flag(3) == 0
         fl(4) = 0;
      elseif flag(3) < 0
         fl(4) = -1;
      elseif flag(3) > 0
         fl(4) = flag(3);
      end
   end
end
%
if ni >= 5
   if fl(4) == 0
      [X,F,ev,rcond1] = aresolc(1,A,Q,R,B,L,fl);
   elseif fl(4) == 2
      [X,ev,rcond1,acc] = aresolc(1,A,Q,R,B,L,fl);
   elseif fl(4) < 0
      [X,scale,ev,rcond1] = aresolc(1,A,Q,R,B,L,fl);
   else
      [X,ev,rcond1] = aresolc(1,A,Q,R,B,L,fl);
   end
else
   if fl(4) == 2
      [X,ev,rcond1,acc] = aresolc(1,A,Q,R,fl);
   elseif fl(4) < 0
      [X,scale,ev,rcond1] = aresolc(1,A,Q,R,fl);
   else
      [X,ev,rcond1] = aresolc(1,A,Q,R,fl);
   end
end
%
if fl(4) > 0 || ( fl(4) == 0 && ni <= 4 )
   F = ev;
   ev = rcond1;
   if fl(4) == 2
      rcond1 = acc;
   end
elseif fl(4) < 0
   n = length(A);
   F = scale*X(n+1:n+n,:);
   X = X(1:n,:);
end
%
% end slcaresc
