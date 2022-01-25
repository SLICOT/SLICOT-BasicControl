function [X] = slsylv(A,B,C,flag,trans,schur)
%SLSYLV  Solve continuous-time Sylvester equations.
%
%        X = SLSYLV(A,B,C,FLAG,TRANS,SCHUR)  computes the unique
%        solution X of a continuous-time Sylvester equation
%
%             op(A)*X + X*op(B) = C,
%
%        where op(M) = M or M'.
%
%        FLAG is a vector characterizing the structures of A and B:
%        FLAG(1) =  1 : A is quasi upper triangular;
%                   2 : A is upper Hessenberg;
%                   otherwise, A is in general form.
%        FLAG(2) =  1 : B is quasi upper triangular;
%                   2 : B is upper Hessenberg;
%                   otherwise, B is in general form.
%        Default: FLAG = [0,0].
%
%        TRANS specifies if op(M) = M or M':
%        TRANS = 0 : op(A) = A;   op(B) = B;
%                1 : op(A) = A';  op(B) = B';
%                2 : op(A) = A';  op(B) = B;
%                3 : op(A) = A;   op(B) = B'.
%        Default: TRANS = 0.
%
%        SCHUR specifies whether the Hessenberg-Schur or Schur method
%              should be used:
%        SCHUR = 1 : Hessenberg-Schur method (one matrix is reduced
%                    to Schur form);
%        SCHUR = 2 : Schur method (two matrices are reduced to Schur
%                    form).
%        If one or both matrices are already reduced to Schur/Hessenberg
%        forms, this could be specified by FLAG(1) and FLAG(2).
%        For general matrices, the Hessenberg-Schur method could be 
%        significantly more efficient than the Schur method.
%        Default: SCHUR = 1.
%
%        See also LINMEQ, SLDSYL, SLLYAP, SLSTEI, SLSTLY, SLSTST
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        H. Xu 7-12-1998.
%        Revised V. Sima 30-04-1999.
%                        29-04-2000.
%                        12-09-2000.
%                        23-02-2009.
%

ni = nargin;
fl = [0;0;0];
%
if ni >= 4
   [mf,nf] = size(flag);
   if mf > 1 && nf > 1
      error('flag must be a vector.')
   end
   %
   m = max(mf,nf);
   %
   if m >= 1
      fl(2) = flag(1);
   end
   if m > 1
      fl(3) = flag(2);
   end
end
%
if ni < 5
   trans = 0;
end
%
if ni < 6
   schur = 1;
end
%
X = linmeq(1,A,B,C,fl,trans,schur);
%
% end slsylv.m
