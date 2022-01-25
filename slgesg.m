function [X,Y,dif] = slgesg(A,D,B,E,C,F,flag,trans)
%SLGESG  Solve generalized linear matrix equation pairs.
%
%        [X,Y] = SLGESG(A,D,B,E,C,F,FLAG,TRANS)  computes the unique
%        solutions (X,Y) of the generalized linear matrix equation pairs
%
%             (  A*X - Y*B = C,
%             <
%             (  D*X - Y*E = F,
%
%        if TRANS = 0, or the "transposed" equation pairs
%
%             (  A'*X + D'*Y = C,
%             <
%             (  X*B' + Y*E' = -F.
%
%        if TRANS <> 0. Default: TRANS = 0.
%
%        FLAG is a row or column vector characterizing the structure
%        of the matrix pairs:
%        FLAG(1) = 1 : matrix pair (A,D) is in generalized Schur form,
%                      otherwise, (A,D) is in general form.
%        FLAG(2) = 1 : matrix pair (B,E) is in generalized Schur form,
%                      otherwise, (B,E) is in general form.
%        Default: FLAG = [0,0], i.e., both pairs (A,D), (B,E) are in
%        general forms.
%
%        [X,Y,DIF] = SLGESG(A,D,B,E,C,F,FLAG,TRANS)  computes the
%        solution (X,Y), as well as Diff[(A,D),(B,E)].
%
%        See also GENLEQ, SLGELY, SLGEST, SLGSLY, SLGSST
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        H. Xu 7-12-1998.
%        Revised V. Sima 30-04-1999.
%

ni = nargin;
%
if ni < 7
   flag = zeros(2,1);
end
%
if ni < 8
   trans = 0;
end
%
if nargout <= 2
   [X,Y] = genleq(1,A,D,B,E,C,F,flag,trans);
else
   [X,Y,dif] = genleq(1,A,D,B,E,C,F,flag,trans);
end
%
% end slgesg
