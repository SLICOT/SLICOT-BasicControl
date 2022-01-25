function [sysn,sysm] = lcf(sys,alpha,tol)
%LCF    Left coprime factorization with prescribed stability degree.
%       [SYSN,SYSM] = LCF(SYS,ALPHA,TOL) calculates for the transfer
%       function matrix
%                                       -1
%               G(lambda) =  C(lambdaI-A) B + D
%
%       of a given system SYS = (A,B,C,D) a left coprime factorization
%                                    -1
%               G(lambda) = M(lambda)  *N(lambda)
% 
%       where N(lambda) and M(lambda) are the transfer function 
%       matrices of two stable systems with stability degree ALPHA.
%       Default values: ALPHA = -0.05  for a continuous-time system and
%       ALPHA = 0.95  for a discrete-time system.
%       TOL is an optional observability tolerance
%       Default value: TOL = epsilon_machine*max(norm(A),norm(C)) 

%       RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%       Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%       Interface M-function to the SLICOT-based MEX-function SYSCF.
%       A. Varga 28-07-1998. Revised 24-02-1999.
%       Revised, V. Sima 12-01-2002.
%

if ~isa(sys,'lti')
   error('The input system SYS must be an LTI object')
end
 
ni = nargin;
discr = double(sys.ts > 0);
if ni < 2
   if discr
      alpha = 0.95;
   else
      alpha = -0.05;
   end
end
if ni < 3 
   tol = 0; 
end

[a,b,c,d,TS]=ssdata(sys);

[ar,br,cr,dr,nr]=syscf(4,a,b,c,d,tol,discr,alpha);
[n,m] = size(br); p = size(cr,1);
j1 = 1:m-p; i2 = 1:nr; j3 =m-p+1:m; 

sysn = ss(ar,br(:,j1),cr,dr(:,j1),TS);
sysm = ss(ar(i2,i2),br(i2,j3),cr(:,i2),dr(:,j3),TS);


% end lcf
