function [sysn,sysm] = rcfid(sys,tol)
%RCFID  Right coprime factorization with inner denominator.
%       [SYSN,SYSM] = RCFID(SYS,TOL)  calculates for the transfer
%       function matrix
%                                       -1
%               G(lambda) =  C(lambdaI-A) B + D
%
%       of a given system SYS = (A,B,C,D) a right coprime factorization
%                                              -1
%               G(lambda) = N(lambda)*M(lambda)
% 
%       where N(lambda) and M(lambda) are the transfer function 
%       matrices of two stable systems and M(lambda) is inner.
%       TOL is an optional controllability tolerance
%       Default value: TOL = epsilon_machine*max(norm(A),norm(B)) 

%       RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%       Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%       Interface M-function to the SLICOT-based MEX-function SYSCF.
%       A. Varga 22-07-1998. Revised 25-02-1999.
%       Revised, V. Sima 12-01-2002.
%

if ~isa(sys,'lti')
   error('The input system SYS must be an LTI object')
end

ni = nargin;
discr = double(sys.ts > 0);
if ni < 2 
   tol = 0; 
end

[a,b,c,d,TS]=ssdata(sys);

[ar,br,cr,dr,nr]=syscf(1,a,b,c,d,tol,discr);
[n,m] = size(br); p = size(cr,1);
i1 = 1:p-m; i2 = n-nr+1:n; i3 =p-m+1:p; 

sysn = ss(ar,br,cr(i1,:),dr(i1,:),TS);
sysm = ss(ar(i2,i2),br(i2,:),cr(i3,i2),dr(i3,:),TS);


% end rcfid
