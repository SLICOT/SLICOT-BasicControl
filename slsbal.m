function [sysb,U] = slsbal(sys,job,maxred)
%SLSBAL  Balance the system matrix for a state-space system SYS = (A,B,C).
%
%        [SYSB,U] = SLSBAL(SYS,JOB,MAXRED)  computes a diagonal 
%        similarity transformation U and applies it on the system
%        matrix S = [A, B; C, 0], so that [inv(U)*A*U, inv(U)*B; C*U, 0] 
%        has approximately equal row and column norms.
%
%        SYSB = SLSBAL(SYS,JOB,MAXRED)  performs the same computations,
%        but U is not returned.
%
%        JOB is an integer parameter specifying the part of the system
%        matrix used to compute the scaling factors for balancing:
%           JOB = 1  : use the matrix A only.
%           JOB = 2  : use the matrices A and B.
%           JOB = 3  : use the matrices A and C.
%           otherwise, use all matrices A, B, and C.
%        Default:  JOB = 0.
%
%        MAXRED is a real parameter specifying the maximum allowed
%        reduction in the 1-norm of S, if zero rows or columns are 
%        encountered.
%        Default: MAXRED = 10.
%
%        This m function uses the LTI object of Control System Toolbox. 
%        If this toolbox is not available, the mex function SYSTRA 
%        can be directly called in the form
%        [A,B,C,U] = SYSTRA(1,A,B,C,JOB,MAXRED);
%
%        See also SLSDEC, SLSORSF, SLSRSF, SYSTRA
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        H. Xu 8-12-1998.
%        Revised V. Sima 30-04-1999, 08-06-1999.
%

ni = nargin;
%
if ni <= 1 
   job = 0;
end
if ni <= 2
   maxred = 10.0;
end
%
[A,B,C,D] = ssdata(sys);
%
if nargout <= 1
   [A,B,C] = systra(1,A,B,C,job,maxred);
else  
   [A,B,C,u] = systra(1,A,B,C,job,maxred);
   U = diag(u);
end 
%
sysb = ss(A,B,C,D,sys);
%
% end slsbal
