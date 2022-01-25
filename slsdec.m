function [sys1,sys2] = slsdec(sys,sdeg)
%SLSDEC  Additive spectral decomposition of a system SYS = (A,B,C) with
%        respect to a given stability domain defined by SDEG.
%
%        [SYS1,SYS2] = SLSDEC(SYS,SDEG)  applies a similarity
%        transformation with a nonsingular matrix U on the system matrix
%        S = [A, B; C, 0], of a state-space system SYS = (A,B,C), 
%        so that the state matrix A is reduced to a block diagonal form
%        with eigenvalues of the top left diagonal block satisfying a
%        stability condition.
%
%        SDEG is a real number specifying the required stability margin 
%        for the eigenvalues to be moved to the top left diagonal block.
%        These eigenvalues will have real parts less than SDEG, for a
%        continuous-time system, and moduli less than SDEG, for a
%        discrete-time system.
%        Default: SDEG =    -sqrt(epsilon_machine) for continuous-time;
%                 SDEG = 1.0-sqrt(epsilon_machine) for discrete-time.
%
%        This m function uses the LTI object of Control System Toolbox. 
%        If this toolbox is not available, the mex function SYSTRA 
%        can be directly called in the form 
%        [A,B,C,NUM,EV,U] = SYSTRA(4,A,B,C,FLAG,SDEG);
%        Type HELP SYSTRA to see the definitions for the additional
%        arguments.
%
%        See also SLSBAL, SLSORSF, SLSRSF, SYSTRA
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        H. Xu 2-11-1998.
%        Revised V. Sima 30-04-1999, 08-06-1999.
%
%        Comments
%        In discrete-time case, sdeg must be nonnegative.
%

ni  = nargin;
%
n = size(sys.a,1);
discr = sys.ts > 0;
%
if ni < 2
   sdeg = -sqrt(eps);
   if discr
      sdeg = 1 + sdeg;
   end
end
% 
flag = [0;0];
if discr
   flag(1) = 1;
end
%
[A,B,C,D] = ssdata(sys);
%
[A,B,C,num] = systra(4,A,B,C,flag,sdeg);
%
sys1 = ss(A(1:num,1:num),B(1:num,:),C(:,1:num),D,sys);
sys2 = ss(A(num+1:n,num+1:n),B(num+1:n,:),C(:,num+1:n),0,sys);
%
% end slsdec
