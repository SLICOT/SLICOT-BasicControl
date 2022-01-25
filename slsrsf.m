function [syssch,U] = slsrsf(sys)
%SLSRSF  Transform the state matrix A to a real Schur form.
%
%        [SYSSCH,U] = SLSRSF(SYS)  applies an orthogonal similarity
%        transformation U on the system matrix S = [A, B; C, 0], of 
%        a state-space system SYS = (A,B,C), so that the transformed
%        state matrix is in a real Schur form.
%
%        SYSSCH = SLSRSF(SYS)  performs the same computations, but U is
%        not returned.
%
%        This m function uses the LTI object of Control System Toolbox. 
%        If this toolbox is not available, the mex function SYSTRA 
%        can be directly called in the form 
%        [A,B,C] = SYSTRA(2,A,B,C);
%
%        See also SLSBAL, SLSDEC, SLSORSF, SYSTRA
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        H. Xu 8-12-1998.
%        Revised V. Sima 30-04-1999, 08-06-1999.
%

%
[A,B,C,D] = ssdata(sys);
%
if nargout <= 1
   [A,B,C] = systra(2,A,B,C);
else  
   [A,B,C,ev,U] = systra(2,A,B,C);
end
%
syssch = ss(A,B,C,D,sys);
%
% end slsrsf
