function [sysosch,nd,U] = slsorsf(sys,stab,sdeg)
%SLSORSF Transform the state matrix of a state space system
%        to a specified eigenvalue-ordered real Schur form.
%
%        [SYSOSCH,ND,U] = SLSORSF(SYS,STAB,SDEG)  performs a similarity
%        transformation with an orthogonal matrix U on the system 
%        SYS = (A,B,C,D) such that the state matrix of the transformed
%        system SYSOSCH = (U'*A*U,U'*B,C*U,D) is reduced to a real Schur 
%        form with the eigenvalues of the leading ND x ND diagonal
%        block satisfying the stability condition specified by the
%        parameters STAB and SDEG.
%
%        STAB is a scalar specifying the stability condition put
%        on the eigenvalues of the leading block:
%        STAB = 0 : the stable eigenvalues of the state matrix U'*A*U 
%                   should appear in the leading diagonal block;
%                   otherwise, the unstable eigenvalues should appear
%                   in the leading diagonal block.
%        Default: STAB = 0.
%
%        SDEG is a real number specifying the stability or instability 
%        degree for the eigenvalues of the leading diagonal block.
%        Default: SDEG =    -sqrt(epsilon_machine) for continuous-time;
%                 SDEG = 1.0-sqrt(epsilon_machine) for discrete-time.
%
%        This m function uses the LTI object of Control System Toolbox. 
%        If this toolbox is not available, the mex function SYSTRA 
%        can be directly called in the form  
%        [A,B,C,ND,EV,U] = SYSTRA(3,A,B,C,STAB,SDEG);        
%        Type HELP SYSTRA to see the definitions for the additional
%        arguments.
%
%        See also SLSBAL, SLSDEC, SLSRSF, SYSTRA
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
no  = nargout;
%
discr = sys.ts > 0;
% 
if discr
   fl(1) = 1;
else
   fl(1) = 0;
end
%
if ni < 2
   fl(2) = 0;
else
   fl(2) = stab;
end
%
if ni < 3
   sdeg = -sqrt(eps);
   if discr
      sdeg = 1 + sdeg;
   end
end
%
[A,B,C,D] = ssdata(sys);
%
if no >= 3
   [A,B,C,nd,ev,U] = systra(3,A,B,C,fl,sdeg);
else
   [A,B,C,nd] = systra(3,A,B,C,fl,sdeg);
end
%
sysosch = ss(A,B,C,D,sys);
%
% end slsorsf
