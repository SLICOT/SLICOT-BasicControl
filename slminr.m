function sysr = slminr(sys,tol,bal)
%SLMINR  Minimal realization of a system SYS = (A,B,C).
%
%        SYSR = SLMINR(SYS,TOL,BAL)  computes a minimal realization
%        SYSR = (AR,BR,CR) of an original state-space system SYS = (A,B,C).
%
%        TOL is a lower bound on the reciprocal condition numbers, used
%        to decide matrix ranks in the reductions.
%        Default: TOL = N*N*machine_epsilon, where N = size(A,1).
%
%        BAL is an integer specifying whether or not balancing should be
%        performed before reductions. If BAL = 0, balancing is performed,
%        otherwise, balancing is not performed.
%        Default: BAL = 0.
%
%        This m function uses the LTI object of Control System Toolbox. 
%        If this toolbox is not available, the mex function SYSCOM 
%        can be directly called in the form  
%        [AR,BR,CR] = SYSCOM(3,A,B,C,TOL,BAL);
%
%        See also SLCONF, SLOBSF, SLSBAL, SYSCOM, SYSTRA
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        H. Xu 8-12-98.
%        Revised V. Sima 30-04-1999, 08-06-1999.
%

ni = nargin;
%
% Default tol = N*N*machine_epsilon is set by SLICOT routine.
%
if ni <= 1
   tol = 0.0;
end
%
if ni <= 2
   bal = 0;
end
%
[A,B,C,D] = ssdata(sys);
%
[Ar,Br,Cr] = syscom(3,A,B,C,tol,bal);
%
if size(Ar,1) < size(A,1)
   sysr = ss(Ar,Br,Cr,D,sys);
else
   sysr = sys;
end
%
% end slminr
