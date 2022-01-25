function [syscf,Nc,U,s] = slconf(sys,tol)
%SLCONF  Controllability staircase form of a system SYS = (A,B,C).
%
%        [SYSCF,NC,U,S] = SLCONF(SYS,TOL)  computes the controllable
%        staircase form corresponding to (A,B), using an orthogonal
%        transformation matrix U, and updates the matrix C.
%        The transformed system is returned in the LTI object SYSCF.
%        The size of the controllable subsystem, and the block sizes in
%        the staircase form of the matrix A, are also returned in the
%        scalar NC, and in a vector S, respectively.
%
%        [SYSCF,NC] = SLCONF(SYS,TOL)  computes the transformed system
%        SYSCF, without accumulating the orthogonal matrix U.
%
%        TOL is a lower bound on the reciprocal condition numbers, used
%        to decide matrix ranks in the reductions.
%        Default: TOL = N*N*machine_epsilon, where N = size(A,1).
%
%        This m function uses the LTI object of Control System Toolbox. 
%        If this toolbox is not available, the mex function SYSCOM 
%        can be directly called in the form 
%        [A,B,C,NC,U,S] = SYSCOM(1,A,B,C,TOL);
%
%        See also SLMINR, SLOBSF, SYSCOM
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        H. Xu 8-12-1998.
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
[A,B,C,D] = ssdata(sys);
%
if nargout <= 2
   [A,B,C,Nc] = syscom(1,A,B,C,tol);
else
  [A,B,C,Nc,U,s] = syscom(1,A,B,C,tol);
end  
% 
syscf = ss(A,B,C,D,sys);
%
% end slconf
