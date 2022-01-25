function [syso,No,U,s] = slobsf(sys,tol)
%SLOBSF  Observability staircase form of a system SYS = (A,B,C).
%
%        [SYSO,NO,U,S] = SLOBSF(SYS,TOL)  computes the observable
%        staircase form corresponding to (A,C), using an orthogonal
%        transformation matrix U, and updates the matrix B. 
%        The transformed system is returned in the LTI object SYSO.
%        The size of the observable subsystem, and the block sizes in
%        the staircase form of the matrix A, are also returned in the
%        scalar NO, and in a vector S, respectively.
%
%        [SYSO,NO] = SLOBSF(SYS,TOL)  computes the transformed system
%        SYSO, without accumulating the orthogonal matrix U.
%
%        TOL is a lower bound on the reciprocal condition numbers, used
%        to decide matrix ranks in the reductions.
%        Default: TOL = N*N*machine_epsilon, where N = size(A,1).
%
%        This m function uses the LTI object of Control System Toolbox. 
%        If this toolbox is not available, the mex function SYSCOM 
%        can be directly called in the form 
%        [A,B,C,NO,U,S] = SYSCOM(2,A,B,C,TOL);
%
%        See also SLCONF, SLMINR, SYSCOM
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
   [A,B,C,No] = syscom(2,A,B,C,tol);
else
   [A,B,C,No,U,s] = syscom(2,A,B,C,tol);
end
%
syso = ss(A,B,C,D,sys);
%
% end slobsf
