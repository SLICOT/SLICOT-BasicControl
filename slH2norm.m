function [Hnrm,sysq] = slH2norm(sys,jobn,tol)
%SLH2NORM  H2/L2-norm of a system.
%
%        HNRM = SLH2NORM(SYS)  returns the H2-norm of a stable
%        system SYS = (A,B,C,D), i.e., sqrt(trace(B'*X*B)),
%        where X is the solution of a Lyapunov equation.
%        If SYS is a continuous-time system, D is assumed 0.
%        The system must not have poles on the boundary of the
%        stability domain (the imaginary axis, for a continuous-time
%        system, or the unit circle, for a discrete-time system).
%
%        HNRM = SLH2NORM(SYS,JOBN,TOL)  has additional input
%        arguments:
%
%        JOBN specifies the norm to be computed:
%        JOBN = 1 :  the H2-norm;
%        JOBN = 2 :  the L2-norm. In this case, SYS could be unstable.
%        Default: JOBN = 1.
%
%        TOL specifies an absolute tolerance level below which the
%        elements of B are considered zero (used for controllability
%        tests). If tol <= 0, the default value is used.
%        Default: TOL = n*norm(B,1)*eps.
%
%        [HNRM,SYSQ] = SLH2NORM(SYS,JOBN)  also returns a state-space
%        representation (AQ,BQ,CQ,DQ) of the left factor Q of a
%        right coprime factorization with inner denominator of the
%        transfer-function G of the system, G = Q*inv(R), where Q
%        and R are stable transfer-function matrices and R is inner.
%        If G is stable, then Q = G and R = I. The returned norm is
%        sqrt(trace(BQ'*X*BQ)),  for a continuous-time system, and
%        sqrt(trace(BQ'*X*BQ + DQ'*DQ)), for a discrete-time system.
%
%        See also SLHKNORM, SLINORM
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, Dec. 2002.
%
%        Revisions:
%        V. Sima, April 2003, Jan. 2007.
%

if nargin < 1,
    error('Usage: HNRM = SLH2NORM(SYS,JOBN,TOL)')
elseif nargin < 2,
   jobn = 1;
   tol  = 0;
elseif nargin < 3,
   tol = 0;
end

Ts = sys.Ts;
if Ts ~= 0
   dico = 2;
else
   dico = 1;
end
[A,B,C,D] = ssdata(sys);

[Hnrm,nq,Ao,Bo,Co,Do] = Hnorm(2,A,B,C,D,dico,jobn,tol);
if nargout > 1,
    sysq = ss(Ao,Bo,Co,Do,Ts);
end
%
% end slH2norm
