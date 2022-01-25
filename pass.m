function [F,split,Po] = pass(A,B,P,tol,discr,alpha)
%PASS    Assigns part or all of the controllable poles of a system. 
% 
%        F = PASS(A,B,P)  computes a state feedback matrix F for the 
%        system (A,B) such that the closed-loop state matrix A+B*F has
%        specified eigenvalues P (closed-loop system poles).
%        Complex conjugate pairs must appear consecutively in P.
%
%        F = PASS(A,B,P,TOL,DISCR,ALPHA)  has additional input arguments:
%
%        TOL is an absolute tolerance level below which the elements
%        of A or B are considered zero (used for controllability tests).
%        If TOL <= 0, then a default tolerance is used.
%        Default:    n*epsilon_machine*max(norm(A),norm(B)), where 
%                    epsilon_machine is the relative machine precision, 
%                    and norm denotes the 1-norm.
%
%        DISCR indicates the type of system:
%        DISCR = 0 : continuous-time (default);
%        DISCR = 1 : discrete-time.
%
%        ALPHA is a scalar specifying the maximum admissible value,
%        either for real parts, if DISCR = 0, or for moduli, if DISCR = 1,
%        of the eigenvalues of A which will not be modified by the
%        eigenvalue assignment algorithm (ALPHA >= 0 if DISCR = 1).
%        Default:    -sqrt(epsilon_machine)  for continuous-time;
%                 1.0-sqrt(epsilon_machine)  for discrete-time.
%        To assign all controllable poles, set ALPHA = -Inf, if DISCR = 0, 
%        or ALPHA = 0, if DISCR = 1.
%
%        [F,SPLIT,PO] = PASS(A,B,P,TOL,DISCR,ALPHA)  also returns details 
%        on the assigned poles:
%
%        SPLIT = [nfp; nap; nup], where
%                 nfp - the number of fixed poles, not modified by the
%                       algorithm (poles having real parts, if DISCR = 0,
%                       or moduli, if DISCR = 1, less than ALPHA);
%                 nap - the number of assigned poles, nap = n - nfp - nup;
%                 nup - the number of uncontrollable poles detected
%                       by the algorithm.
%
%        PO is a vector with the closed-loop poles of interest: 
%        the leading nap elements contain the assigned poles, and the
%        trailing np-nap elements contain the unassigned poles. 
%
%        See also POLASS
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, June 2002.
%        Revised: V. Sima, July 2002, March 2009.
%

ni   = nargin;
nout = nargout;
%
if ni < 3,
    error('Usage: F = pass(A,B,P)')
elseif ni == 3,
    tol   = -1;
    discr = 0;
    alpha = -sqrt( eps );
elseif ni == 4,
    discr = 0;
    alpha = -sqrt( eps );
elseif ni == 5,
    alpha = -sqrt( eps );
    if discr == 1,  alpha = 1 + alpha;  end
end
%
if nout == 1,
    F = polass(A,B,real(P),imag(P),tol,discr,alpha);
elseif nout == 2,
    [F,split] = polass(A,B,real(P),imag(P),tol,discr,alpha);
elseif nout == 3,
    [F,split,Por,Poi] = polass(A,B,real(P),imag(P),tol,discr,alpha);  Po = Por + 1i*Poi;
end
%
% end pass
