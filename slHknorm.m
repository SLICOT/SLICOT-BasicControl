function [Hnrm,syso,hsv] = slHknorm(sys,equil,alpha)
%SLHKNORM  Hankel-norm of a system.
%
%        HNRM = SLHKNORM(SYS)  returns the Hankel-norm of the
%        stable projection of the transfer-function matrix
%        of the system SYS = (A,B,C,0).
%
%        HNRM = SLHKNORM(SYS,EQUIL,ALPHA)  has additional input
%        arguments:
%
%        EQUIL specifies whether the system (A,B,C) should be
%        preliminarily equilibrated:
%        EQUIL = 1 :  do not perform equilibration;
%        EQUIL = 2 :  perform equilibration (scaling).
%        Default: EQUIL = 1.
%
%        ALPHA specifies the ALPHA-stability boundary for the
%        eigenvalues of the state dynamics matrix A. For a
%        continuous-time system, ALPHA <= 0 is the boundary value
%        for the real parts of eigenvalues, while for a discrete-
%        time system, 0 <= ALPHA <= 1 is the boundary value for the
%        moduli of eigenvalues. The ALPHA-stability domain is defined
%        as the open half complex plane left to ALPHA, or the
%        interior of the ALPHA-radius circle centered in the origin,
%        for a continuous-time or discrete-time system, respectively.
%        Default: ALPHA =  -sqrt(eps), for a continuous-time system;
%                 ALPHA = 1-sqrt(eps), for a discrete-time system.
%
%        [HNRM,SYSO,HSV] = SLHKNORM(SYS,EQUIL,ALPHA)  also returns a state-
%        space representation SYSO = (AO,BO,CO,DO) of a transformed system,
%        whose state dynamics matrix AO is in a block diagonal real Schur
%        form with its eigenvalues reordered and separated, as well as the 
%        ns-vector HSV containing the Hankel singular values (ordered
%        decreasingly) of the ALPHA-stable part of the system. AO has two
%        diagonal blocks: the leading ns-by-ns part has eigenvalues in 
%        the ALPHA-stability domain and the trailing (n-ns)-by-(n-ns)
%        part has eigenvalues outside the ALPHA-stability domain.
%
%        See also SLH2NORM, SLINORM
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, Dec. 2002.
%
%        Revisions:
%        V. Sima, April 2003, Jan. 2007, Mar. 2009.
%

if nargin < 1,
    error('Usage: HNRM = SLHKNORM(SYS,EQUIL,ALPHA)')
elseif nargin < 2,
   equil = 1;
   alpha = 0;
elseif nargin < 3,
   alpha = 0;
end

Ts = sys.ts;
if Ts ~= 0
   dico = 2;
   if nargin < 3,  alpha = 1;  end
else
   dico = 1;
end

[A,B,C] = ssdata(sys);

[Hnrm,ns,Ao,Bo,Co,hsv] = Hnorm(1,A,B,C,dico,equil,alpha);
if nargout > 1,
    D = zeros( size( Co,1 ), size( Bo,2 ) );  syso = ss(Ao,Bo,Co,D,Ts);
end
%
% end slHknorm
