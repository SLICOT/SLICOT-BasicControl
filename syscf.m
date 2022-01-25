%SYSCF  MEX-function based on SLICOT coprime factorization routines.
%
%   [Af,Bf,Cf,Df,NR] = SYSCF(METH,A,B,C,D,TOL,DISCR,ALPHA,BETA)
%
%   SYSCF computes for a state-space system (A,B,C,D) with the
%   transfer-function matrix G a left coprime factorization (LCF) or
%   a right coprime factorization (RCF)
%                   -1                    -1
%         LCF: G = M  N,  or  RCF:  G = NM  ,
%
%   respectively, where (Af,Bf,Cf,Df) is the state space system with
%   the transfer-function matrix [ N M ] for a LCF or [ N ] for a RCF.
%                                                     [ M ]
%
%   NR is the order of a minimal realization of M which can be recovered
%   from the resulting state-space matrices (see LCF, RCF, LCFID or RCFID).
%
%   Description of other input parameters:
%   METH  - type of right/left coprime factorization:
%           = 1 : RCF with inner denominator;
%           = 2 : LCF with inner denominator;
%           = 3 : RCF with ALPHA stability degree;
%           = 4 : LCF with ALPHA stability degree.
%   TOL   - (optional) controllability/observability tolerance
%            for computing right/left coprime factorizations.
%            Default controllability tolerance (for computing RCF)
%                TOL = epsilon_machine*max(norm(A),norm(B));
%            default  observability tolerance  (for computing LCF)
%                TOL = epsilon_machine*max(norm(A),norm(C)).
%   DISCR - (optional) type of system:
%              = 0 : continuous-time (default);
%              = 1 : discrete-time.
%   ALPHA - (optional) stability degree for LCF/RCF with prescribed
%            stability degree.
%            Default: -0.05  for continuous-time;
%                      0.95  for discrete-time.
%   BETA  - (optional) stability margin for the LCF/RCF with prescribed
%            stability degree.
%            Default: beta = alpha.


%   RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   A. Varga 2-24-99.
%   Revised: A. Varga, 20-9-2005.

%   Reference:
%   [1] A. Varga
%       Computation of coprime factorizations of rational matrices.
%       Linear Algebra and Its Applications, vol. 271, pp.83--115, 1998.

