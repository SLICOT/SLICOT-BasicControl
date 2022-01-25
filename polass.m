% POLASS.F - MEX-function for performing (partial) pole assignment using
%            SLICOT routine SB01BD.
%
%   [F(,split,WRo,WIo,Z,Ao)] = polass(A,B,WR,WI(,tol,discr,alpha))
%
%  POLASS determines the state feedback matrix F for a given system (A,B)
%  such that the closed-loop state matrix A+B*F has specified eigenvalues
%  (closed-loop system poles).
%
%   Description of input parameters:
%   A     - the n-by-n system state matrix A.
%   B     - the n-by-m system input matrix B.
%   WR    - the np-vector of real parts of the desired system poles, 
%           np <= n.
%   WI    - the np-vector of imaginary parts of the desired system
%           poles. Complex conjugate pairs must appear consecutively.
%   tol   - (optional) absolute tolerance level below which the elements
%           of A or B are considered zero (used for controllability
%           tests). If tol <= 0, then a default tolerance is used.
%           Default:    n*epsilon_machine*max(norm(A),norm(B)), where 
%                       epsilon_machine is the relative machine 
%                       precision, and norm denotes the 1-norm.
%   discr - (optional) scalar indicating the type of system:
%           = 0 : continuous-time (default);
%           = 1 : discrete-time.
%   alpha - (optional) scalar specifying the maximum admissible value,
%           either for real parts, if discr = 0, or for moduli,
%           if discr = 1, of the eigenvalues of A which will not be
%           modified by the eigenvalue assignment algorithm
%           (alpha >= 0 if discr = 1).
%           Default:    -sqrt(epsilon_machine)  for continuous-time;
%                    1.0-sqrt(epsilon_machine)  for discrete-time.
%
%   Description of output parameters:
%   F     - the m-by-n state feedback matrix, which assigns nap
%           closed-loop poles and keeps unaltered n-nap open-loop poles.
%   split - (optional) a 3-vector containing the pole splitting details;
%           split = [nfp; nap; nup], where
%                    nfp - the number of fixed poles, not modified by
%                          the assignment algorithm (poles having real
%                          parts, if discr = 0, or moduli, if discr = 1,
%                          less than alpha);
%                    nap - the number of assigned poles,
%                          nap = n - nfp - nup;
%                    nup - the number of uncontrollable poles detected
%                          by the algorithm.
%   WRo,  - (optional) the np-vector of real and imaginary parts,
%   WIo     respectively, of the closed-loop poles of interest;
%           the leading nap elements contain the real/imaginary parts 
%           of the assigned poles, and the trailing np-nap elements
%           contain the unassigned poles. 
%   Z     - (optional) the n-by-n orthogonal matrix which reduces the 
%           closed-loop system state matrix A+B*F to real Schur form.
%   Ao    - (optional) the n-by-n matrix Z'*(A+B*F)*Z in a real Schur
%           form. The leading nfp-by-nfp diagonal block of Ao
%           corresponds to the fixed poles. The trailing nup-by-nup 
%           diagonal block of A corresponds to the uncontrollable poles.
%
%   Comments
%   1. Not all uncontrollable poles of the pair (A,B) are necessarily
%      detected by the eigenvalue assignment algorithm. Undetected 
%      uncontrollable poles may exist if nfp > 0 and/or np < n-nfp.
%   2. To assign all controllable poles, set alpha = -Inf, for discr = 0,
%      and alpha = 0, for discr = 1.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, June 2002.
%
% Revisions:
%   V. Sima, Research Institute for Informatics, Bucharest, July 2002.
%
