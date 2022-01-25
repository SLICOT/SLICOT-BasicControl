% MUHOPT.F   - MEX-function to compute the mu optimal or H_inf
%              controller using SLICOT routines SB10AD, SB10MD, AB04MD,
%              AB05MD, and AB07ND.
%
%   If mu optimal controller is desired (job > 0):
%
%   [AK,BK,CK,DK(,mju,RCOND)] = muHopt(job,discr,A,B,C,D,ncon,nmeas,
%                gamma,omega,nblock,itype,ord(,qutol(,gtol(,actol))))
%
%   If H_inf optimal controller is only desired (job <= 0):
%
%   [AK,BK,CK,DK(,gammin,RCOND)] = muHopt(job,discr,A,B,C,D,ncon,
%                nmeas,gamma(,gtol(,actol)))
%
%   MUHOPT computes the matrices of the mu optimal or H_inf optimal
%   controller given the model in a state space. It also outputs the
%   mu norm of the closed loop system, if mu optimal controller is
%   desired, or the value of gamma reached in the H_inf synthesis
%   problem, if H_inf controller is only desired.
%   The discrete-time systems are handled via bilinear transformation
%   to continuous-time and then the controller obtained is discretised.
%   For the K step the SB10AD subroutine is employed, and the SB10MD
%   subroutine performs the D step.
%
%   Description of input parameters:
%   job    - indicates the type of the controller as well as the strategy
%            for reducing the gamma value:
%            >  0 : mu optimal controller is desired;
%            <= 0 : H_inf optimal controller only is desired.
%            Specifically, job
%            = 0 : find suboptimal controller only;
%            and abs(job) specifies the strategy for reducing gamma:
%            = 1 : use bisection method for decreasing gamma from gamma
%                  to gammin until the closed-loop system leaves
%                  stability;
%            = 2 : scan from gamma to 0 trying to find the minimal gamma
%                  for which the closed-loop system retains stability;
%            = 3 : first bisection, then scanning.
%   discr  - indicates the type of the system, as follows:
%            = 0 : continuous-time system;
%            = 1 : discrete-time system.
%   A      - the n-by-n system state matrix A of the plant.
%   B      - the n-by-m system input matrix B of the plant.
%   C      - the p-by-n system output matrix C of the plant.
%   D      - the p-by-m system input/output matrix D of the plant.
%   ncon   - the number of control inputs.
%            p-nmeas >= ncon >= 0.
%   nmeas  - the number of measurements.
%            p-nmeas = m-ncon >= nmeas >= 0.
%   gamma  - the initial value of gamma on input. It is assumed that
%            gamma is sufficiently large so that the controller is
%            admissible.  gamma >= 0.
%   omega  - the vector of length lendat >= 2 with the frequencies.
%            They must be nonnegative, in increasing order, and
%            for discrete-time systems between 0 and pi.
%   nblock - the vector with the block structure of the uncertainty.
%            nblock(I) is the size of each block.
%   itype  - the vector of the same length as nblock indicating
%            the type of each block.
%            itype(I) = 1 indicates that the corresponding block is a
%            real block. THIS OPTION IS NOT SUPPORTED NOW.
%            itype(I) = 2 indicates that the corresponding block is a
%            complex block. THIS IS THE ONLY ALLOWED VALUE NOW!
%   ord    - the maximum order of each block in the D-fitting procedure.
%            1 <= ord <= lendat-1.
%   qutol  - (optional) the acceptable mean relative error between
%            the D(jw) and the frequency response of the estimated block
%            [ADi,BDi;CDi,DDi]. When it is reached, the result is
%            taken as good enough.
%            Default: qutol = 2.
%   gtol   - (optional) tolerance used for controlling the accuracy
%            of gamma and its distance to the estimated minimal possible
%            value of gamma.
%            If gtol <= 0, then sqrt(EPS) is used, where EPS is the
%            relative machine precision.
%            Default: gtol = 0.01.
%   actol  - (optional) upper bound for the poles of the closed-loop
%            system used for determining if it is stable.
%            actol <= 0 for stable systems.
%            Default: actol = 0.
%
%   Description of output parameters:
%   AK     - the n-by-n controller state matrix AK.
%   BK     - the n-by-nmeas controller input matrix BK.
%   CK     - the ncon-by-n controller output matrix CK.
%   DK     - the ncon-by-nmeas controller input/output matrix DK.
%   mju    - (optional) the vector with the estimated upper bound of
%            the structured singular value for each frequency in omega
%            for the closed-loop system.
%   gammin - (optional) the estimated minimal admissible gamma.
%   RCOND  - (optional) for each successful J-th K step: 
%            RCOND(J) contains the reciprocal condition number of the
%                     control transformation matrix;
%            RCOND(J+1) contains the reciprocal condition number of the
%                     measurement transformation matrix;
%            RCOND(J+2) contains an estimate of the reciprocal condition
%                     number of the X-Riccati equation;
%            RCOND(J+3) contains an estimate of the reciprocal condition
%                     number of the Y-Riccati equation.
%            If job = 0, only RCOND(1:4) are set by the first K step.
%

% RELEASE 2.0 of SLICOT Robust Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   A. Markovski, Technical University of Sofia, October 2003.
%
% Revisions:
%   V. Sima, April 2004.
%
