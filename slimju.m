function [sysk,mju,rcnd] = slimju(sys,ncon,nmeas,gmax,w,nblock,itype,...
                                  maxord,job,discr,qutol,gtol,actol)
%SLIMJU  Compute the mu optimal controller given a state space
%        model. Also, compute the mu norm of the closed loop system.
%
%        SYSK = SLIMJU(SYS,NCON,NMEAS,GMAX,W,NBLOCK,ITYPE,MAXORD)  
%        computes the mu optimal controller SYSK, containing the matrices
%        Ak,Bk,Ck,Dk in a packed form, for a continuous-time system SYS
%        (given in packed form), with NCON control inputs and NMEAS
%        measurements, for the initial value GMAX of GAMMA on input.
%        It is assumed that GMAX is sufficiently large so that the
%        controller is admissible. GMAX >= 0. The other parameters are
%        described below.
%
%        W is a vector of frequencies, which must be nonnegative, and
%        in increasing order.
%
%        NBLOCK is a vector with the block structure of the uncertainty.
%        NBLOCK(i) is the size of each block.
%
%        ITYPE is a vector of the same length as NBLOCK indicating 
%        the type of each block.
%        ITYPE(i) = 1 indicates that the corresponding block is a
%        real block. THIS OPTION IS NOT SUPPORTED NOW.
%        ITYPE(i) = 2 indicates that the corresponding block is a 
%        complex block. THIS IS THE ONLY ALLOWED VALUE NOW!
%
%        MAXORD is the maximum order of each block in the D-fitting 
%        procedure.
%
%        [SYSK,MJU,RCND] = SLIMJU(SYS,NCON,NMEAS,GMAX,W,NBLOCK,ITYPE,...
%                                 MAXORD,JOB,DISCR,QUTOL,GTOL,ACTOL)
%        has additional input and output arguments:
%
%        JOB indicates the strategy for reducing the GAMMA value:
%        JOB = 1 : Use bisection method for decreasing GAMMA from
%                  GMAX to GMIN until the closed-loop system leaves
%                  stability.
%        JOB = 2 : Scan from GMAX to 0 trying to find the minimal GAMMA
%                  for which the closed-loop system retains stability.
%        JOB = 3 : First bisection, then scanning.
%        Default:  JOB = 3.
%
%        DISCR indicates the type of system:
%        DISCR = 0 : continuous-time (default);
%        DISCR = 1 : discrete-time.
%
%        QUTOL is the acceptable mean relative error between D(jw) and
%        the frequency response of the estimated block [ADi,BDi; CDi,DDi].
%        When it is reached, the result is taken as good enough.
%        Default:  QUTOL = 2.
%
%        GTOL is a tolerance used for controlling the accuracy of GMIN
%        and its distance to the estimated minimal possible value
%        of GAMMA. If GTOL <= 0, then sqrt(eps) is used, where eps is 
%        the relative machine precision.
%        Default:  GTOL = 0.01.
%
%        ACTOL is an upper bound for the poles of the closed-loop 
%        system used for determining if it is stable.
%        ACTOL <= 0 for stable systems.
%        Default:  ACTOL = 0.
%                                                                     
%        MJU is the vector with the estimated upper bound of the structured 
%        singular value for each W for the closed-loop system.
%
%        RCND contains the reciprocal condition numbers. Specifically, 
%        for each successful j-th K step: 
%        RCND(j) contains the reciprocal condition number of the
%                  control transformation matrix;
%        RCND(j+1) contains the reciprocal condition number of the
%                  measurement transformation matrix;
%        RCND(j+2) contains an estimate of the reciprocal condition
%                  number of the X-Riccati equation;
%        RCND(j+3) contains an estimate of the reciprocal condition
%                  number of the Y-Riccati equation.
%
%   Comments:
%   1) The discrete-time systems are handled via bilinear transformation
%   to continuous-time and then the controller obtained is discretized.
%   For discrete-time systems the frequencies must be between 0 and pi.
%   2) It is assumed that 
%      NP - NMEAS >= NCON  >= 0 and 
%      M  - NCON  >= NMEAS >= 0
%   where M and NP are the number of inputs and outputs, respectively,
%   of the system SYS.
%   3) The commands pck and unpck from the Mu-Analysis and Synthesis 
%   Toolbox are used.
% 
%   See also MUHOPT, SLIHINF

% RELEASE 2.0 of SLICOT Robust Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   A. Markovski, Technical University of Sofia, October 2003.
%
% Revisions:
%   V. Sima, Research Institute for Informatics, Bucharest, March 2004,
%            March 2009.
%

ni   = nargin;
nout = nargout;
%
if ( ni < 8 || ni > 13 || nout < 1 || nout > 3 )
   error( [ 'Usage: [SYSK]          = SLIMJU(SYS,NCON,NMEAS,GMAX,W,NBLOCK,ITYPE,MAXORD)',  sprintf('\n'),...
            '       [SYSK,MJU]      = SLIMJU(SYS,NCON,NMEAS,GMAX,W,NBLOCK,ITYPE,MAXORD,JOB,DISCR)',  sprintf('\n'),...
            '       [SYSK,MJU,RCND] = SLIMJU(SYS,NCON,NMEAS,GMAX,W,NBLOCK,ITYPE,MAXORD,JOB,DISCR,QUTOL,GTOL,ACTOL)' ] )
end
%
if ni == 8,
   job   = 3;
   discr = 0;
elseif ni == 9,
   discr = 0;
end
%
[a,b,c,d] = unpck( sys );
%
if ni == 10
   [ak,bk,ck,dk,mju,rcnd] = muHopt(job,discr,a,b,c,d,ncon,nmeas,gmax,w,nblock,itype,maxord);
elseif ni == 11
   [ak,bk,ck,dk,mju,rcnd] = muHopt(job,discr,a,b,c,d,ncon,nmeas,gmax,w,nblock,itype,maxord,qutol);
elseif ni == 12
   [ak,bk,ck,dk,mju,rcnd] = muHopt(job,discr,a,b,c,d,ncon,nmeas,gmax,w,nblock,itype,maxord,qutol,gtol);
elseif ni == 13
   [ak,bk,ck,dk,mju,rcnd] = muHopt(job,discr,a,b,c,d,ncon,nmeas,gmax,w,nblock,itype,maxord,qutol,gtol,actol);
end
%
sysk = pck(ak,bk,ck,dk);
%
% end slimju
