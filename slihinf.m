function [sysk,gmin,rcnd] = slihinf(sys,ncon,nmeas,gmax,job,discr,gtol,actol)
%SLIHINF Compute the H_inf optimal controller given a state
%        space model. 
%
%        SYSK = SLIHINF(SYS,NCON,NMEAS,GMAX)  computes the H_inf optimal 
%        controller SYSK, containing the matrices Ak,Bk,Ck,Dk in a
%        packed form, for a continuous-time system SYS (given in packed
%        form), with NCON control inputs and NMEAS measurements, for
%        the initial value GMAX of GAMMA on input. It is assumed that
%        GMAX is sufficiently large so that the controller is
%        admissible. GMAX >= 0.
%
%        [SYSK,GMIN,RCND] = SLIHINF(SYS,NCON,NMEAS,GMAX,JOB,DISCR,GTOL,ACTOL)
%        has additional input and output arguments:
%
%        JOB indicates the strategy for reducing the GAMMA value:
%        JOB = 1 : Use bisection method for decreasing GAMMA from
%                  GMAX to GMIN until the closed-loop system leaves
%                  stability.
%        JOB = 2 : Scan from GMAX to 0 trying to find the minimal GAMMA
%                  for which the closed-loop system retains stability.
%        JOB = 3 : First bisection, then scanning.
%        JOB = 0 : Suboptimal H_inf controller only.
%        Default:  JOB = 3.
%
%        DISCR indicates the type of system:
%        DISCR = 0 : continuous-time (default);
%        DISCR = 1 : discrete-time.
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
%        GMIN is the minimal GAMMA found.
%
%        RCND contains the reciprocal condition numbers. Specifically,
%        RCND(1) contains the reciprocal condition number of the
%                control transformation matrix;
%        RCND(2) contains the reciprocal condition number of the
%                measurement transformation matrix;
%        RCND(3) contains an estimate of the reciprocal condition
%                number of the X-Riccati equation;
%        RCND(4) contains an estimate of the reciprocal condition
%                number of the Y-Riccati equation.
%
%   Comments:
%   1) The discrete-time systems are handled via bilinear transformation
%   to continuous-time and then the controller obtained is discretized.
%   2) It is assumed that 
%      NP - NMEAS >= NCON  >= 0 and 
%      M  - NCON  >= NMEAS >= 0
%   where M and NP are the number of inputs and outputs, respectively,
%   of the system SYS.
%   3) The commands pck and unpck from the Mu-Analysis and Synthesis 
%   Toolbox are used.
% 
%   See also MUHOPT, SLIMJU

% RELEASE 2.0 of SLICOT Robust Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   A. Markovski, Technical University of Sofia, October 2003.
%
% Revisions:
%   V. Sima, Research Institute for Informatics, Bucharest, March 2004;
%   May 2005, March 2009.
%

ni   = nargin;
nout = nargout;
%
if ( ni < 4 || ni > 8 || nout < 1 || nout > 3 )
   error( [ 'Usage: [SYSK]           = SLIHINF(SYS,NCON,NMEAS,GMAX)',  sprintf('\n'),...
            '       [SYSK,GMIN]      = SLIHINF(SYS,NCON,NMEAS,GMAX,JOB,DISCR)',  sprintf('\n'),...
            '       [SYSK,GMIN,RCND] = SLIHINF(SYS,NCON,NMEAS,GMAX,JOB,DISCR,GTOL,ACTOL)' ] )
end
%
if ni == 4,
   job   = 3;
   discr = 0;
elseif ni == 5,
   discr = 0;
end
%
[a,b,c,d] = unpck( sys );
job = -job;
%
if ni <= 6
   [ak,bk,ck,dk,gmin,rcnd] = muHopt(job,discr,a,b,c,d,ncon,nmeas,gmax);
elseif ni == 7
   [ak,bk,ck,dk,gmin,rcnd] = muHopt(job,discr,a,b,c,d,ncon,nmeas,gmax,gtol);
elseif ni == 8
   [ak,bk,ck,dk,gmin,rcnd] = muHopt(job,discr,a,b,c,d,ncon,nmeas,gmax,gtol,actol);
end
%
sysk = pck( ak,bk,ck,dk );
%
% end slihinf
