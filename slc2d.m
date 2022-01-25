function [sysd] = slc2d(sys,alpha,beta)
% SLC2D  Performs a bilinear transformation of a continuous-time system 
%        to a discrete-time system.
%
%        SYSD = SLC2D(SYS,ALPHA,BETA)  computes the transformed 
%        system SYSD = (Ao,Bo,Co,Do) corresponding to the original
%        system SYS = (A,B,C,D),
%
%        Ao = ALPHA*inv(BETA*I - A) * (BETA*I + A)
%        Bo = sqrt(2*ALPHA*BETA) * inv(BETA*I - A) * B
%        Co = sqrt(2*ALPHA*BETA) * C * inv(BETA*I - A)
%        Do = D + C * inv(BETA*I - A) * B
% 
%        which is equivalent to the bilinear transformation
% 
%                       BETA + s
%        s -> z = ALPHA -------- 
%                       BETA - s
% 
%        of one transfer matrix onto the other.
%
%        The parameters ALPHA and BETA must be non-zero.
%        Recommended values for stable systems: ALPHA = 1, BETA = 1.
%        These values are used by default if ALPHA and BETA are missing.
%
%        See also SLD2C, CONDIS
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, July 2003.
%
%        Revisions:
%        March 2009.
%

ni   = nargin;
nout = nargout;
%
if min( ni, nout ) < 1,
    error(['Usage: SYSD = SLC2D(SYS)', sprintf('\n'),...
           '       SYSD = SLC2D(SYS,ALPHA,BETA))'])
end
% 
if ni == 1,
    alpha = 1;
    beta  = 1;
elseif ni == 2,
    beta  = 1;
end
% 
[ A,B,C,D ] = ssdata( sys );
[ A,B,C,D ] = condis( 2,A,B,C,D,alpha,beta );
sysd = ss( A,B,C,D,-1 );
%
% end slc2d
