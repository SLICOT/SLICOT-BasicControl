function [sysc] = sld2c(sys,alpha,beta)
% SLD2C  Performs a bilinear transformation of a discrete-time system 
%        to a continuous-time system.
%
%        SYSC = SLD2C(SYS,ALPHA,BETA)  computes the transformed 
%        system SYSC = (Ao,Bo,Co,Do) corresponding to the original
%        system SYS = (A,B,C,D),
%
%        Ao = BETA*inv(ALPHA*I + A) * (A - ALPHA*I)
%        Bo = sqrt(2*ALPHA*BETA) * inv(ALPHA*I + A) * B
%        Co = sqrt(2*ALPHA*BETA) * C * inv(ALPHA*I + A)
%        Do = D - C * inv(ALPHA*I + A) * B
%
%        which is equivalent to the bilinear transformation
%
%                      z - ALPHA
%        z -> s = BETA --------- ,
%                      z + ALPHA
%
%        of one transfer matrix onto the other.
%
%        The parameters ALPHA and BETA must be non-zero.
%        Recommended values for stable systems: ALPHA = 1, BETA = 1.
%        These values are used by default if ALPHA and BETA are missing.
%
%        See also SLC2D, CONDIS
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, July 2003.
%
%        Revisions: March 2009.
%        -
%

ni   = nargin;
nout = nargout;
%
if min( ni, nout ) < 1,
    error(['Usage: SYSC = SLD2C(SYS)', sprintf('\n'),...
           '       SYSC = SLD2C(SYS,ALPHA,BETA))'])
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
[ A,B,C,D ] = condis( 1,A,B,C,D,alpha,beta );
sysc = ss( A,B,C,D,0 );
%
% end sld2c
