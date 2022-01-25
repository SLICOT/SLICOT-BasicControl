function [sysi,rcondD] = slinv(sys)
% SLINV  Computes the inverse of a standard system.
%
%        SYSI = SLINV(SYS)  computes the inverse SYSI = (Ai,Bi,Ci,Di) 
%        of a standard system SYS = (A,B,C,D), i.e.,
%
%        Ai = A - B*inv(D)*C,  Bi = -B*inv(D),
%        Ci = inv(D)*C,        Di = inv(D). 
%
%        [SYSI,RCONDD] = SLINV(SYS)  also computes RCONDD, the estimated
%        reciprocal condition number of the input/output matrix D of the
%        original system SYS.
%
%        See also SLDUAL, INVERT
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, July 2003, March, 2009.
%
%        Revisions:
%        -
%

ni   = nargin;
nout = nargout;
%
if min( ni, nout ) < 1,
    error(['Usage: SYSI          = SLINV(SYS)', sprintf('\n'),...
           '       [SYSI,RCONDD] = SLINV(SYS)'])
end
% 
[ A,B,C,D ] = ssdata( sys );
[ A,B,C,D,rcondD ] = invert( 2,A,B,C,D );
sysi = ss( A,B,C,D,sys );
%
% end slinv
