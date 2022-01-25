function [sysd] = sldual(sys)
% SLDUAL Computes the dual of a standard system.
%
%        SYSD = SLDUAL(SYS)  computes the dual system SYSD = (Ao,Bo,Co,Do) 
%        of a standard system SYS = (A,B,C,D), i.e.,
%
%        Ao = A';  Bo = C';  Co = B';  Do = D'.
%
%        See also SLINV, INVERT
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
    error('Usage: SYSD = SLDUAL(SYS)')
end
% 
[ A,B,C,D ] = ssdata( sys );
[ A,B,C,D ] = invert( 1,A,B,C,D );
sysd = ss( A,B,C,D,sys );
%
% end sldual
