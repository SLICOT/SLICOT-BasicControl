function [sys] = slapp(sys1,sys2)
% SLAPP  Appends two systems in state-space form.
%
%        SYS = SLAPP(SYS1,SYS2)  appends two systems 
%        SYS1 = (A1,B1,C1,D1) and SYS2 = (A2,B2,C2,D2), 
%        obtaining the system SYS = (A,B,C,D), i.e.,
%
%            ( A1   0  )         ( B1  0  )
%        A = (         ) ,   B = (        ) ,
%            ( 0    A2 )         ( 0   B2 )
%      
%            ( C1   0  )         ( D1  0  )
%        C = (         ) ,   D = (        ) .
%            ( 0    C2 )         ( 0   D2 )
%
%        See also SLFEED, SLPAR, SLSER, SLSPAR, SYSCONN
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, July 2003.
%
%        Revisions:
%        March 2009.
%

if nargin < 2 || nargout < 1,
    error('Usage: SYS = SLAPP(SYS1,SYS2)')
end
% 
[ A1,B1,C1,D1 ] = ssdata( sys1 );
[ A2,B2,C2,D2 ] = ssdata( sys2 );
[ A,B,C,D ] = sysconn( 5,A1,B1,C1,D1,A2,B2,C2,D2 );
sys = ss( A,B,C,D,sys1 );
%
% end slapp
