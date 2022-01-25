function [sys] = slpar(sys1,sys2,alpha)
% SLPAR  Computes the parallel inter-connection of two systems 
%        in state-space form.
%
%        SYS = SLPAR(SYS1,SYS2,ALPHA)  computes the state-space 
%        model SYS = (A,B,C,D) corresponding to the sum 
%        G = G1 + ALPHA*G2, where G, G1, and G2 are the 
%        transfer-function matrices of the corresponding 
%        state-space models SYS, SYS1 = (A1,B1,C1,D1), and 
%        SYS2 = (A2,B2,C2,D2), respectively, i.e.,
%
%            ( A1   0  )             ( B1 )
%        A = (         ) ,       B = (    ) ,
%            ( 0    A2 )             ( B2 )
%
%        C = ( C1  ALPHA*C2 ) ,  D = D1 + ALPHA*D2 .
%
%        Default:  ALPHA = 1.
%
%        See also SLAPP, SLFEED, SLSER, SLSPAR, SYSCONN
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, July 2003.
%
%        Revisions:
%        March 2009.
%

ni = nargin;
%
if ni < 2 || nargout < 1,
    error(['Usage: SYS = SLPAR(SYS1,SYS2)', sprintf('\n'),...
           '       SYS = SLPAR(SYS1,SYS2,ALPHA)'])
end
% 
if ni == 2,
    alpha = 1;
end
%
[ A1,B1,C1,D1 ] = ssdata( sys1 );
[ A2,B2,C2,D2 ] = ssdata( sys2 );
[ A,B,C,D ] = sysconn( 4,A1,B1,C1,D1,A2,B2,C2,D2,alpha );
sys = ss( A,B,C,D,sys1 );
%
% end slpar
