function [sys] = slspar(sys1,sys2,alpha)
% SLSPAR Computes rowwise concatenation (i.e., parallel
%        inter-connection on outputs, with separate inputs) 
%        of two systems in state-space form.
%
%        SYS = SLSPAR(SYS1,SYS2,ALPHA)  computes the state-space 
%        model SYS = (A,B,C,D) corresponding to the rowwise 
%        concatenation of two systems, SYS1 = (A1,B1,C1,D1), and 
%        SYS2 = (A2,B2,C2,D2), with the output equation for the 
%        second system multiplied by the scalar ALPHA, i.e.,
%
%            ( A1   0  )              ( B1  0  )
%        A = (         ) ,        B = (        ) ,
%            ( 0    A2 )              ( 0   B2 )
%      
%        C = ( C1   ALPHA*C2 ),   D = ( D1   ALPHA*D2 ).
%
%        Default:  ALPHA = 1.
%
%        See also SLAPP, SLFEED, SLPAR, SLSER, SYSCONN
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, July 2003.
%
%        Revisions: March 2009.
%        -
%

ni = nargin;
%
if ni < 2 || nargout < 1,
    error(['Usage: SYS = SLSPAR(SYS1,SYS2)', sprintf('\n'),...
           '       SYS = SLSPAR(SYS1,SYS2,ALPHA)'])
end
% 
if ni == 2,
    alpha = 1;
end
%
[ A1,B1,C1,D1 ] = ssdata( sys1 );
[ A2,B2,C2,D2 ] = ssdata( sys2 );
[ A,B,C,D ] = sysconn( 3,A1,B1,C1,D1,A2,B2,C2,D2,alpha );
sys = ss( A,B,C,D,sys1 );
%
% end slspar
