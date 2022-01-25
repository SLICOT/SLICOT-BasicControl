function [sys] = slfeed(sys1,sys2,alpha)
% SLFEED Computes the feedback inter-connection of two 
%        systems in state-space form.
%
%        SYS = SLFEED(SYS1,SYS2,ALPHA)  computes the feedback 
%        inter-connection of the systems SYS1 = (A1,B1,C1,D1) 
%        and SYS2 = (A2,B2,C2,D2), obtaining the system 
%        SYS = (A,B,C,D), defined by
%
%        A = ( A1  -  ALPHA*B1*E12*D2*C1       -  ALPHA*B1*E12*C2    ),
%            (        B2*E21*C1            A2  -  ALPHA*B2*E21*D1*C2 )
%
%        B = (  B1*E12    ),
%            (  B2*E21*D1 )
%
%        C = (  E21*C1     -  ALPHA*E21*D1*C2 ),
%
%        D = (  E21*D1 ),
%
%        E21  =  inv( I + ALPHA*D1*D2 ) and
%        E12  =  inv( I + ALPHA*D2*D1 ) = I - ALPHA*D2*E21*D1.
%
%        ALPHA = +1 corresponds to positive feedback, and
%        ALPHA = -1 corresponds to negative feedback.
%        Default:  ALPHA = -1.
%
%        Taking the orders of A1 and/or A2 equal to 0 will solve the
%        constant plant and/or constant feedback cases.
%
%        See also SLAPP, SLPAR, SLSER, SLSPAR, SYSCONN
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, July 2003.
%
%        Revisions: 
%        03-03-2009.
%

ni = nargin;
%
if ni < 2 || nargout < 1,
    error(['Usage: SYS = SLFEED(SYS1,SYS2)', sprintf('\n'),...
           '       SYS = SLFEED(SYS1,SYS2,ALPHA)'])
end
% 
if ni == 2,
    alpha = -1;
end
%
[ A1,B1,C1,D1 ] = ssdata( sys1 );
[ A2,B2,C2,D2 ] = ssdata( sys2 );
[ A,B,C,D ] = sysconn( 2,A1,B1,C1,D1,A2,B2,C2,D2,alpha );
sys = ss( A,B,C,D,sys1 );
%
% end slfeed
