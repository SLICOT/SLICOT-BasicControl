function [sys] = slser(sys1,sys2,uplo)
% SLSER  Computes the series inter-connection of two systems 
%        in state-space form.
%
%        SYS = SLSER(SYS1,SYS2)  computes the series (cascaded)
%        inter-connection of two systems SYS1 = (A1,B1,C1,D1) 
%        and SYS2 = (A2,B2,C2,D2), obtaining the system 
%        SYS = (A,B,C,D), i.e.,
%
%        A = ( A1     0  ),   B = (  B1   ),
%            ( B2*C1  A2 )        ( B2*D1 )                     (1)
%
%        C = ( D2*C1  C2 ),   D = ( D2*D1 ).
%
%        SYS = SLSER(SYS1,SYS2,UPLO)    has the additional
%        input argument UPLO.
%
%        UPLO is an integer indicating whether the matrix A should
%        be obtained in the upper or lower block diagonal form:
%        UPLO = 1 :  Obtain A in the lower block diagonal form, see (1);
%        UPLO = 2 :  Obtain A in the upper block diagonal form, see (2).
%        Default:  UPLO = 1.
%
%        If UPLO = 2, then
%
%        A = ( A2  B2*C1 ),   B = ( B2*D1 ),
%            ( 0     A1  )        (  B1   )                     (2)
%
%        C = ( C2  D2*C1 ),   D = ( D2*D1 ).
%
%        Therefore, when A1 and A2 are block upper triangular
%        (for instance, in the real Schur form), the resulting
%        state matrix is also block upper triangular.
%
%        See also SLAPP, SLFEED, SLPAR, SLSPAR, SYSCONN
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
    error(['Usage: SYS = SLSER(SYS1,SYS2)', sprintf('\n'),...
           '       SYS = SLSER(SYS1,SYS2,UPLO)'])
end
% 
if ni == 2,
    uplo = 1;
end
%
[ A1,B1,C1,D1 ] = ssdata( sys1 );
[ A2,B2,C2,D2 ] = ssdata( sys2 );
[ A,B,C,D ] = sysconn( 1,A1,B1,C1,D1,A2,B2,C2,D2,uplo );
sys = ss( A,B,C,D,sys1 );
%
% end slser
