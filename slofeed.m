function [sysc,rcondi] = slofeed(sys,alpha,fbtype,F)
% SLOFEED Constructs the closed-loop system corresponding to 
%        the output feedback control law.
%
%        SYSC = SLOFEED(SYS,ALPHA,FBTYPE,F)  computes the
%        closed-loop system SYSC = (Ac,Bc,Cc,Dc) corresponding
%        to the system SYS = (A,B,C,D) and the output feedback
%        control law
%
%        u = ALPHA*F*y + v.
%
%        The output matrices are given by
%
%        Ac = A + ALPHA*B*F*E*C,  Bc = B + ALPHA*B*F*E*D, 
%        Cc = E*C,                Dc = E*D,
%     
%        where E = (I - ALPHA*D*F)**-1.
%
%        FBTYPE is an integer specifying the type of the 
%        feedback law:
%        FBTYPE = 1 :  Unitary output feedback (F = I); 
%        FBTYPE = 2 :  General output feedback.  
%        Default:  FBTYPE = 1.
%
%        The parameter F need not be specified if ALPHA = 0
%        and/or FBTYPE = 1.
%
%        [SYSC,RCONDI] = SLOFEED(SYS,ALPHA,FBTYPE,F)  also
%        returns the reciprocal condition number of the matrix 
%        I - ALPHA*D*F.
%
%        See also SLOSFEED, SYSFCONN
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
    error(['Usage: [SYSC,RCONDI] = SLOFEED(SYS,ALPHA,FBTYPE,F)', sprintf('\n'),...
           '        SYSC         = SLOFEED(SYS,ALPHA,FBTYPE)',   sprintf('\n'),...
           '        SYSC         = SLOFEED(SYS,ALPHA)'])
end
% 
if ni == 2,
    fbtype = 1;
end
% 
[ A,B,C,D ] = ssdata( sys );
[ p,m ] = size( D );
%
if fbtype == 1 && m ~= p,
    error('Matrix D must be square for unitary feedback')
end
%
if norm( D,1 ) == 0,
    jobd = 0;
else
    jobd = 1;
end
%
if jobd == 0,
    if fbtype == 1 || alpha == 0,
        [ Ac,Bc,Cc,rcondi ] = sysfconn( 2,jobd,fbtype,alpha,A,B,C );
    else
        [ Ac,Bc,Cc,rcondi ] = sysfconn( 2,jobd,fbtype,alpha,A,B,C,F );
    end
    Dc = zeros( p,m );
else
    if fbtype == 1 || alpha == 0,
        [ Ac,Bc,Cc,Dc,rcondi ] = sysfconn( 2,jobd,fbtype,alpha,A,B,C,D );
    else
        [ Ac,Bc,Cc,Dc,rcondi ] = sysfconn( 2,jobd,fbtype,alpha,A,B,C,D,F );
    end
end
sysc = ss( Ac,Bc,Cc,Dc,sys );
%
% end slofeed
