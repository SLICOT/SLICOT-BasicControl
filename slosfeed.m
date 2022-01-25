function [sysc,rcondi] = slosfeed(sys,G,H,alpha,fbtype,F,beta,K)
% SLOSFEED Constructs the closed-loop system corresponding to 
%        the mixed output and state feedback control law.
%
%        SYSC = SLOSFEED(SYS,G,H,ALPHA,FBTYPE,F,BETA,K)  computes
%        the closed-loop system SYSC = (Ac,Bc,Cc,Dc) corresponding
%        to the system SYS = (A,B,C,D) and the mixed output and
%        state feedback control law
%
%        u = ALPHA*F*y + BETA*K*x + G*v,
%        z = H*y.
%
%        The output matrices are given by
%
%        Ac = AC + BETA*BC*K,      Bc = BC*G,  
%        Cc = H*(CC + BETA*DC*K),  Dc = H*DC*G,
%     
%        where
%        E  = (I - ALPHA*D*F)**-1,
%        AC = A + ALPHA*B*F*E*C,  BC = B + ALPHA*B*F*E*D, 
%        CC = E*C,                DC = E*D.
%
%        FBTYPE is an integer specifying the type of the 
%        feedback law:
%        FBTYPE = 1 :  Unitary output feedback (F = I); 
%        FBTYPE = 2 :  General output feedback.  
%        Default:  FBTYPE = 1.
%
%        The parameter K need not be specified if BETA = 0
%        (the default value for BETA).
%
%        [SYSC,RCONDI] = SLOSFEED(SYS,G,H,ALPHA,FBTYPE,F,BETA,K) 
%        also returns the reciprocal condition number of the matrix 
%        I - ALPHA*D*F.
%
%        See also SLOFEED, SYSFCONN
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
if ni < 4 || nargout < 1,
    error(['Usage: [SYSC,RCONDI] = SLOSFEED(SYS,G,H,ALPHA,FBTYPE,F,BETA,K)', sprintf('\n'),...
           '        SYSC         = SLOSFEED(SYS,G,H,ALPHA,FBTYPE,F)',   sprintf('\n'),...
           '        SYSC         = SLOSFEED(SYS,G,H,ALPHA)'])
end
% 
if ni == 4,
    fbtype = 1;
    beta   = 0;
elseif ni <= 6,
    beta   = 0;
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
        if beta == 0,
            [ Ac,Bc,Cc,rcondi ] = sysfconn( 1,jobd,fbtype,alpha,beta, ...
                                            A,B,C,G,H );
        else
            [ Ac,Bc,Cc,rcondi ] = sysfconn( 1,jobd,fbtype,alpha,beta, ...
                                            A,B,C,G,H,K );
        end
    else
        if beta == 0,
            [ Ac,Bc,Cc,rcondi ] = sysfconn( 1,jobd,fbtype,alpha,beta, ...
                                            A,B,C,G,H,F );
        else
            [ Ac,Bc,Cc,rcondi ] = sysfconn( 1,jobd,fbtype,alpha,beta, ...
                                            A,B,C,G,H,F,K );
        end
    end
    Dc = zeros( size( H,1 ),size( G,2 ) );
else
    if fbtype == 1 || alpha == 0,
        if beta == 0,
            [ Ac,Bc,Cc,Dc,rcondi ] = sysfconn( 1,jobd,fbtype,alpha,beta, ...
                                               A,B,C,D,G,H );
        else
            [ Ac,Bc,Cc,Dc,rcondi ] = sysfconn( 1,jobd,fbtype,alpha,beta, ...
                                               A,B,C,D,G,H,K );
        end
    else
        if beta == 0,
            [ Ac,Bc,Cc,Dc,rcondi ] = sysfconn( 1,jobd,fbtype,alpha,beta, ...
                                               A,B,C,D,G,H,F );
        else
            [ Ac,Bc,Cc,Dc,rcondi ] = sysfconn( 1,jobd,fbtype,alpha,beta, ...
                                               A,B,C,D,G,H,F,K );
        end
    end
end
sysc = ss( Ac,Bc,Cc,Dc,sys );
%
% end slosfeed
