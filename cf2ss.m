function [ sys, rcnd ] = cf2ss( sysQ, sysR, right )
%CF2SS   Constructs the state-space representation of a system from
%        the factors of its left or right coprime factorization.
%
%        SYS = CF2SS(SYSQ,SYSR)  computes the state-space representation
%        of a system, SYS = (A,B,C,D), from the factors Q and R of its 
%        left coprime factorization, where Q and R are given by their
%        state-space representations, SYSQ = (AQR,BQ,CQR,DQ) and
%        SYSR = (AQR,BR,CQR,DR), respectively. Specifically, if G, Q, and
%        R are the transfer-function matrices corresponding to the systems
%        SYS, SYSQ, and SYSR, then G = R^(-1)*Q.
%        
%        [SYS,RCND] = CF2SS(SYSQ,SYSR,RIGHT)  has additional input and 
%        output parameters.
%        
%        RIGHT specifies the computations to be performed.
%        RIGHT = 0 :  compute the state-space representation from the
%                     factors Q and R of its left coprime factorization
%                     (see above);
%        RIGHT = 1 :  compute the state-space representation from the
%                     factors Q and R of its right coprime factorization,
%                     where Q and R are given by their state-space
%                     representations, SYSQ = (AQR,BQR,CQ,DQ) and
%                     SYSR = (AQR,BQR,CR,DR), and SYS is state-space
%                     representation of the transfer-function matrix
%                     G = Q*R^(-1).
%        Default: RIGHT = 0.
%        
%        RCND is an estimate of the reciprocal condition number
%        of the matrix DR.
%
%        See also CFSYS
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Dec-22-2003.
%        Revised Mar-02-2009.
%

ni = nargin;
if ni < 2 || nargout < 1,  
    error( ['Usage: SYS = CF2SS(SYSQ,SYSR)', sprintf('\n'), ...
            '       [SYS,RCND] = CF2SS(SYSQ,SYSR,RIGHT)' ] );
end
%
if ni == 2,  right = 0;  end
task = right + 1;
%
[ A,  B,  C,  D  ] = ssdata( sysQ );
[ AR, BR, CR, DR ] = ssdata( sysR );
%
if norm( sysQ.Ts - sysR.Ts ) ~= 0,  
    error( 'SYSQ and SYSR must be both continuous- or discrete-time systems' ),  
end
if norm( A - AR, 1 ) ~= 0,  
    error( 'SYSQ and SYSR must have the same state matrix' ),  
end
if task == 1 && norm( C - CR, 1 ) ~= 0,  
    error( 'SYSQ and SYSR must have the same state/output matrix' ),  
end
if task == 2 && norm( B - BR, 1 ) ~= 0,  
    error( 'SYSQ and SYSR must have the same input/state matrix' ),  
end
%
if task == 1,  
    [ A, B, C, D, rcnd ] = cfsys( task, A, B, C, D, BR, DR );
else  
    [ A, B, C, D, rcnd ] = cfsys( task, A, B, C, D, CR, DR );
end
sys = ss( A, B, C, D, sysQ.Ts ); 
%
% end cf2ss
