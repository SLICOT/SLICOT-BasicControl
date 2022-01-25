function F = deadbt( sys, tol, bal )
%DEADBT  Constructs the minimum norm feedback matrix 
%        performing "deadbeat control" on a system 
%        SYS = (A,B,C,D).
%
%        F = DEADBT(SYS)  constructs the minimum norm deadbeat 
%        feedback matrix F for the system SYS. 
%        
%        F = DEADBT(SYS,TOL,BAL)  has additional input parameters. 
%        
%        TOL is a lower bound on the reciprocal condition numbers,
%        used to decide matrix ranks in the reductions.
%        Default:  TOL = N*N*machine_epsilon, where N = size(A,1).
%
%        BAL is an integer scalar indicating whether the
%        (A,B)-pair should be balanced, as follows:
%        BAL = 0 :  use balancing;
%        BAL = 1 :  do not use balancing.
%        Default:  BAL = 0.
%
%        Note: If the system is already in the controllability
%              (upper) staircase form, use DEADBEAT instead,
%              which provides more options.
% 
%        See also DEADBEAT
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Sep-26-2003.
%        Revised 
%        V. Sima Mar-2-2009.
%

%
ni = nargin;
if ni < 1 || nargout ~= 1,  
    error( ['Usage: F = DEADBT(SYS)', sprintf('\n'), ...
            '       F = DEADBT(SYS,TOL,BAL)' ] );
end
%
if ni == 1,
    tol = 0;
    bal = 0;
elseif ni == 2,
    bal = 0;
end
%
[ A, B, C, D ] = ssdata( sys );  
%
IStair = 0;  WithU = 1;  WithV = 1;
F = deadbeat( A, B, IStair, tol, WithU, WithV, bal );
%
% end deadbt
