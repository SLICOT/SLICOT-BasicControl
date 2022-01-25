function [sysb,U,V] = slgsbal(sys,job,thresh)
%SLGSBAL Balance the system matrix for a descriptor system
%        SYS = (A,E,B,C).
%
%        [SYSB,U,V] = SLGSBAL(SYS,JOB,THRESH)  computes a diagonal 
%        equivalence transformation, 
%           (U*A*V - lambda U*E*V, U*B, C*V),
%        applied to the system (A-lambda E,B,C) to make the rows
%        and columns of system pencil
%           S =  ( A  B ) - lambda ( E  0 ) ,
%                ( C  0 )          ( 0  0 )
%        as close in norm as possible.
%
%        SYSB = SLGSBAL(SYS,JOB,THRESH)  performs the same computations,
%        but the diagonal matrices U and V are not returned.
%
%        JOB is an integer parameter specifying the part of the system
%        matrix used to compute the scaling factors for balancing:
%           JOB = 1  : use the matrices A, E and B.
%           JOB = 2  : use the matrices A, E and C.
%           JOB = 3  : use the matrices A and E.
%           otherwise, use all matrices A, E, B, and C.
%        Default:  JOB = 0.
%
%        THRESH is a real parameter specifying a threshold value:
%        elements with magnitude less than or equal to THRESH
%        are ignored for balancing.
%        Default: THRESH = 0.
%
%        Note: This function cannot work with a singular matrix E. 
%              Use GSYSTRA instead.
%        Note that GSYSTRA can work with rectangular matrices A and E.
%
%        See also SLGSHES, SLGSQRQ, SLGSRSF, SLGSSVD, GSYSTRA
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 30-04-2003.
%

ni = nargin;
%
if ni <= 1 
   job = 0;
end
if ni <= 2
   thresh = 0;
end
flag(1) = job;
flag(2) = thresh;
%
[A,B,C,D,E] = dssdata( sys );
%
if nargout <= 1
   [A,E,B,C] = gsystra( 0,A,E,B,C,flag );
else  
   [A,E,B,C,u,v] = gsystra( 0,A,E,B,C,flag );
   U = diag( u );
   if nargout == 3,  V = diag( v );  end
end 
%
sysb = dss( A,B,C,D,E,sys );
%
% end slgsbal
