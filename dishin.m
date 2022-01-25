% DISHIN.F - Gateway function for H_infinity or H_2 design of
%            discrete-time systems based on SLICOT routines
%
%   [AK,BK,CK,DK,(RCOND)] = dishin(task,A,B,C,D,ncon,nmeas,(gamma))
%
%   task = 1 :  [AK,BK,CK,DK,(RCOND)] = dishin(1,A,B,C,D,ncon,nmeas,
%                                              gamma)
%
%     To compute the matrices of an H-infinity (sub)optimal n-state
%     controller
%
%              | AK | BK |
%          K = |----|----|,
%              | CK | DK |
%
%     for the discrete-time system
%
%              | A  | B1  B2  |   | A | B |
%          P = |----|---------| = |---|---|,
%              | C1 | D11 D12 |   | C | D | 
%              | C2 | D21 D22 |
%
%     and for a given value of gamma, where B2 has column size of the
%     number of control inputs (ncon) and C2 has row size of the number
%     of measurements (nmeas) being provided to the controller.
%
%   task = 2 :  [AK,BK,CK,DK,(RCOND)] = dishin(2,A,B,C,D,ncon,nmeas)
%
%     To compute the matrices of the H2 optimal n-state controller
%
%              | AK | BK |
%          K = |----|----|,
%              | CK | DK |
%
%     for the discrete-time system
%
%              | A  | B1  B2  |   | A | B |
%          P = |----|---------| = |---|---|,
%              | C1 |  0  D12 |   | C | D | 
%              | C2 | D21 D22 |
%
%     where B2 has column size of the number of control inputs (ncon)
%     and C2 has row size of the number of measurements (nmeas) being
%     provided to the controller. 
%
% Input parameters: 
%   task  - integer option to determine the type of the design:
%           = 1 : H_infinity design;
%           = 2 : H_2 design.
%   A     - the n-by-n system state matrix A.
%   B     - the n-by-m system input matrix B.
%   C     - the p-by-n system output matrix C.
%   D     - the p-by-m system matrix D.
%   ncon  - the number of control inputs. m >= ncon >= 0,
%           p-nmeas >= ncon.
%   nmeas - the number of measurements. p >= nmeas >= 0,
%           m-ncon >= nmeas.
%   gamma - (task 1 only) the parameter gamma used in H_infinity design.
%           It is assumed that gamma is sufficiently large so that the
%           controller is admissible. gamma >= 0.
%
% Output parameters:
%   AK    - the n-by-n controller state matrix AK.
%   BK    - the n-by-nmeas controller input matrix BK.
%   CK    - the ncon-by-n controller output matrix CK.
%   DK    - the ncon-by-nmeas controller matrix DK.
%   RCOND - (optional) a vector containing estimates of the reciprocal
%           condition numbers of the matrices which are to be inverted
%           and estimates of the reciprocal condition numbers of the
%           Riccati equations which have to be solved during the
%           computation of the controller. (See the description of
%           the algorithm in [1].)
%
%           If task = 1 then
%             RCOND(1) contains the reciprocal condition number of the 
%                      matrix R3,
%             RCOND(2) contains the reciprocal condition number of the 
%                      matrix R1 - R2'*inv(R3)*R2,
%             RCOND(3) contains the reciprocal condition number of the 
%                      matrix V21,
%             RCOND(4) contains the reciprocal condition number of the 
%                      matrix St3,
%             RCOND(5) contains the reciprocal condition number of the 
%                      matrix V12, 
%             RCOND(6) contains the reciprocal condition number of the 
%                      matrix Im2 + DKHAT*D22,
%             RCOND(7) contains the reciprocal condition number of the 
%                      X-Riccati equation,
%             RCOND(8) contains the reciprocal condition number of the 
%                      Z-Riccati equation.
%
%           If task = 2 then
%             RCOND(1) contains the reciprocal condition number of the 
%                      control transformation matrix TU,
%             RCOND(2) contains the reciprocal condition number of the 
%                      measurement transformation matrix TY,
%             RCOND(3) contains the reciprocal condition number of the 
%                      matrix Im2 + B2'*X2*B2,
%             RCOND(4) contains the reciprocal condition number of the 
%                      matrix Ip2 + C2*Y2*C2',
%             RCOND(5) contains the reciprocal condition number of the 
%                      X-Riccati equation,
%             RCOND(6) contains the reciprocal condition number of the 
%                      Y-Riccati equation,
%             RCOND(7) contains the reciprocal condition number of the
%                      matrix Im2 + DKHAT*D22 .
%
% References
%   [1] P.Hr. Petkov, D.W. Gu and M.M. Konstantinov. Fortran 77 routines
%       for Hinf and H2 design of linear discrete-time control systems.
%       Report99-8, Department of Engineering, Leicester University,
%       April 1999.

% RELEASE 2.0 of SLICOT Robust Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
