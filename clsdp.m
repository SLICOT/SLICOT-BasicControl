% CLSDP.F - Gateway function for Loop Shaping Design of
%           continuous-time systems using SLICOT routine SB10ID.
%
% Matlab call:  
%   [AK,BK,CK,DK,(RCOND)] = clsdp(A,B,C,D,factor)
%
% Purpose:
%     To compute the matrices of the positive feedback controller
%
%              | Ak | Bk |
%          K = |----|----|
%              | Ck | Dk |
%
%     for the shaped plant
%
%              | A | B |
%          G = |---|---|
%              | C | D |
%
%     in the McFarlane/Glover Loop Shaping Design Procedure.
%
% Input parameters: 
%   A      - the n-by-n system state matrix A.
%   B      - the n-by-m system input matrix B.
%   C      - the p-by-n system output matrix C.
%   D      - the p-by-m system matrix D.
%   factor - = 1 implies that an optimal controller is required
%            > 1 implies that a suboptimal controller is required
%                achieving a performance FACTOR less than optimal.
%
% Output parameters:
%   AK    - the nk-by-nk controller state matrix Ak.
%   BK    - the nk-by-np controller input matrix Bk.
%   CK    - the m-by-nk controller output matrix Ck.
%   DK    - the m-by-np controller matrix Dk.
%   RCOND - (optional) a vector containing estimates of the reciprocal
%           condition numbers of the Riccati equations which have to be
%           solved during the computation of the controller.
%           RCOND(1) contains an estimate of the reciprocal condition
%                    number of the X-Riccati equation,
%           RCOND(2) contains an estimate of the reciprocal condition
%                    number of the Z-Riccati equation.
%
% References
%   [1] D.W. Gu, P.Hr. Petkov, D.W. Gu and M.M. Konstantinov. 
%       Hinf Loop Shaping Design procedure routines in SLICOT.
%       NICONET Report 1999-15, November 1999.

% RELEASE 2.0 of SLICOT Robust Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%


