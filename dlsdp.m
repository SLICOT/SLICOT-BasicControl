% DLSDP.F - Gateway function for Loop Shaping Design of
%           discrete-time systems using SLICOT routine SB10KD.
%
% Matlab call:  
%   [AK,BK,CK,DK,(RCOND)] = dlsdp(A,B,C,factor)
%
% Purpose:
%     To compute the matrices of the positive feedback controller
%
%              | Ak | Bk |
%          K = |----|----|
%              | Ck | Dk |
%
%     for the discrete-time shaped plant
%
%              | A | B |
%          G = |---|---|
%              | C | 0 |
%
%     in the McFarlane/Glover Loop Shaping Design Procedure.
%
% Input parameters: 
%   A      - the n-by-n system state matrix A.
%   B      - the n-by-m system input matrix B.
%   C      - the p-by-n system output matrix C.
%   factor - = 1 implies that an optimal controller is required
%            > 1 implies that a suboptimal controller is required
%                achieving a performance FACTOR less than optimal.
%
% Output parameters:
%   AK    - the n-by-n controller state matrix Ak.
%   BK    - the n-by-p controller input matrix Bk.
%   CK    - the m-by-n controller output matrix Ck.
%   DK    - the m-by-p controller matrix Dk.
%   RCOND - (optional) a vector containing estimates of the reciprocal
%           condition numbers of the matrices which have to be inverted
%           during the computation of the controller.
%           RCOND(1) contains an estimate of the reciprocal condition
%                    number of the linear system of equations from
%                    which the solution of the P-Riccati equation is
%                    obtained
%           RCOND(2) contains an estimate of the reciprocal condition
%                    number of the linear system of equations from
%                    which the solution of the Q-Riccati equation is
%                    obtained
%           RCOND(3) contains an estimate of the reciprocal condition
%                    number of the linear system of equations from
%                    which the solution of the X-Riccati equation is
%                    obtained
%           RCOND(4) contains an estimate of the reciprocal condition
%                    number of the matrix Rx + Bx'*X*Bx
%
% References
%   [1] D.W. Gu, P.Hr. Petkov, D.W. Gu and M.M. Konstantinov.
%       Hinf Loop Shaping Design procedure routines in SLICOT.
%       NICONET Report 1999-15, November 1999.

% RELEASE 2.0 of SLICOT Robust Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
