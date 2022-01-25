% DLSDZ.F - Gateway function for Loop Shaping Design of
%           discrete-time systems using SLICOT routine SB10ZD.
%
% Matlab call:
%   [AK,BK,CK,DK(,RCOND)] = dlsdz(A,B,C,D,factor(,tol))
%
% Purpose:
%   To compute the matrices of the positive feedback controller
%
%              | Ak | Bk |
%          K = |----|----|
%              | Ck | Dk |
%
%   for the discrete-time shaped plant
%
%              | A | B |
%          G = |---|---|
%              | C | D |
%
%   in the McFarlane/Glover Loop Shaping Design Procedure.
%
% Input parameters:
%   A      - the n-by-n system state matrix A.
%   B      - the n-by-m system input matrix B.
%   C      - the p-by-n system output matrix C.
%   D      - the p-by-m system matrix D.
%   factor - = 1 implies that an optimal controller is required;
%                (not reccomended);
%            > 1 implies that a suboptimal controller is required
%                achieving a performance FACTOR less than optimal.
%   tol    - (optional) tolerance used for checking the nonsingularity
%            of the matrices to be inverted.
%            Default:  tol = sqrt(epsilon_machine), where
%            epsilon_machine is the relative machine precision.
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
%                    obtained;
%           RCOND(2) contains an estimate of the reciprocal condition
%                    number of the linear system of equations from
%                    which the solution of the Q-Riccati equation is
%                    obtained;
%           RCOND(3) contains an estimate of the reciprocal condition
%                    number of the matrix (gamma^2-1)*In - P*Q;
%           RCOND(4) contains an estimate of the reciprocal condition
%                    number of the matrix Rx + Bx'*X*Bx;
%           RCOND(5) contains an estimate of the reciprocal condition
%                                                ^
%                    number of the matrix Ip + D*Dk;
%           RCOND(6) contains an estimate of the reciprocal condition
%                                              ^
%                    number of the matrix Im + Dk*D.
%
% References
%   [1] Gu, D.-W., Petkov, P.H., and Konstantinov, M.M.
%       On discrete H-infinity loop shaping design procedure routines.
%       Technical Report 00-6, Dept. of Engineering, Univ. of
%       Leicester, UK, 2000.
%

% RELEASE 2.0 of SLICOT Robust Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
