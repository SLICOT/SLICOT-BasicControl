% LDSIMT.F - MEX-function for computing the output response of a
%            linear discrete-time system using SLICOT routines TF01MD
%            and TF01ND.
%
%   [Y(,x)] = ldsimt(A,B,C,D,U(,x,HessA))
%
%   LDSIMT computes the output vector sequence y(1), y(2),..., y(t)
%
%        x(k+1) = A x(k) + B u(k)
%        y(k)   = C x(k) + D u(k),
%
%   given an initial state vector x(1), and the input vector sequence
%   u(1), u(2),..., u(t), where y(k) and u(k) are vectors of length p
%   and m, respectively. The input trajectories are given as
%
%        U = ( u(1)  u(2) ...  u(t) )
%
%   and the output trajectories result in similarly. Optionally, the
%   matrix A can be given in upper or lower Hessenberg form.
%
%   Description of input parameters:
%   A      - the n-by-n state matrix A. If HessA = 1 or 2 only the upper
%            or lower Hessenberg part, respectively, must be defined.
%   B      - the n-by-m input matrix B.
%   C      - the p-by-n output matrix C.
%   D      - the p-by-m input-output matrix D.
%   U      - the m-by-t matrix U.
%   x      - (optional) the initial state x(1).
%            Default: x = 0.
%   HessA  - (optional) scalar indicating whether the matrix A is
%            general or in an upper/lower Hessenberg form:
%            = 0 :  general matrix;
%            = 1 :  upper Hessenberg matrix;
%            = 2 :  lower Hessenberg matrix.
%            Default: HessA = 0.
%
%   Description of output parameters:
%   Y      - the p-by-t output matrix Y.
%   x      - (optional) the final state x(t+1).
%
%   Comments
%   This gateway function stores the input and output trajectories in
%   a transposed form comparing to the LDSIM function. Moreover, A could
%   be given in a Hessenberg form.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003.
%
% Revisions:
%   -
%
