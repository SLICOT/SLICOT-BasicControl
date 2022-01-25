% LINORM.F - Gateway function for computation of the L_infinity norm of
%            a continuous-time or discrete-time system, either standard
%            or in the descriptor form, using SLICOT routine AB13DD.
%
%   [gpeak(,fpeak)] = linorm(A,E,B,C,D(,dico,systype,equil,fpeak0,tol))
%
% Purpose:
%   To compute the L-infinity norm of the continuous-time or 
%   discrete-time system, either standard or in the descriptor form,
%
%                                     -1
%        G(lambda) = C*( lambda*E - A ) *B + D .
%
% Input parameters: 
%   A      - the n-by-n system state matrix A.
%   E      - the n-by-n descriptor matrix E of the system, or an empty
%            matrix (i.e., with 0 rows and/or columns), in which case
%            E is taken as an identity matrix of order n (standard
%            system). The contents of E are ignored if it is non-empty
%            and systype = 1, but its order is then checked out.
%   B      - the n-by-m system input matrix B.
%   C      - the p-by-n system output matrix C.
%   D      - the p-by-m system matrix D.
%   dico   - (optional) specifies the type of the system:
%            = 1, continuous-time system;
%            = 2, discrete-time system.
%            Default: dico = 1.
%   systype- (optional) specifies whether or not the system is of
%            descriptor type: 
%            = 0 : descriptor system;
%            = 1 : standard system (E = I).
%            Default: systype = 0.
%   equil  - (optional) specifies whether the user wishes to 
%            preliminarily equilibrate the system (A,E,B,C) or (A,B,C):
%             = 1:  do not perform equilibration;
%             = 2:  perform equilibration (scaling).
%            Default: equil = 1.
%   fpeak0 - (optional) array of length 2 containing an estimate of the
%            frequency where the gain of the frequency response would
%            achieve its peak value. If fpeak0(2) = 0, the frequency is
%            infinite.
%            Default: fpeak0 = [0; 1].
%   tol    - (optional) tolerance used to set the accuracy in
%            determining the norm.
%            Default: sqrt(epsilon_machine) where epsilon_machine is
%            the relative machine precision.
%
% Output parameters:
%   gpeak  - array of length 2 containing the value of the L_infinity
%            norm of the system. If gpeak(2) = 0, the norm is infinite.
%   fpeak  - (optional) array of length 2 containing the frequency for
%            which the frequency response achieves its peak value gpeak.
%            If fpeak(2) = 0, the frequency is infinite.

% RELEASE 2.0 of SLICOT Robust Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%

