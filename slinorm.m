function [lnorm,fpeak] = slinorm(sys,tol,fpeak0)
%SLINORM  L-infinity system norm.
%
%   SLINORM(SYS) is the L-infinity norm of SYS, i.e., the peak gain
%   of its frequency response (as measured by the largest singular 
%   value in the MIMO case).
%
%   SLINORM(SYS,TOL) specifies a relative accuracy TOL for the 
%   computed infinity norm (TOL = 1.e-4 by default).
%       
%   [NINF,FPEAK] = SLINORM(SYS) also returns the frequency FPEAK
%   where the gain achieves its peak value NINF.
%
%   [NINF,FPEAK] = SLINORM(SYS,TOL,FPEAK0) specifies a relative 
%   accuracy TOL and an initial estimate FPEAK0 of the frequency FPEAK 
%   where the gain achieves its peak value (FPEAK0 = 0 by default).
%       

%  Reference:
%      Bruisma, N.A., and M. Steinbuch, ``A Fast Algorithm to Compute
%      the Hinfinity-Norm of a Transfer Function Matrix,'' Syst. Contr. 
%      Letters, 14 (1990), pp. 287-293.

%   RELEASE 2.0 of SLICOT Robust Control Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Author: A. Varga.
%
%   Revisions:
%   V. Sima, June 2001.
%   A. Varga, Aug. 2001.
%   V. Sima, May 2003, Jan. 2007.
%

if nargin < 3,
   fpeak0 = [0;1];	
elseif isinf(fpeak0),
   fpeak0(2) = 0;  fpeak0(1) = 1;	
end
if nargin < 2,
   tol = 1.e-4;
else
   tol = max(100*eps,tol);
end

Ts = sys.Ts;
if Ts ~= 0
   dico = 2;
else
   dico = 1;
end
[a,b,c,d,e] = dssdata(sys);

[gpeak,fpeak1] = linorm(a,e,b,c,d,dico,0,2,fpeak0,tol);
if ~( gpeak(2) == 0 ),
   lnorm = gpeak(1);
else
   lnorm = Inf;
end
if ~( fpeak1(2) == 0 ),
   fpeak = fpeak1(1);
else
   fpeak = Inf;
end
if Ts ~= 0
   fpeak = fpeak/abs(Ts);
end


