echo on
%   This demonstration shows the use of the SLICOT-based Matlab function 
%   for solving a linear system with a block Toeplitz matrix.

echo off
%   RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000.
%
%   Revisions: V. Sima, March 2009.
%   

echo on
global pause_wait  % This could be used in pause(n) command.
                   % If pause_wait < 0, standard command pause is used (default).
                   % Any key should then be pressed to continue.
                   
if ~exist('pause_wait', 'var') || isempty(pause_wait),  pause_wait = -1;  end

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Consider an n-by-n block Toeplitz matrix, BT, with k-by-k blocks.  
%       Let k = 2, n = 3.  The first block-row of BT has n = 3 blocks:
%       the block i has elements containing the figure i, i = 1, 2, 3.

k = 2;  n = 3;  Tr = [ [ 100 1; 1 100] [26 27; 28 29] [36 37; 38 39]];
Tc = Tr';  BT = btoeplitz(Tc,Tr);

echo off
disp(' '),  disp('The first block-row, Tr, is')
Tr
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The first block-column, Tc, is')
Tc
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The block Toeplitz matrix, BT, is')
BT
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Available SLICOT-based Matlab functions:  fstchol, fstgen, fstupd, fstsol.
%       These functions call the mexfile fstoep.  Use function fstsol here.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

more on
help fstsol
more off

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Use fstsol to solve the systems  X*BT = B, given the first 
%       block-row of BT.
%       First, generate X, compute B, and then solve the equations.
%       The columns of X contain the first n*k*m natural numbers.

m = n;  X0 = 1 : k*n*m;  X0 = reshape(X0,m,k*n);  B = X0*BT;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end
    
X = fstsol(Tr,B);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The solution of X*BT = B is')
X
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  
disp(['The relative error norm:   norm(X - X0)/norm(X0) = ',...
      num2str(norm(X - X0)/norm(X0))])
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Use fstsol to solve the systems  BT*X = B, given the first 
%       block-column of BT.
%       First, generate X, compute B, and then solve the equations.
%       The columns of X contain the first n*k*m natural numbers.

X0 = 1 : k*n*m;  X0 = reshape(X0,k*n,m);  B = BT*X0;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end
    
X = fstsol(Tc,B);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The solution of BT*X = B is')
X
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  
disp(['The relative error norm:   norm(X - X0)/norm(X0) = ',...
      num2str(norm(X - X0)/norm(X0))])
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off



