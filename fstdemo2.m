echo on
%   This demonstration shows the use of the SLICOT-based Matlab function 
%   for computing the generator and the Cholesky factor for the inverse
%   of a positive definite block Toeplitz matrix.

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
%       These functions call the mexfile fstoep.  Use function fstgen here.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

more on
help fstgen
more off

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Call fstgen, given the first block-row of BT.
    
G = fstgen(Tr);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The generator of the inverse of BT, G, is')
G
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       The generator can be used to compute an approximate inverse of BT:

approxinv = btoeplitz(G(1:k,:)', G(1:k,1:k) * eye(k,size(G,2))) ...
          * btoeplitz(G(1:k,:)', G(1:k,1:k) * eye(k,size(G,2)))' ...
          - btoeplitz(G(k+1:2*k,:)', G(k+1:2*k,1:k) * eye(k,size(G,2))) ...
          * btoeplitz(G(k+1:2*k,:)', G(k+1:2*k,1:k) * eye(k,size(G,2)))';

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  
disp(['The error norm:   norm( approxinv * BT - eye(size(BT)) ) = ',...
      num2str(norm( approxinv * BT - eye(size(BT)) ))])
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Compute also the Cholesky factor, L, for the inverse of BT.
    
[G,L] = fstgen(Tr);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The lower Cholesky factor of inv(BT), L, is')
L
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('Check with Matlab function chol, Lm = inv(chol(BT))''')
Lm = inv(chol(BT))';
Lm

disp(' '),  
disp(['The relative error norm:   norm(L - Lm)/norm(Lm) = ',...
      num2str(norm(L - Lm)/norm(Lm))])
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Use fstgen to solve the systems  X*BT = B.
%       First, generate X, compute B, and then solve the equations.
%       The columns of X contain the first n*k*m natural numbers.

m = n;  X0 = 1 : k*n*m;  X0 = reshape(X0,m,k*n);  B = X0*BT;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end
    
[G,L,X] = fstgen(Tr,B);

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

%       Call fstgen, given the first block-column of BT.
    
G = fstgen(Tc);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       The generator can be used to compute an approximate inverse of BT:

approxinv = btoeplitz(G(:,1:k), G(1:k,1:k) * eye(k,size(G,1))) ...
          * btoeplitz(G(:,1:k), G(1:k,1:k) * eye(k,size(G,1)))' ...
          - btoeplitz(G(:,k+1:2*k), G(1:k,k+1:2*k) * eye(k,size(G,1))) ...
          * btoeplitz(G(:,k+1:2*k), G(1:k,k+1:2*k) * eye(k,size(G,1)))';

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  
disp(['The error norm:   norm( approxinv * BT - eye(size(BT)) ) = ',...
      num2str(norm( approxinv * BT - eye(size(BT)) ))])
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Compute also the upper Cholesky factor, L, for the inverse of BT.
    
[G,L] = fstgen(Tc);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The uper Cholesky factor of inv(BT), L, is')
L
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('Check with Matlab function chol, Lm = inv(chol(BT))')
Lm = inv(chol(BT));

disp(' '),  
disp(['The relative error norm:   norm(L - Lm)/norm(Lm) = ',...
      num2str(norm(L - Lm)/norm(Lm))])
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Use fstgen to solve the systems  BT*X = B.
%       First, generate X, compute B, and then solve the equations.
%       The columns of X contain the first n*k*m natural numbers.

X0 = 1 : k*n*m;  X0 = reshape(X0,k*n,m);  B = BT*X0;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end
    
[G,L,X] = fstgen(Tc,B);

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



