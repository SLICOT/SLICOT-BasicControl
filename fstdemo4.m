echo on
%   This demonstration shows the use of the SLICOT-based Matlab function 
%   for updating the factorizations after adding new blocks in a
%   positive definite block Toeplitz matrix.

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

k = 2;  n = 3;  Tr = [ [ 1000 1; 1 1000] [26 27; 28 29] [36 37; 38 39] ];
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
%       These functions call the mexfile fstoep.  Use function fstupd here.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

more on
help fstupd
more off

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Call fstupd, given the first block-row of BT.
    
[R,G,H,CSH] = fstupd(Tr);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The Cholesky factor of BT, R, is')
R
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Now, add one new block (m = 1), and get a 4-by-4 block matrix BTA.
%       The first block-row of BTA has n+m = 4 blocks:
%       the block i has elements containing the figure i, i = 1 : 4.

m = 1;  Tra = [46 47; 48 49];
Tca = Tra';  BTA = btoeplitz([Tc; Tca],[Tr Tra]);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The first block-row, [Tr Tra], is')
[Tr Tra]
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The first block-column, [Tc; Tca], is')
[Tc; Tca]
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The block Toeplitz matrix, BTA, is')
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
BTA
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Call fstupd again for updating, given the first block-row of BTA.
    
[Ru,Gu,Hu,CSHu] = fstupd(Tra,R,G,H,CSH);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       The final results will be checked later.
%       Now, add two new blocks (p = 2), and get a 6-by-6 block matrix BTE.
%       The first block-row of BTE has n+m+p = 6 blocks.

p = 2;  Tre = [ [36 37; 38 39] [26 27; 28 29] ];
Tce = Tre';  BTE = btoeplitz([Tc; Tca; Tce],[Tr Tra Tre]);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The block Toeplitz matrix, BTE, is')
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
BTE
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Call fstupd again for updating, given the first block-row of BTE.
    
[Re,Ge,He,CSHe] = fstupd(Tre,Ru,Gu,Hu,CSHu);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The final Cholesky factor (of BTE), Re, is')
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
Re
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('Check with Matlab function chol, Rm = chol(BTE)')
echo on
Rm = chol(BTE);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
Rm
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  
disp(['The relative error norm:   norm(Re - Rm)/norm(Rm) = ',...
      num2str(norm(Re - Rm)/norm(Rm))])
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Use fstupd to solve the systems  X*BTE = B.
%       First, generate X, compute B, and then solve the equations.
%       The columns of X contain the first k*(n+m+p)*nrhs natural numbers.

nrhs = 2;  X0 = 1 : k*(n+m+p)*nrhs;  X0 = reshape(X0,nrhs,k*(n+m+p));  B = X0*BTE;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end
    
[Re,X] = fstupd(Tre,Ru,Gu,Hu,CSHu,B);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The solution of X*BTE = B is')
X
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  
disp(['The relative error norm:   norm(X - X0)/norm(X0) = ',...
      num2str(norm(X - X0)/norm(X0))])
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Now, call fstupd, given the first block-column of BT.
    
[R,G,H,CSH] = fstupd(Tc);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The lower Cholesky factor of BT, R, is')
R
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Update for one new block (m = 1), i.e., for the 4-by-4 block matrix BTA,
%       given the first block-column of BTE.
    
[Ru,Gu,Hu,CSHu] = fstupd(Tca,R,G,H,CSH);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Update for two new blocks (p = 2), i.e., for the 6-by-6 block matrix BTE,
%       given the first block-column of BTE.
    
[Re,Ge,He,CSHe] = fstupd(Tce,Ru,Gu,Hu,CSHu);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The final lower Cholesky factor (of BTE), Re, is')
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
Re
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('Check with Matlab result, Rm')
disp(' '),  
disp(['The relative error norm:   norm(Re - Rm'')/norm(Rm) = ',...
      num2str(norm(Re - Rm')/norm(Rm))])
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Use fstupd to solve the systems  BTE*X = B.
%       First, generate X, compute B, and then solve the equations.
%       The columns of X contain the first k*(n+m+p)*nrhs natural numbers.

X0 = 1 : k*(n+m+p)*nrhs;  X0 = reshape(X0,k*(n+m+p),nrhs);  B = BTE*X0;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end
    
[Re,X] = fstupd(Tce,Ru,Gu,Hu,CSHu,B);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' '),  disp('The solution of BTE*X = B is')
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



