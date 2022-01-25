% Script for testing the performance of SLICOT generalized Sylvester 
% solver.
% Random data sets of row/column dimensions multiplies of 10, or powers
% of 2, are used and the relative residuals and errors are computed
% and the execution times recorded.
% The number of tests performed for each group of tests is defined 
% by the variables Ntests10 (for dimensions multiplies of 10) and
% Ntests2 (for dimensions powers of 2). Default values 5 and 4, 
% respectively.

%   RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   V. Sima, Katholieke Univ. Leuven, Belgium, June 1999.
%
%   Revisions:
%   V. Sima, May 2005, March 2009.
%   

def_Ntests

% First execution is needed to load the mexfile, without counting time.

n = 1;
A = rand(n,n);  E = rand(n,n);  X0 = rand(n,n);  X0 = X0 + X0';
C = A'*X0*E + E'*X0*A;  C = (C + C')/2.;
slgely(A,E,C);

disp('Performance of slgesg (generalized Sylvester solver)')
disp('----------------------------------------------------')
disp(' ')

disp('Tests with random matrices')
disp('--------------------------')
disp(' ')
disp('Please wait!')
disp(' ')

n  = (sz:sz+Ntests10-1)*10;  
nt = numel(sz:sz+Ntests10-1)*numel(n/Ntests10 : n/Ntests10 : n);
sp = '   ';  space = repmat(sp,nt,1);  res = zeros(nt,8);

l = 0;
for Nt = sz:sz+Ntests10-1,
  n = Nt*10; 
  for mc = n/Ntests10 : n/Ntests10 : n
    m = floor(mc);
    l = l + 1;
    res(l,1) = n;  res(l,2) = m;
    
    A = rand(n,n);    D = rand(n,n);   B = rand(m,m);   E = rand(m,m);
    Xo = rand(n,m);   Yo = rand(n,m);
    C = A*Xo - Yo*B;  F = D*Xo - Yo*E;
    normXo = max(1,norm(Xo,'fro'));
    normYo = max(1,norm(Yo,'fro'));
    
%   SLICOT calculations.
    
    time = cputime;    
    [X,Y,Dif] = slgesg(A,D,B,E,C,F);
    res(l,3) = cputime - time;  
    res(l,4) = norm(A*X - Y*B - C,'fro')/max(normXo,normYo);
    res(l,5) = norm(D*X - Y*E - F,'fro')/max(normXo,normYo);
    res(l,6) = norm(X - Xo,'fro')/normXo;
    res(l,7) = norm(Y - Yo,'fro')/normYo;
    res(l,8) = Dif;
  end
end  

disp('----------------------------------------------------------------------------------')
disp('             Time      Relative residuals         Relative errors           Dif')   
disp('----------------------------------------------------------------------------------')
disp('    n    m  slgesg     in X         in Y         in X         in Y')
disp('----------------------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,int2str(res(:,2)), ...
  space,num2str(res(:,3),'%5.2f'), ...
  space,num2str(res(:,4)),space,num2str(res(:,5)), ...	
  space,num2str(res(:,6)),space,num2str(res(:,7)), ...	
  space,num2str(res(:,8))] )	
disp('----------------------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page
page_plot(res,'slgesg',1)
uiwait
disp(' ')
pause,

rp_slgesg1 = res;
clear res space
n  = fix(2.^(sz2:sz2+Ntests2-1));  
nt = numel(sz2:sz2+Ntests2-1)*numel(n/Ntests2 : n/Ntests2 : n);
space = repmat(sp,nt,1);  res = zeros(nt,7);

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  for mc = n/Ntests2 : n/Ntests2 : n
    m = floor(mc);
    l = l + 1;
    res(l,1) = n;  res(l,2) = m;
    
    A = rand(n,n);    D = rand(n,n);   B = rand(m,m);   E = rand(m,m);
    Xo = rand(n,m);   Yo = rand(n,m);
    C = A*Xo - Yo*B;  F = D*Xo - Yo*E;
    normXo = max(1,norm(Xo,'fro'));
    normYo = max(1,norm(Yo,'fro'));
    
%   SLICOT calculations.
    
    time = cputime;    
    [X,Y] = slgesg(A,D,B,E,C,F);
    res(l,3) = cputime - time;  
    res(l,4) = norm(A*X - Y*B - C,'fro')/max(normXo,normYo);
    res(l,5) = norm(D*X - Y*E - F,'fro')/max(normXo,normYo);
    res(l,6) = norm(X - Xo,'fro')/normXo;
    res(l,7) = norm(Y - Yo,'fro')/normYo;
  end
end  

disp('----------------------------------------------------------------------')
disp('             Time      Relative residuals         Relative errors    ')   
disp('----------------------------------------------------------------------')
disp('    n    m  slgesg     in X         in Y         in X         in Y')
disp('----------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,int2str(res(:,2)), ...
  space,num2str(res(:,3),'%5.2f'), ...
  space,num2str(res(:,4)),space,num2str(res(:,5)), ...	
  space,num2str(res(:,6)),space,num2str(res(:,7))] )	
disp('----------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slgesg',1)
uiwait
disp(' ')
pause,

rp_slgesg2 = res;
clear res space
n  = fix(2.^(sz2:sz2+Ntests2-1));  
nt = numel(sz2:sz2+Ntests2-1)*numel(n/Ntests2 : n/Ntests2 : n);
space = repmat(sp,nt,1);  res = zeros(nt,7);

% (A,D) generalized Schur, (B,E) general

disp('Now, (A,D) is in a generalized real Schur form')
disp(' ')
pause,

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  for mc = n/Ntests2 : n/Ntests2 : n
    m = floor(mc);
    l = l + 1;
    res(l,1) = n;  res(l,2) = m;
    
    A = rand(n,n);    A = schur(A);   D = triu(rand(n,n)) + n*eye(n);   
    B = rand(m,m);    E = rand(m,m);
    Xo = rand(n,m);   Yo = rand(n,m);
    C = A*Xo - Yo*B;  F = D*Xo - Yo*E;
    normXo = max(1,norm(Xo,'fro'));
    normYo = max(1,norm(Yo,'fro'));
    
%   SLICOT calculations.
    
    time = cputime;    
    [X,Y] = slgesg(A,D,B,E,C,F,1);
    res(l,3) = cputime - time;  
    res(l,4) = norm(A*X - Y*B - C,'fro')/max(normXo,normYo);
    res(l,5) = norm(D*X - Y*E - F,'fro')/max(normXo,normYo);
    res(l,6) = norm(X - Xo,'fro')/normXo;
    res(l,7) = norm(Y - Yo,'fro')/normYo;
  end
end  

disp('----------------------------------------------------------------------')
disp('             Time      Relative residuals         Relative errors    ')   
disp('----------------------------------------------------------------------')
disp('    n    m  slgesg     in X         in Y         in X         in Y')
disp('----------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,int2str(res(:,2)), ...
  space,num2str(res(:,3),'%5.2f'), ...
  space,num2str(res(:,4)),space,num2str(res(:,5)), ...	
  space,num2str(res(:,6)),space,num2str(res(:,7))] )	
disp('----------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slgesg',1)
uiwait
disp(' ')
pause,

rp_slgesg3 = res;
clear res space
n  = fix(2.^(sz2:sz2+Ntests2-1));  
nt = numel(sz2:sz2+Ntests2-1)*numel(n/Ntests2 : n/Ntests2 : n);
space = repmat(sp,nt,1);  res = zeros(nt,7);

% (A,D), (B,E) generalized Schur

disp('Now, both (A,D) and (B,E) are in generalized real Schur forms')
disp(' ')
pause,

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  for mc = n/Ntests2 : n/Ntests2 : n
    m = floor(mc);
    l = l + 1;
    res(l,1) = n;  res(l,2) = m;
    
    A = rand(n,n);    A = schur(A);   D = triu(rand(n,n)) + n*eye(n); 
    B = rand(m,m);    B = schur(B);   E = triu(rand(m,m)) + m*eye(m); 
    Xo = rand(n,m);   Yo = rand(n,m);
    C = A*Xo - Yo*B;  F = D*Xo - Yo*E;
    normXo = max(1,norm(Xo,'fro'));
    normYo = max(1,norm(Yo,'fro'));
    
%   SLICOT calculations.
    
    time = cputime;    
    [X,Y] = slgesg(A,D,B,E,C,F,[1,1]);
    res(l,3) = cputime - time;  
    res(l,4) = norm(A*X - Y*B - C,'fro')/max(normXo,normYo);
    res(l,5) = norm(D*X - Y*E - F,'fro')/max(normXo,normYo);
    res(l,6) = norm(X - Xo,'fro')/normXo;
    res(l,7) = norm(Y - Yo,'fro')/normYo;
  end
end  

disp('----------------------------------------------------------------------')
disp('             Time      Relative residuals         Relative errors    ')   
disp('----------------------------------------------------------------------')
disp('    n    m  slgesg     in X         in Y         in X         in Y')
disp('----------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,int2str(res(:,2)), ...
  space,num2str(res(:,3),'%5.2f'), ...
  space,num2str(res(:,4)),space,num2str(res(:,5)), ...	
  space,num2str(res(:,6)),space,num2str(res(:,7))] )	
disp('----------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slgesg',2)
uiwait
disp(' ')
pause,

rp_slgesg4 = res;
clear res A B C D E F X Xo Y Yo normXo normYo

