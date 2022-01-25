% Script for testing the performance of SLICOT Sylvester solvers.
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
A = rand(n,n);  X0 = rand(n,n);  X0 = X0 + X0';
C = A'*X0 + X0*A;  C = (C + C')/2.;
sllyap(A,C);

disp('Comparison between slsylv and lyap')
disp('----------------------------------')
disp(' ')

disp('Tests with random matrices')
disp('--------------------------')
disp(' ')
disp('Please wait!')
disp(' ')

n  = (sz:sz+Ntests10-1)*10;  
nt = numel(sz:sz+Ntests10-1)*numel(n/Ntests10 : n/Ntests10 : 0.7*n);
sp = '   ';  space = repmat(sp,nt,1);  res = zeros(nt,8);

l = 0;
for Nt = sz:sz+Ntests10-1,
  n = Nt*10; 
  for mc = n/Ntests10 : n/Ntests10 : 0.7*n
    m = floor(mc);
    l = l + 1;
    res(l,1) = n;  res(l,2) = m;
    
    A = rand(n,n);  B = rand(m,m);  X0 = rand(n,m);
    C = A*X0 + X0*B;
    normX0 = max(1,norm(X0,'fro'));
    
%   SLICOT calculations.
    
    time = cputime;    
    x = slsylv(A,B,C);
    res(l,3) = cputime - time;  
    res(l,5) = norm(A*x + x*B - C,'fro')/normX0;
    res(l,7) = norm(x - X0,'fro')/normX0;
    
%   MATLAB calculations.
    
    time = cputime;    
    X = lyap(A,B,-C);
    res(l,4) = cputime - time;  
    res(l,6) = norm(A*X + X*B - C,'fro')/normX0;
    res(l,8) = norm(X - X0,'fro')/normX0;
    
  end
end  

disp('----------------------------------------------------------------------------')
disp('                 Time         Relative residual          Relative error    ')   
disp('----------------------------------------------------------------------------')
disp('    n    m  slsylv  lyap     slsylv        lyap        slsylv        lyap')
disp('----------------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,int2str(res(:,2)), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4),'%5.2f'), ...
  space,num2str(res(:,5)),space,num2str(res(:,6)), ...	
  space,num2str(res(:,7)),space,num2str(res(:,8))] )	
disp('----------------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page
page_plot(res,'slsylv',1)
uiwait
disp(' ')
pause,

rp_slsylv1 = res;
clear res space
n  = fix(2.^(sz2:sz2+Ntests2-1));  
nt = numel(sz2:sz2+Ntests2-1)*numel(n/Ntests2 : n/Ntests2 : 0.7*n);
space = repmat(sp,nt,1);  res = zeros(nt,8);

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  for mc = n/Ntests2 : n/Ntests2 : 0.7*n
    m = floor(mc);
    l = l + 1;
    res(l,1) = n;  res(l,2) = m;
    
    A = rand(n,n);  B = rand(m,m);  X0 = rand(n,m);
    C = A*X0 + X0*B;
    normX0 = max(1,norm(X0,'fro'));
    
%   SLICOT calculations.
    
    time = cputime;    
    x = slsylv(A,B,C);
    res(l,3) = cputime - time;  
    res(l,5) = norm(A*x + x*B - C,'fro')/normX0;
    res(l,7) = norm(x - X0,'fro')/normX0;
    
%   MATLAB calculations.
    
    time = cputime;    
    X = lyap(A,B,-C);
    res(l,4) = cputime - time;  
    res(l,6) = norm(A*X + X*B - C,'fro')/normX0;
    res(l,8) = norm(X - X0,'fro')/normX0;
  end
end  

disp('-----------------------------------------------------------------------------')
disp('                 Time         Relative residual          Relative error    ')   
disp('-----------------------------------------------------------------------------')
disp('    n    m   slsylv  lyap     slsylv        lyap        slsylv        lyap')
disp('-----------------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,int2str(res(:,2)), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4),'%5.2f'), ...
  space,num2str(res(:,5)),space,num2str(res(:,6)), ...	
  space,num2str(res(:,7)),space,num2str(res(:,8))] )	
disp('-----------------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slsylv',1)
uiwait
disp(' ')
pause,

rp_slsylv2 = res;
clear res space
n  = fix(2.^(sz2:sz2+Ntests2-1));  
nt = numel(sz2:sz2+Ntests2-1)*numel(n/Ntests2/2 : n/Ntests2/2);
space = repmat(sp,nt,1);  res = zeros(nt,8);

% A is in a real Schur form, and B in Hessenberg form

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  for mc = n/Ntests2/2 : n/Ntests2/2
    m = floor(mc);
    l = l + 1;
    res(l,1) = n;  res(l,2) = m;
    
    A = rand(n,n);  A = schur(A);  B = triu(rand(m,m),-1);  
    X0 = rand(n,m);
    C = A*X0 + X0*B;
    normX0 = max(1,norm(X0,'fro'));
    
%   SLICOT calculations.
    
    time = cputime;    
    x = slsylv(A,B,C,[1 2]);
    res(l,3) = cputime - time;  
    res(l,5) = norm(A*x + x*B - C,'fro')/normX0;
    res(l,7) = norm(x - X0,'fro')/normX0;
    
%   MATLAB calculations.
    
    time = cputime;    
    X = lyap(A,B,-C);
    res(l,4) = cputime - time;  
    res(l,6) = norm(A*X + X*B - C,'fro')/normX0;
    res(l,8) = norm(X - X0,'fro')/normX0;
  end
end  

if exist('res','var')
disp('Now, A is in a real Schur form, and B in Hessenberg form')
disp(' ')
pause,

disp('----------------------------------------------------------------------------')
disp('                 Time         Relative residual          Relative error    ')   
disp('----------------------------------------------------------------------------')
disp('    n    m   slsylv  lyap     slsylv        lyap        slsylv        lyap')
disp('----------------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,int2str(res(:,2)), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4),'%5.2f'), ...
  space,num2str(res(:,5)),space,num2str(res(:,6)), ...	
  space,num2str(res(:,7)),space,num2str(res(:,8))] )	
disp('----------------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slsylv',1)
uiwait
disp(' ')
pause,

rp_slsylv3 = res;
clear res space
end
n  = fix(2.^(sz2:sz2+Ntests2-1));  
nt = numel(sz2:sz2+Ntests2-1)*numel(n/Ntests2 : n/Ntests2 : 0.7*n);
space = repmat(sp,nt,1);  res = zeros(nt,8);

% Both A and B are in real Schur form

disp('Now, both A and B are in real Schur form')
disp(' ')
pause,

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  for mc = n/Ntests2 : n/Ntests2 : 0.7*n
    m = floor(mc);
    l = l + 1;
    res(l,1) = n;  res(l,2) = m;
    
    A = rand(n,n);  A = schur(A);  B = rand(m,m);  B = schur(B);  
    X0 = rand(n,m);
    C = A*X0 + X0*B;
    normX0 = max(1,norm(X0,'fro'));
    
%   SLICOT calculations.
    
    time = cputime;    
    x = slsylv(A,B,C,[1 1]);
    res(l,3) = cputime - time;  
    res(l,5) = norm(A*x + x*B - C,'fro')/normX0;
    res(l,7) = norm(x - X0,'fro')/normX0;
    
%   MATLAB calculations.
    
    time = cputime;    
    X = lyap(A,B,-C);
    res(l,4) = cputime - time;  
    res(l,6) = norm(A*X + X*B - C,'fro')/normX0;
    res(l,8) = norm(X - X0,'fro')/normX0;
  end
end  

disp('-----------------------------------------------------------------------------')
disp('                 Time         Relative residual          Relative error    ')   
disp('-----------------------------------------------------------------------------')
disp('    n    m   slsylv  lyap     slsylv        lyap        slsylv        lyap')
disp('-----------------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,int2str(res(:,2)), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4),'%5.2f'), ...
  space,num2str(res(:,5)),space,num2str(res(:,6)), ...	
  space,num2str(res(:,7)),space,num2str(res(:,8))] )	
disp('-----------------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slsylv',2)
uiwait
disp(' ')
pause,

rp_slsylv4 = res;
clear res A B C X X0 x

