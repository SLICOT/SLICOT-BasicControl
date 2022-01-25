% Script for testing the performance of SLICOT discrete-time 
% Lyapunov solver.
% Random data sets of row/column dimensions multiplies of 10, or powers
% of 2, are used and the relative residuals and errors are computed
% and the execution times recorded. Symmetry is forced.
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

disp('Comparison between slstei and dlyap')
disp('-----------------------------------')
disp(' ')

disp('Tests with random matrices')
disp('--------------------------')
disp(' ')
disp('Please wait!')
disp(' ')

sp = '   ';  space = repmat(sp,Ntests10,1);  res = zeros(Ntests10,7);

l = 0;
for Nt = sz:sz+Ntests10-1,
  n = Nt*10; 
  l = l + 1;
  res(l,1) = n;
  
  A = rand(n,n);  X0 = rand(n,n);  X0 = X0 + X0';
  C = A'*X0*A - X0;  C = (C + C')/2.;
  normX0 = max(1,norm(X0,'fro'));
  
% SLICOT calculations.

  time = cputime;    
  x = slstei(A,C);
  res(l,2) = cputime - time;  
  res(l,4) = norm(A'*x*A - x - C,'fro')/normX0;
  res(l,6) = norm(x - X0,'fro')/normX0;
  
% MATLAB calculations.
  
  time = cputime;    
  X = dlyap(A',-C);
  res(l,3) = cputime - time;  
  res(l,5) = norm(A'*X*A - X - C,'fro')/normX0;
  res(l,7) = norm(X - X0,'fro')/normX0;

end  

disp('------------------------------------------------------------------------')
disp('            Time         Relative residual          Relative error      ')   
disp('------------------------------------------------------------------------')
disp('    n  slstei  dlyap    slstei       dlyap        slstei       dlyap    ')
disp('------------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,num2str(res(:,2),'%5.2f'), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4)), ...
  space,num2str(res(:,5)),space,num2str(res(:,6)), ...	
  space,num2str(res(:,7))] )	
disp('------------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page
page_plot(res,'slstei',1)
uiwait
disp(' ')
pause,

rp_slstei1 = res;
clear res space
space = repmat(sp,Ntests2,1);  res = zeros(Ntests2,7);

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  l = l + 1;
  res(l,1) = n;
  
  A = rand(n,n);  X0 = rand(n,n);  X0 = X0 + X0';
  C = A'*X0*A - X0;  C = (C + C')/2.;
  normX0 = max(1,norm(X0,'fro'));
  
% SLICOT calculations.

  time = cputime;    
  x = slstei(A,C);
  res(l,2) = cputime - time;  
  res(l,4) = norm(A'*x*A - x - C,'fro')/normX0;
  res(l,6) = norm(x - X0,'fro')/normX0;
  
% MATLAB calculations.
  
  time = cputime;    
  X = dlyap(A',-C);
  res(l,3) = cputime - time;  
  res(l,5) = norm(A'*X*A - X - C,'fro')/normX0;
  res(l,7) = norm(X - X0,'fro')/normX0;

end  

disp('-------------------------------------------------------------------------')
disp('            Time         Relative residual          Relative error      ')   
disp('-------------------------------------------------------------------------')
disp('    n   slstei  dlyap    slstei       dlyap        slstei       dlyap    ')
disp('-------------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,num2str(res(:,2),'%5.2f'), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4)), ...
  space,num2str(res(:,5)),space,num2str(res(:,6)), ...	
  space,num2str(res(:,7))] )	
disp('-------------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slstei',1)
uiwait
disp(' ')
pause,

rp_slstei2 = res;
clear res space
space = repmat(sp,Ntests2,1);  res = zeros(Ntests2,7);

% A is in a real Schur form

disp('Now, A is in a real Schur form')
disp(' ')
pause,

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  l = l + 1;
  res(l,1) = n;
  
  A = rand(n,n);  A = schur(A);  X0 = rand(n,n);  X0 = X0 + X0';
  C = A'*X0*A - X0;  C = (C + C')/2.;
  normX0 = max(1,norm(X0,'fro'));
  
% SLICOT calculations.

  time = cputime;    
  x = slstei(A,C,1);
  res(l,2) = cputime - time;  
  res(l,4) = norm(A'*x*A - x - C,'fro')/normX0;
  res(l,6) = norm(x - X0,'fro')/normX0;
  
% MATLAB calculations.
  
  time = cputime;    
  X = dlyap(A',-C);
  res(l,3) = cputime - time;  
  res(l,5) = norm(A'*X*A - X - C,'fro')/normX0;
  res(l,7) = norm(X - X0,'fro')/normX0;

end  

disp('-------------------------------------------------------------------------')
disp('            Time         Relative residual          Relative error      ')   
disp('-------------------------------------------------------------------------')
disp('    n   slstei  dlyap    slstei       dlyap        slstei       dlyap    ')
disp('-------------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,num2str(res(:,2),'%5.2f'), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4)), ...
  space,num2str(res(:,5)),space,num2str(res(:,6)), ...	
  space,num2str(res(:,7))] )	
disp('-------------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slstei',1)
uiwait
disp(' ')
pause,

rp_slstei3 = res;
clear res A C X x X0 normX0

perfb_slstei