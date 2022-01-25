% Script for testing the performance of SLICOT generalized 
% discrete-time Lyapunov solver.
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

disp('Performance of slgest (generalized discrete-time Lyapunov solver)')
disp('-----------------------------------------------------------------')
disp(' ')

disp('Tests with random matrices')
disp('--------------------------')
disp(' ')
disp('Please wait!')
disp(' ')

sp = '    ';  space = repmat(sp,Ntests10,1);  res = zeros(Ntests10,5);

l = 0;
for Nt = sz:sz+Ntests10-1,
  n = Nt*10; 
  l = l + 1;
  res(l,1) = n;
  
  A = rand(n,n);  E = rand(n,n);  Xo = rand(n,n);  Xo = Xo + Xo';
  C = A'*Xo*A - E'*Xo*E;  C = (C + C')/2.;
  normXo = max(1,norm(Xo,'fro'));
  
% SLICOT calculations.

  time = cputime;    
  [X,sep] = slgest(A,E,C);
  res(l,2) = cputime - time;  
  res(l,3) = norm(A'*X*A - E'*X*E - C,'fro')/normXo;
  res(l,4) = norm(X - Xo,'fro')/normXo;
  res(l,5) = sep;
end  

disp('--------------------------------------------------------')
disp('    n     Time   Rel. residual  Rel. error      sep')
disp('--------------------------------------------------------')
	                    
disp( [space,int2str(res(:,1)), ...
  space,num2str(res(:,2),'%5.2f'), ...
  space,num2str(res(:,3)),space,num2str(res(:,4)) ...
  space,num2str(res(:,5))] )	
disp('--------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page
page_plot(res,'slgest',1)
uiwait
disp(' ')
pause,

rp_slgest1 = res;
clear res space
space = repmat(sp,Ntests2,1);  res = zeros(Ntests2,4);

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  l = l + 1;
  res(l,1) = n;
  
  A = rand(n,n);  E = rand(n,n);  Xo = rand(n,n);  Xo = Xo + Xo';
  C = A'*Xo*A - E'*Xo*E;  C = (C + C')/2.;
  normXo = max(1,norm(Xo,'fro'));
  
% SLICOT calculations.

  time = cputime;    
  X = slgest(A,E,C);
  res(l,2) = cputime - time;  
  res(l,3) = norm(A'*X*A - E'*X*E - C,'fro')/normXo;
  res(l,4) = norm(X - Xo,'fro')/normXo;
end  

disp('-------------------------------------------')
disp('    n     Time   Rel. residual  Rel. error')
disp('-------------------------------------------')
	                    
disp( [space,int2str(res(:,1)), ...
  space,num2str(res(:,2),'%5.2f'), ...
  space,num2str(res(:,3)),space,num2str(res(:,4))] )	
disp('-------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slgest',1)
uiwait
disp(' ')
pause,

rp_slgest2 = res;
clear res space
space = repmat(sp,Ntests2,1);  res = zeros(Ntests2,4);

% (A,E) generalized Schur

disp('Now, (A,E) is in a generalized real Schur form')
disp(' ')
pause,

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  l = l + 1;
  res(l,1) = n;
  
  A = rand(n,n);  A = schur(A);  E = triu(rand(n,n)) + n*eye(n);  
  Xo = rand(n,n);  Xo = Xo + Xo';
  C = A'*Xo*A - E'*Xo*E;  C = (C + C')/2.;
  normXo = max(1,norm(Xo,'fro'));
  
% SLICOT calculations.

  time = cputime;    
  X = slgest(A,E,C,1);
  res(l,2) = cputime - time;  
  res(l,3) = norm(A'*X*A - E'*X*E - C,'fro')/normXo;
  res(l,4) = norm(X - Xo,'fro')/normXo;
end  

disp('-------------------------------------------')
disp('    n     Time   Rel. residual  Rel. error')
disp('-------------------------------------------')
	                    
disp( [space,int2str(res(:,1)), ...
  space,num2str(res(:,2),'%5.2f'), ...
  space,num2str(res(:,3)),space,num2str(res(:,4))] )	
disp('-------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slgest',1)
uiwait
disp(' ')
pause,

rp_slgest3 = res;
clear res A C E X Xo normXo

perfb_slgest