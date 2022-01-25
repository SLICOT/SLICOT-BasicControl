% Script for testing the performance of SLICOT generalized 
% discrete-time Lyapunov solver using dtlex Benchmack collection.
% Data sets of row/column dimensions multiplies of 10, or powers
% of 2, are used and the relative residuals and errors are
% computed and the execution times recorded.
% The number of tests performed for each group of tests is defined 
% by the variables Ntests10 (for dimensions multiplies of 10) and
% Ntests2 (for dimensions powers of 2). Default values 4 and 2, 
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

Ntests10 = 4;
Ntests2 = 2;

def_Ntests

disp('Performance of slgest (generalized discrete-time Lyapunov solver)')
disp('-----------------------------------------------------------------')
disp(' ')

disp('Tests with dtlex Benchmack collection matrices')
disp('----------------------------------------------')
disp(' ')
disp('Please wait!')
disp(' ')

sp = '    ';  space = repmat(sp,Ntests10,1);  res = zeros(Ntests10,4);

nr = [4;3];
disp('Example 4.3 (default values)')
disp(' ')

for l = 1:Ntests10,
  parin = [l*10,1.5,1.5]; 
  [E,A,C,B,Xo] = dtlex(nr,parin);
  n = length(A);
  res(l,1) = n;
  
  normXo = max(1,norm(Xo,'fro'));
  
% SLICOT calculations.

  time = cputime;    
  X = slgest(A,E,C);
  res(l,2) = cputime - time;  
  res(l,3) = norm(A'*X*A - E'*X*E - C,'fro')/normXo;
  res(l,4) = norm(X - Xo,'fro')/normXo;
end  

disp('-------------------------------------------')
disp('    n     Time   Rel. residual  Rel. error ')
disp('--------------------------------------------')
	                    
disp( [space,int2str(res(:,1)), ...
  space,num2str(res(:,2),'%5.2f'), ...
  space,num2str(res(:,3)),space,num2str(res(:,4))] )	
disp('--------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slgest',1)
uiwait
disp(' ')
pause,

rpb_slgest1 = res;
clear res space
space = repmat(sp,Ntests2,1);  res = zeros(Ntests2,4);

for l = 1:Ntests2,
  j = 2^l;
  parin = [j*8,1.5,1.5]; 
  [E,A,C,B,Xo] = dtlex(nr,parin);
  n = length(A);
  res(l,1) = n;
  
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
page_plot(res,'slgest',2)
uiwait
disp(' ')
pause,

rpb_slgest2 = res;
clear res A B C E X Xo normXo
