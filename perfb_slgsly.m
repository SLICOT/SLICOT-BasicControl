% Script for testing the performance of SLICOT stable generalized 
% continuous-time Lyapunov solver using ctlex Benchmack collection.
% Data sets of row/column dimensions multiplies of 10, or powers
% of 2, are used and the relative residuals are computed and the
% execution times recorded.
% The number of tests performed for each group of tests is defined 
% by the variables Ntests10 (for dimensions multiplies of 10) and
% Ntests2 (for dimensions powers of 2). Default values 3 and 2, 
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

Ntests10 = 3;
Ntests2 = 2;

def_Ntests

disp('Performance of slgsly (stable generalized continuous-time Lyapunov solver)')
disp('--------------------------------------------------------------------------')
disp(' ')

disp('Tests with ctlex Benchmack collection matrices')
disp('----------------------------------------------')
disp(' ')
disp('Please wait!')
disp(' ')

sp = '    ';  space = repmat(sp,Ntests10,1);  res = zeros(Ntests10,3);

nr = [4;4];
disp('Example 4.4 (default values)')
disp(' ')

for l = 1:Ntests10,
  parin = [l*4,1.5]; 
  [E,A,C,B] = ctlex(nr,parin);
  n = length(A);
  res(l,1) = n;
  
% SLICOT calculations.

  time = cputime;    
  X = slgsly(A,E,B);
  res(l,2) = cputime - time;  
  res(l,3) = norm(A'*X'*X*E + E'*X'*X*A + B'*B,'fro') ...
             /max(1,norm(X'*X,'fro'));
end  

disp('-------------------------------')
disp('    n     Time   Rel. residual ')
disp('-------------------------------')
	                    
disp( [space,int2str(res(:,1)), ...
  space,num2str(res(:,2),'%5.2f'), ...
  space,num2str(res(:,3))] )	
disp('-------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slgsly',1)
uiwait
disp(' ')
pause,

rpb_slgsly1 = res;
clear res space
space = repmat(sp,Ntests2,1);  res = zeros(Ntests2,3);

for l = 1:Ntests2,
  j = 2^l;
  parin = [j*3,1.5]; 
  [E,A,C,B] = ctlex(nr,parin);
  n = length(A);
  res(l,1) = n;
  
% SLICOT calculations.

  time = cputime;    
  X = slgsly(A,E,B);
  res(l,2) = cputime - time;  
  res(l,3) = norm(A'*X'*X*E + E'*X'*X*A + B'*B,'fro') ...
             /max(1,norm(X'*X,'fro'));
end  

disp('-------------------------------')
disp('    n     Time   Rel. residual ')
disp('-------------------------------')
	                    
disp( [space,int2str(res(:,1)), ...
  space,num2str(res(:,2),'%5.2f'), ...
  space,num2str(res(:,3))] )	
disp('-------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slgsly',2)
uiwait
disp(' ')
pause,

rpb_slgsly2 = res;
clear res A B C E X
