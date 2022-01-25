% Script for testing the performance of SLICOT stable generalized
% discrete-time Lyapunov solver.
% Random data sets of row/column dimensions multiplies of 10, or powers
% of 2, are used and the relative residuals and errors are computed
% and the execution times recorded. Stability and symmetry are forced.
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

disp('Performance of slgsst (stable generalized discrete-time Lyapunov solver)')
disp('------------------------------------------------------------------------')
disp(' ')

disp('Tests with random matrices')
disp('--------------------------')
disp(' ')
disp('Please wait!')
disp(' ')

n  = (sz:sz+Ntests10-1)*10;  
nt = numel(sz:sz+Ntests10-1)*numel(n/Ntests10 : n/Ntests10 : n);
sp = '    ';  space = repmat(sp,nt,1);  res = zeros(nt,4);

l = 0;
for Nt = sz:sz+Ntests10-1,
  n = Nt*10; 
  for mc = n/Ntests10 : n/Ntests10 : n
    m = floor(mc);
    l = l + 1;
    res(l,1) = n;  res(l,2) = m;
    
    A = rand(n,n)/n;  E = rand(n,n) + n*eye(n);  C  = rand(m,n);
  
%   SLICOT calculations.

    time = cputime;    
    XC = slgsst(A,E,C);
    res(l,3) = cputime - time;  
    X = XC'*XC;
    res(l,4) = norm(A'*X*A - E'*X*E + C'*C,'fro')/max(1,norm(X,'fro'));
  end  
end  

disp('------------------------------------')
disp('    n     m     Time   Rel. residual')
disp('------------------------------------')
	                    
disp( [space,int2str(res(:,1)),space,int2str(res(:,2)), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4))] )	
disp('------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page
page_plot(res,'slgsst',1)
uiwait
disp(' ')
pause,

rp_slgsst1 = res;
clear res space
n  = fix(2.^(sz2:sz2+Ntests2-1));  
nt = numel(sz2:sz2+Ntests2-1)*numel(n/Ntests2 : n/Ntests2 : n);
space = repmat(sp,nt,1);  res = zeros(nt,4);

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  for mc = n/Ntests2 : n/Ntests2 : n
    m = floor(mc);
    l = l + 1;
    res(l,1) = n;  res(l,2) = m;
  
    A = rand(n,n)/n;  E = rand(n,n) + n*eye(n);  C  = rand(m,n);
  
%   SLICOT calculations.

    time = cputime;    
    XC = slgsst(A,E,C);
    res(l,3) = cputime - time;  
    X = XC'*XC;
    res(l,4) = norm(A'*X*A - E'*X*E + C'*C,'fro')/max(1,norm(X,'fro'));
  end  
end  

disp('------------------------------------')
disp('    n     m     Time   Rel. residual')
disp('------------------------------------')
	                    
disp( [space,int2str(res(:,1)),space,int2str(res(:,2)), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4))] )	
disp('------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slgsst',1)
uiwait
disp(' ')
pause,

rp_slgsst2 = res;
clear res space
space = repmat(sp,nt,1);  res = zeros(nt,4);

% (A,E) generalized Schur

disp('Now, (A,E) is in a generalized real Schur form')
disp(' ')
pause,

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  for mc = n/Ntests2 : n/Ntests2 : n
    m = floor(mc);
    l = l + 1;
    res(l,1) = n;  res(l,2) = m;
  
    A = triu(rand(n,n))/n;  E = triu(rand(n,n)) + n*eye(n);   
    C  = rand(m,n);
  
%   SLICOT calculations.

    time = cputime;    
    XC = slgsst(A,E,C,1);
    res(l,3) = cputime - time;  
    X = XC'*XC;
    res(l,4) = norm(A'*X*A - E'*X*E + C'*C,'fro')/max(1,norm(X,'fro'));
  end  
end  

disp('------------------------------------')
disp('    n     m     Time   Rel. residual')
disp('------------------------------------')
	                    
disp( [space,int2str(res(:,1)),space,int2str(res(:,2)), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4))] )	
disp('------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slgsst',1)
uiwait
disp(' ')
pause,

rp_slgsst3 = res;
clear res A C E X XC

perfb_slgsst