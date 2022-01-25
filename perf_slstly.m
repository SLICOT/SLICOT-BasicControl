% Script for testing the performance of SLICOT stable continuous-time 
% Lyapunov solver.
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

ver7 = version;  ver7 = ver7(1) >= '7';

if ver7,
  disp('Comparison between slstly and lyapchol')
  disp('--------------------------------------')
  disp(' ')
  disp('Note:  The errors are evaluated for X = U''*U')
else
  disp('Comparison between slstly and lyap')
  disp('----------------------------------')
  disp(' ')
  disp('Note:  The errors for slstly are evaluated for X = U''*U')
end
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
  
  A = rand(n,n)/n - n*eye(n);  
  X0 = rand(n,n);  X0 = X0 + X0' + 2*n*eye(n);
  C = -(A'*X0 + X0*A);  C = (C + C')/2.;  CC = chol(C);
  normX0 = max(1,norm(X0,'fro'));
  
% SLICOT calculations.

  time = cputime;    
  XC = slstly(A,CC);  x = XC'*XC;
  res(l,2) = cputime - time;  
  res(l,4) = norm(A'*x + x*A + C,'fro')/normX0;
  res(l,6) = norm(x - X0,'fro')/normX0;
  
% MATLAB calculations.
  
  time = cputime;    
  if ver7,
    Xc = lyapchol(A',CC');  X = Xc'*Xc;
  else
    X = lyap(A',C);
  end
  res(l,3) = cputime - time;  
  res(l,5) = norm(A'*X + X*A + C,'fro')/normX0;
  res(l,7) = norm(X - X0,'fro')/normX0;

end  

disp('------------------------------------------------------------------------')
disp('            Time         Relative residual          Relative error      ')   
disp('------------------------------------------------------------------------')
if ver7,
  disp('    n  slstly lyapchol  slstly      lyapchol      slstly      lyapchol  ')
else
  disp('    n  slstly  lyap     slstly        lyap        slstly        lyap    ')
end
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
page_plot(res,'slstly',1)
uiwait
disp(' ')
pause,

rp_slstly1 = res;
clear res space
space = repmat(sp,Ntests2,1);  res = zeros(Ntests2,7);

l = 0;
for Nt = sz2:sz2+Ntests2-1,
  n = 2^Nt;
  l = l + 1;
  res(l,1) = n;
  
  A = rand(n,n) - n*eye(n);  
  X0 = rand(n,n);  X0 = X0 + X0' + 2*n*eye(n);
  C = -(A'*X0 + X0*A);  C = (C + C')/2.;  CC = chol(C);
  normX0 = max(1,norm(X0,'fro'));
  
% SLICOT calculations.

  time = cputime;    
  XC = slstly(A,CC);  x = XC'*XC;
  res(l,2) = cputime - time;  
  res(l,4) = norm(A'*x + x*A + C,'fro')/normX0;
  res(l,6) = norm(x - X0,'fro')/normX0;
  
% MATLAB calculations.
  
  time = cputime;    
  if ver7,
    Xc = lyapchol(A',CC');  X = Xc'*Xc;
  else
    X = lyap(A',C);
  end
  res(l,3) = cputime - time;  
  res(l,5) = norm(A'*X + X*A + C,'fro')/normX0;
  res(l,7) = norm(X - X0,'fro')/normX0;

end  

disp('-------------------------------------------------------------------------')
disp('            Time         Relative residual          Relative error      ')   
disp('-------------------------------------------------------------------------')
if ver7,
  disp('    n  slstly lyapchol  slstly      lyapchol      slstly      lyapchol  ')
else
  disp('    n  slstly  lyap     slstly        lyap        slstly        lyap    ')
end
disp('-------------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,num2str(res(:,2),'%5.2f'), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4)), ...
  space,num2str(res(:,5)),space,num2str(res(:,6)), ...	
  space,num2str(res(:,7))] )	
disp('-------------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slstly',1)
uiwait
disp(' ')
pause,

rp_slstly2 = res;
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
  
  A = rand(n,n)/n - n*eye(n);  
  X0 = rand(n,n);  X0 = X0 + X0' + 2*n*eye(n);
  C = -(A'*X0 + X0*A);  C = (C + C')/2.;  CC = chol(C);
  normX0 = max(1,norm(X0,'fro'));
  
% SLICOT calculations.

  time = cputime;    
  XC = slstly(A,CC);  x = XC'*XC;
  res(l,2) = cputime - time;  
  res(l,4) = norm(A'*x + x*A + C,'fro')/normX0;
  res(l,6) = norm(x - X0,'fro')/normX0;
  
% MATLAB calculations.
  
  time = cputime;    
  if ver7,
    Xc = lyapchol(A',CC');  X = Xc'*Xc;
  else
    X = lyap(A',C);
  end
  res(l,3) = cputime - time;  
  res(l,5) = norm(A'*X + X*A + C,'fro')/normX0;
  res(l,7) = norm(X - X0,'fro')/normX0;

end  

disp('-------------------------------------------------------------------------')
disp('            Time         Relative residual          Relative error      ')   
disp('-------------------------------------------------------------------------')
if ver7,
  disp('    n  slstly lyapchol  slstly      lyapchol      slstly      lyapchol  ')
else
  disp('    n  slstly  lyap     slstly        lyap        slstly        lyap    ')
end
disp('-------------------------------------------------------------------------')
		                    
disp( [space,int2str(res(:,1)),space,num2str(res(:,2),'%5.2f'), ...
  space,num2str(res(:,3),'%5.2f'),space,num2str(res(:,4)), ...
  space,num2str(res(:,5)),space,num2str(res(:,6)), ...	
  space,num2str(res(:,7))] )	
disp('-------------------------------------------------------------------------')
		                    
%disp('Press any key to continue')
echo on 
echo off
page_plot(res,'slstly',1)
uiwait
disp(' ')
pause,

rp_slstly3 = res;
clear res A C X x CC X0 XC normX0
space = [];

perfb_slstly