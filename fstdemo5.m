echo on
%   This demonstration shows the performance of the fast Toeplitz solver 
%   incorporated in the SLICOT Toolbox (fstoep), in comparison with
%   standard Matlab functions and the corresponding fast Matlab solver,
%   developed under the NICONET project. 
%   Random data sets are used and the relative errors or relative 
%   residuals are checked.  Positivity is forced.  The coefficient matrices
%   are nxk-by-nxk symmetric block Toeplitz matrices, with k taking values 
%   in [1 2 20 30], and n obtained by dividing 300 to k.  
%

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

echo off

disp('Random test of fast Toeplitz solvers')
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

if ~exist('tol', 'var'),  tol = 1.e-11;  end

if ~exist('Details', 'var'),  Details = 0;  end

total_time = cputime;

count = 0;  l = 0;  lx = 0;

max_err = 0.0;  max_err1 = 0.0;  max_err2 = 0.0;

dim = zeros(1,3);  dimx = zeros(1,4);

listk = [1 2 20 30];  listn = 300./listk;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Main loop for computations.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

nt = numel(listk);
time1 = zeros(nt,4);  err1 = zeros(nt,3);
time2 = zeros(nt,4);  err2 = zeros(nt,4);

timex1 = zeros(nt^2,3);  errx1 = zeros(nt^2,4);
timex2 = zeros(nt^2,3);  errx2 = zeros(nt^2,4);

for k = listk
   if k == listk(2),  echo off,  end
   l = l + 1;  n = listn(l);  dim(l,:) = [l k n];
        
   %  Factorization performance tests (MB02CD-based and Matlab-based).
        
   Tr = rand(k,n*k);   
   if n > 0,  Tr(1:k,1:k) = Tr(1:k,1:k) + Tr(1:k,1:k)' + k*n*eye(k);  end   
   Tc = Tr';  BT = btoeplitz(Tc,Tr);
         
   %  Check with TYPET = 'R'.
   %  Check that different options produce the same results.
    
   time           = cputime;    
   R              = fstoep(1,Tr);
   time1(l,1)     = cputime - time;  
   time           = cputime;    
   Rm             = toepinv(Tr,'O');
   time2(l,1)     = cputime - time;  
   time           = cputime;    
   G              = fstoep(2,Tr);
   time1(l,2)     = cputime - time;  
   time           = cputime;    
   Gm             = toepinv(Tr,'G');
   time2(l,2)     = cputime - time;  
   time           = cputime;    
   [G1,R1]        = fstoep(3,Tr);
   time1(l,3)     = cputime - time;  
   time           = cputime;    
   [Gm1,Rm1]      = toepinv(Tr,'R');
   time2(l,3)     = cputime - time;  
   
   Rm  = Rm(1:n*k,1:n*k);
   err = max( [ norm(R  - R1,  1), norm(G  - G1,  1), ...
                norm(Rm - Rm1, 1), norm(Gm - Gm1, 1) ] );
   clear G1 Gm1 R1 Rm1
          
   time           = cputime;    
   [G2,Li]        = fstoep(4,Tr);
   time1(l,4)     = cputime - time;  
   time           = cputime;    
   [G3,Li1,R2]    = fstoep(5,Tr);
   time1(l,5)     = cputime - time;  
   
   err = max( [ norm(R - R2, 1), norm(G  - G2,  1), ...
                norm(G - G3, 1), norm(Li - Li1, 1), err ] );
   clear G2 G3 Li1 R2
          
   time           = cputime;    
   [Gm2,Lim]      = toepinv(Tr,'L');
   time2(l,4)     = cputime - time;  
   time           = cputime;    
   [Gm3,Lim1,Rm2] = toepinv(Tr,'A');
   time2(l,5)     = cputime - time;  
   
   err = max( [ norm(Rm - Rm2, 1), norm(Gm - Gm2, 1), ...
                norm(Gm - Gm3, 1), norm(Lim - Lim1, 1), err ] );
   clear Gm2 Gm3 Lim1 Rm2
          
   if ~( err == 0 ),       
      max_err2 = max( max_err2, 1/eps );
   end
        
   %  Check the factors.
        
   if n*k > 0,
      approxinv = btoeplitz(G(1:k,:)', G(1:k,1:k) * eye(k,size(G,2))) ...
                * btoeplitz(G(1:k,:)', G(1:k,1:k) * eye(k,size(G,2)))' ...
                - btoeplitz(G(k+1:2*k,:)', G(k+1:2*k,1:k) * eye(k,size(G,2))) ...
                * btoeplitz(G(k+1:2*k,:)', G(k+1:2*k,1:k) * eye(k,size(G,2)))';
      nrmT = norm(BT,1);
      nrmG = norm(G, 1);
      nrmR = norm(R, 1);
      nrmL = norm(Li,1);
   else
      approxinv = [];
      nrmT = 1;
      nrmG = 1;
      nrmR = 1;
      nrmL = 1;
   end
   
   err1(l,1:3) = [ norm(R'*R - BT,1) / nrmT, ...
                   norm(Li * BT * Li'  - eye(n*k),1), ...
                   norm(approxinv * BT - eye(n*k),1) ];
   clear approxinv
        
   if n*k > 0,
      approxinm = btoeplitz(Gm(1:k,:)', Gm(1:k,1:k) * eye(k,size(Gm,2))) ...
                * btoeplitz(Gm(1:k,:)', Gm(1:k,1:k) * eye(k,size(Gm,2)))' ...
                - btoeplitz(Gm(k+1:2*k,:)', Gm(k+1:2*k,1:k) * eye(k,size(Gm,2))) ...
                * btoeplitz(Gm(k+1:2*k,:)', Gm(k+1:2*k,1:k) * eye(k,size(Gm,2)))';
   else
      approxinm = [];
   end
   err2(l,1:3) = [ norm(Rm'*Rm - BT,1) / nrmT, ...
                   norm(Lim * BT * Lim' - eye(n*k),1), ...
                   norm(approxinm * BT  - eye(n*k),1) ];
   clear approxinm

   Rc        = chol(BT);
   err2(l,4) = norm(Rc'*Rc - BT,1) / nrmT;
   clear Rc

   max_err1 = max( [ max_err1, err1(l,1:3) ] );
   max_err2 = max( [ max_err2, err2(l,1:4) ] );
   max_err  = max( [ max_err, norm(R - Rm,1) / nrmR, ...
                     norm(Li - Lim,1) / nrmL, ...
                     norm(G - Gm,1) / nrmG ] );
   clear G Gm Li Lim
        
   %  Update the counter.
        
   count = count + 1;
    
   %  Factorization and/or solution tests (MB02CD, MB02ED).
        
   for nrhs = listk
      if nrhs == listk(2),  echo off,  end
           
      %  Check with TYPET = 'R'.
      %  Check that some options produce the same results.
      %  (Not possible for X.)
           
      lx = lx + 1;
      dimx(lx,:) = [lx k n nrhs];
           
      Xr = rand(nrhs,n*k);  Xc = Xr';
      Br = Xr*BT;  Bc = Br';
      if n*k*nrhs > 0,
         nrmB = norm(Br);
         nrmX = norm(Xr);
      else
         nrmB = 1;
         nrmX = 1;
      end
      time           = cputime;    
      [R,X]          = fstoep(1,Tr,Br);
      timex1(lx,1)   = cputime - time;  
      time           = cputime;    
      Rm             = chol(BT);
      Xm             = Br/BT;
      timex2(lx,1)   = cputime - time;  
      time           = cputime;    
      X4             = fstoep(11,Tr,Br);
      timex1(lx,2)   = cputime - time;  
      time           = cputime;    
      Xm4            = Br/BT;
      timex2(lx,2)   = cputime - time;  
      max_err = max( [ max_err, norm(R - Rm, 1) / nrmR, ...
                       norm(X - Xm) / nrmX, norm(X - X4) / nrmX ] );
           
      %  Check the solution.
           
      errx1(lx,1:2) = [ norm(X*BT  - Br) / nrmB, norm(X  - Xr) / nrmX ];
      errx2(lx,1:2) = [ norm(Xm*BT - Br) / nrmB, norm(Xm - Xr) / nrmX ];
      max_err1 = max( [ max_err1, errx1(lx,1:2) ] );
      max_err2 = max( [ max_err2, errx2(lx,1:2) ] );
      max_err  = max( [ max_err, norm(X - Xm) / nrmX, ...
                        norm(R - Rm,1) / nrmR ] );
          
      if ~( n == 1 || n*k*nrhs == 0 ), 
              
         %  Check with TYPET = 'C'.
              
         time         = cputime;    
         X            = fstoep(11,Tc,Bc);
         timex1(lx,3) = cputime - time;  
         time         = cputime;    
         Xm           = BT\Bc;
         timex2(lx,3) = cputime - time;  
         
         errx1(lx,3:4) = [ norm(BT*X  - Bc) / nrmB, norm(X  - Xc) / nrmX ];
         errx2(lx,3:4) = [ norm(BT*Xm - Bc) / nrmB, norm(Xm - Xc) / nrmX ];
         max_err1 = max( [ max_err1, errx1(lx,3:4) ] );
         max_err2 = max( [ max_err2, errx2(lx,3:4) ] );
         max_err  = max( max_err, norm(X - Xm) / nrmX );
              
         %  Update the counter.
              
         count = count + 1;
      end              
      
      %  Update the counter.
           
      count = count + 1;
   end              
end
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' ')
if max_err < tol
   disp(['fstoep :    passed  -- relative error norm = ', num2str(max_err)])
   disp(['            Number of problems solved  = ', num2str(count)])
else
   disp(['fstoep :    failed  -- relative error norm = ', num2str(max_err)])
end
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Assess the results.

nrme1 = [ norm(err1(:,1)) norm(err1(:,2)) norm(err1(:,3)) ];
nrme2 = [ norm(err2(:,1)) norm(err2(:,2)) norm(err2(:,3)) norm(err2(:,4)) ];

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' ')
disp('Cumulative relative errors in factorizations, norm(e_?(:,j))')
disp(' ')
disp('Legend:  e_R = norm(R''*R - BT,1)/norm(BT,1);')
disp('         e_L = norm(L*BT*L'' - eye(size(BT,1)),1);')
disp('         e_I = norm(BTi*BT  - eye(size(BT,1)),1);  with BTi - approximate inverse')
disp('         e_C = norm(U''*U - BT,1)/norm(BT,1);       with U = chol(BT);')
disp(' ')
disp('-------------------------------------------------------------------------------------------')
disp('                        fstoep                                 Matlab                      ')
disp('-------------------------------------------------------------------------------------------')
disp('      e_R          e_L          e_I          e_R          e_L          e_I          e_C    ')
disp('-------------------------------------------------------------------------------------------')
disp([sprintf('  %0.4e',nrme1),sprintf('  %0.4e',nrme2)])
disp('-------------------------------------------------------------------------------------------')
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

nrmex1 = [ norm(errx1(:,1)) norm(errx1(:,2)) norm(errx1(:,3)) norm(errx1(:,4)) ];
nrmex2 = [ norm(errx2(:,1)) norm(errx2(:,2)) ];

echo off
disp(' ')
disp('Cumulative relative errors in solutions, norm(res_X(:,j)), norm(e_X(:,j)), j = 1 : 4')
disp(' ')
disp('Legend:  res_X = norm(Y*BT - C) / norm(C) or res_X = norm(BT*X - B) / norm(B);')
disp('           e_X = norm(X - X_true) / norm(X_true);')
disp('           C = B'',  Y = X''. Slightly different algorithms are used for X and Y.')
disp(' ')
disp('------------------------------------------------------------------------------')
disp('                        fstoep                                 Matlab         ')
disp('      Y*BT = C ([R,Y])           BT*X = B ([X])         Y*BT = C,  BT*X = B   ')
disp('------------------------------------------------------------------------------')
disp('     res_X         e_X         res_X         e_X         res_X         e_X    ')
disp('------------------------------------------------------------------------------')
disp([sprintf('  %0.4e',nrmex1),sprintf('  %0.4e',nrmex2)])
disp('------------------------------------------------------------------------------')
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' ')
disp('Cumulative execution times')
disp(' ')
disp('                               R        G      [G,R]    [G,L]   [G,L,R]' )
disp(['fstoep - factorization : ', sprintf('%9.2f',sum(time1))])
disp(['Matlab - factorization : ', sprintf('%9.2f',sum(time2))])
disp(' ')
disp('                              [R,Y=C*inv(BT)]  Y=C*inv(BT)   X=inv(BT)*B' )
disp(['fstoep - solution      : ', sprintf('%15.2f',sum(timex1))])
disp(['Matlab - solution      : ', sprintf('%15.2f',sum(timex2))])
disp(' ')
disp('NOTE: The Matlab function inv is never actually used in this demo.')
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Perturb 0 values in the reference timing to 0.01, to avoid dividing by 0.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

time1m = time1;
indx   = find( time1m == 0 );
time1m(indx) = 0.01;

timex1m = timex1;
indxm   = find( timex1m == 0 );
timex1m(indxm) = 0.01;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' ')
disp('Speed-up factors: fstoep versus Matlab')
disp(' ')
disp('                               R        G      [G,R]    [G,L]   [G,L,R]' )
disp(['Factorization          : ', sprintf('%9.2f',sum(time2)./sum(time1m))])
disp(' ')
disp('                               [R,Y=C*inv(BT)]  Y=C*inv(BT)   X=inv(BT)*B' )
disp(['Solution               : ', sprintf('%15.2f',sum(timex2)./sum(timex1m))])
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot detailed results.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
set(axes,'FontSize',12)
bar([err1 err2])
title('Relative errors in factorizations, norm(e_?(:,j))')
legend('e_R(fstoep)','e_L(fstoep)','e_I(fstoep)',...
       'e_R(toepinv)','e_L(toepinv)','e_I(toepinv)','e_C(chol)')
xlabel('Problem # (determined by k and n)')
ylabel('Relative errors')

shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end
close(gcf)
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
set(axes,'FontSize',12)
bar([errx1(:,1:2) errx2(:,1:2)])
title('Relative errors/residuals in solutions, Y*BT = C, ([R,Y])')
legend('res_X(fstoep)','e_X(fstoep)','res_X(chol)','e_X(chol)')
xlabel('Problem # (determined by k, n and nrhs)')
ylabel('Relative errors')

shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end
close(gcf)
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
set(axes,'FontSize',12)
bar([errx1(:,3:4) errx2(:,3:4)])
title('Relative errors/residuals in solutions, BT*X = B ([X])')
legend('res_X(fstoep)','e_X(fstoep)','res_X(chol)','e_X(chol)')
xlabel('Problem # (determined by k, n and nrhs)')
ylabel('Relative errors')

shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end
close(gcf)
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
set(axes,'FontSize',12)
bar([time1 time2])
title('Execution times for factorization')
legend('R(fstoep)','G(fstoep)','[G,R](fstoep)','[G,L](fstoep)','[G,L,R](fstoep)',...
       'R(toepinv)','G(toepinv)','[G,R](toepinv)','[G,L](toepinv)','[G,L,R](toepinv)')
xlabel('Problem # (determined by k and n)')
ylabel('Execution times')

shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end
close(gcf)
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
set(axes,'FontSize',12)
bar([timex1 timex2])
title('Execution times for solution')
legend('[R,Y=C*inv(BT)](fstoep)','Y=C*inv(BT)(fstoep)','X=inv(BT)*B(fstoep)',...
       '[R,Y=C*inv(BT)](chol)','Y=C*inv(BT)(chol)','X=inv(BT)*B(chol)')
xlabel('Problem # (determined by k, n and nrhs)')
ylabel('Execution times')

shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end
close(gcf)
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
set(axes,'FontSize',12)
bar(time2./time1m)
title('Speed-up factors for factorization: fstoep versus Matlab toepinv')
legend('R','G','[G,R]','[G,L]','[G,L,R]')
xlabel('Problem # (determined by k and n)')
ylabel('Speed-up factors fstoep/toepinv')

shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end
close(gcf)
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
set(axes,'FontSize',12)
bar(timex2./timex1m)
title('Speed-up factors for solution: fstoep versus Matlab')
legend('[R,Y=C*inv(BT)]','Y=C*inv(BT)','X=inv(BT)*B')
xlabel('Problem # (determined by k, n and nrhs)')
ylabel('Speed-up factors fstoep/Matlab')

shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end
close(gcf)
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

if Details 
   format short e
   echo off
   disp(' ')
   disp('Relative errors in the factors computed by fstoep and Matlab toepinv')
   disp(' ')
   disp('Legend:  e_R = norm(R''*R - BT,1)/norm(BT,1);')
   disp('         e_L = norm(L*BT*L'' - eye(size(BT,1)),1);')
   disp('         e_I = norm(BTi*BT  - eye(size(BT,1)),1);  with BTi - approximate inverse')
   disp('         e_C = norm(U''*U - BT,1)/norm(BT,1);       with U = chol(BT);')
   disp(' ')
   
   disp('-----------------------------------------------------------------------------------------------')
   disp('    Ex.    Dim.                                 Relative errors                                ')
   disp('                                fstoep                                 Matlab                  ')
   disp('-----------------------------------------------------------------------------------------------')
   disp('     #     k    n      e_R        e_L        e_I        e_R        e_L        e_I        e_C   ')
   disp('-----------------------------------------------------------------------------------------------')
   for L = 1 : l
       disp([sprintf('%6d',dim(L,:)),sprintf('  %0.2e',err1(L,1:3)),sprintf('  %0.2e',err2(L,1:4))])
   end   
   disp('-----------------------------------------------------------------------------------------------')
   echo on

   if pause_wait < 0,  pause,  else  pause(pause_wait),  end

   echo off
   disp(' ')
   disp('Relative errors in the solutions computed by fstoep and Matlab')
   disp(' ')
   disp('Legend:  res_X = norm(Y*BT - C) / norm(C) or res_X = norm(BT*X - B) / norm(B);')
   disp('           e_X = norm(X - X_true) / norm(X_true);')
   disp('           C = B'',  Y = X''. Slightly different algorithms are used for X and Y.')
   disp(' ')
   
   disp('----------------------------------------------------------------------------------------------')
   disp('   Ex.       Dimensions                               Relative errors                         ')
   disp('                                              fstoep                           Matlab         ')
   disp('                              Y*BT = C ([R,Y])       BT*X = B ([X])      Y*BT = C, BT*X = B   ')
   disp('----------------------------------------------------------------------------------------------')
   disp('    #      k     n    nrhs    res_X       e_X       res_X       e_X       res_X       e_X     ')
   disp('----------------------------------------------------------------------------------------------')
   for L = 1 : lx
       disp([sprintf('%6d',dimx(L,:)),'  ',sprintf('  %0.2e',errx1(L,1:4)),sprintf('  %0.2e',errx2(L,1:2))])
   end   
   disp('----------------------------------------------------------------------------------------------')
   echo on

   if pause_wait < 0,  pause,  else  pause(pause_wait),  end

   echo off
   disp(' ')
   disp('Execution times for factorizations computed by fstoep and Matlab toepinv')
   disp(' ')
   
   disp('-------------------------------------------------------------------------------------------------------------')
   disp('    Ex.     Dim.                                         Execution times                                     ')
   disp('                                       fstoep                                       toepinv                  ')
   disp('-------------------------------------------------------------------------------------------------------------')
   disp('     #     k    n       R        G      [G,R]    [G,L]   [G,L,R]     R        G      [G,R]    [G,L]   [G,L,R]')
   disp('-------------------------------------------------------------------------------------------------------------')
   for L = 1 : l
       disp([sprintf('%6d',dim(L,:)), sprintf('%9.2f',time1(L,:)), sprintf('%9.2f',time2(L,:))])
   end   
   disp('-------------------------------------------------------------------------------------------------------------')   
   echo on

   if pause_wait < 0,  pause,  else  pause(pause_wait),  end

   echo off
   disp(' ')
   disp('Speed-up factors for factorization: fstoep versus Matlab toepinv')
   disp(' ')
   disp('-----------------------------------------------------------------')   
   disp('     #     k    n       R        G      [G,R]    [G,L]   [G,L,R] ')
   disp('-----------------------------------------------------------------')   
   for L = 1 : l
       disp([sprintf('%6d',dim(L,:)), sprintf('%9.2f',time2(L,:)./time1m(L,:))])
   end   
   disp('-----------------------------------------------------------------')   
   disp(' ')
   echo on

   if pause_wait < 0,  pause,  else  pause(pause_wait),  end

   echo off
   disp(' ')
   disp('Execution times for solutions computed by fstoep and Matlab')
   disp(' ')
   
   disp('-----------------------------------------------------------------------------------------------------')
   disp('    Ex.      Dimensions                                 Execution times                              ')
   disp('                                          fstoep                               Matlab                ')
   disp('-----------------------------------------------------------------------------------------------------')
   disp('     #     k    n    nrhs      [R,X]    Y=C*inv(BT)  X=inv(BT)*B   [R,X]    Y=C*inv(BT)  X=inv(BT)*B ')
   disp('-----------------------------------------------------------------------------------------------------')
   for L = 1 : lx
       disp([sprintf('%6d',dimx(L,:)), sprintf('%12.2f',timex1(L,:)), sprintf('%12.2f',timex2(L,:))])
   end   
   disp('-----------------------------------------------------------------------------------------------------')
   echo on

   if pause_wait < 0,  pause,  else  pause(pause_wait),  end

   echo off
   disp(' ')
   disp('Speed-up factors for solution: fstoep versus Matlab')
   disp(' ')
   disp('-----------------------------------------------------------------')   
   disp('     #     k    n    nrhs      [R,X]    Y=C*inv(BT)  X=inv(BT)*B ')
   disp('-----------------------------------------------------------------')   
   for L = 1 : lx
       disp([sprintf('%6d',dimx(L,:)), sprintf('%12.2f',timex2(L,:)./timex1m(L,:))])
   end   
   disp('-----------------------------------------------------------------')   
   disp(' ')
end
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

total_time = cputime - total_time; 

echo off
disp(['total_time = ', num2str(total_time)])
