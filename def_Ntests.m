% Script for defining default values for the number of tests to be
% performed for testing the performance of SLICOT routines.
% Random data sets of row/column dimensions multiplies of 10, or powers
% of 2, but not larger than 128, are used.
% The number of tests performed is for each group of tests is defined 
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
%   V. Sima, May 2005, Jan. 2007, Mar. 2009.
%   

global sz sz2

if ~exist( 'sz', 'var' ) || isempty( sz ),
  % Choose resonable dimensions for matrices.
  sz = 100;
  while (1)
    A = rand(sz);    B = rand(sz);  
    time = cputime;  C = A*B;  timing = cputime - time;  
     if timing > 0,  break,  end
  end
  sz = sz/20;  sz2 = fix(log2(sz*10));
  % Dummy calls for loading the solvers.
  A = rand(sz);  E = rand(sz);  X0 = rand(sz);  X0 = X0 + X0';
  C = A'*X0 + X0*A;   C = (C + C')/2;
  sllyap(A,C);  lyap(A',-C);  dlyap(A',-C);
  C = A'*X0*E + E'*X0*A;  C = (C + C')/2;
  slgely(A,E,C);  x = lyap(A,E,C);
end

if ~exist('Ntests10','var')
  Ntests10 = 5;
elseif Ntests10 > 20  
  disp('Warning: Number of tests (Ntests10) to perform for each ')
  disp(['         "linear" group of tests is ',int2str(Ntests10),'.'])
  disp('The computation time could be high.')
  pause
  disp(' ')
end    

if ~exist('Ntests2','var')
  Ntests2 = 4; 
elseif Ntests2 > 4
  disp('Warning: Number of tests (Ntests2) to perform for each ')
  disp(['         "exponential" group of tests is ',int2str(Ntests2),'.'])
  disp('The computation time could be high.')
  pause
  disp(' ')
end    

