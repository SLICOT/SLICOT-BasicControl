%   The SLICOT Toolbox for Structured Matrices currently contains 
%   algorithms for factoring symmetric positive definite (block)  
%   Toeplitz matrices and/or solving associated linear systems.
%
%   SLICOT Toolbox for Structured Matrices demonstrations:
%
%   1) Compute the Cholesky factor for block Toeplitz matrices.
%   2) Compute the generator and the Cholesky factor for the inverse
%      of a block Toeplitz matrix.
%   3) Solve a linear system with a block Toeplitz matrix.
%   4) Update the factorizations after adding new blocks in a
%      Toeplitz matrix.
%   5) Compare SLICOT performance with equivalent Matlab algorithms.
%
%   0) Quit.

%   RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   V. Sima 30-12-2000.
%
%   Revisions:
%   V. Sima 03-07-2020.
%

global pause_wait  % This could be used in pause(n) command.
                   % If pause_wait < 0, standard command pause is used (default).
                   % Any key should then be pressed to continue.
                   % It may need to use the command "pause on" before
                   % calling fstdemo.

K = 0;
while (1)
   disp(' ')
   help fstdemo
   K = input('Select a demonstration example number: ');
   disp(' ')
   if isempty(K), K = 20; end
   if K == 0,  break,     end
   if K == 1,  fstdemo1,  end
   if K == 2,  fstdemo2,  end
   if K == 3,  fstdemo3,  end
   if K == 4,  fstdemo4,  end
   if K == 5,  fstdemo5,  end
end
