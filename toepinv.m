function [Out1, Out2, Out3] = toepinv(T, job)
% TG denotes the first row of a Toeplitz matrix T. Then this function
% will compute
%
% I)  G = toepinv(TG)
%        or
%     G = toepinv(TG, 'G')
%     the generator G of T^(-1)
%
% II) [G, R] = toepinv(TG, 'R')
%     the generator G and the Cholesky factor R of T
%
% III)[G, LI] = toepinv(TG)
%        or
%     [G, LI] = toepinv(TG, 'L')
%     the generator G and the Cholesky factor LI of T^(-1)
%
% IV) [G, LI, R] = toepinv(TG)
%        or
%     [G, LI, R] = toepinv(TG, 'A')
%     the generator G, the lower Cholesky factor LI of T^(-1)
%     and the upper Cholesky factor R of T
%
% V)  R = toepinv(TG, 'O');
%     the upper Cholesky factor R of T.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Daniel Kressner, June 2000.
%
% Revised: V. Sima, Nov. 2000, Jan. 2001, Mar. 2009. 

% Computation options.
cgen = 0;
cli  = 0;
cr   = 0;

if nargin == 0,
  error('Not enough input arguments.');
end

if nargin < 2,
  if nargout <= 1,
    job = 'G';
  elseif nargout == 2,
    job = 'L';
  elseif nargout == 3,
    job = 'A';
  else
    error('Too many output arguments.');
  end  
end

if (job == 'G')
  cgen = 1;
elseif (job == 'L')
  cgen = 1;
  cli  = 1;
elseif (job == 'A')
  cgen = 1;
  cli  = 1;
  cr   = 1;
elseif (job == 'O')
  cr   = 1;
elseif (job == 'R')
  cgen = 1;
  cr   = 1;    
else
  error('Unknown job option.');
end
  

[k, m] = size(T);
n = fix(m/k);
if (n ~= m/k)
  error('T does not have block structure.');
end

% Bring T to proper form.
T(1:k, 1:k) = chol(T(1:k, 1:k));
T(1:k, (k+1):m) = T(1:k, 1:k)' \ T(1:k, (k+1):m);


if cr,
  R = zeros(m);
  R(1:k, :) = T;
end

if cli,
  LI = zeros(m);
end

if cgen,
  % We can use the generator as working array.
  GI = zeros(2*k,n*k);
  GI(k+1:2*k, 1:k) = T(1:k, 1:k)' \ eye(k);
  GI(k+1:2*k, k+1:n*k) = T(:, 1:(n-1)*k);
  GI(1:k, 1:k) = GI(k+1:2*k, 1:k);
  
  if cli,
    LI(1:k, 1:k) = GI(k+1:2*k, 1:k);
  end  

  for i = 2:n,
    startr = (i-1)*k+1;
    endr = i*k;

    % Do the normal stuff.
    [GI(k+1:2*k, k+1:2*k), T(:, startr:endr), rho, q] =  ...
       createg(GI(k+1:2*k, k+1:2*k), T(:, startr:endr));
    [GI(k+1:2*k, 2*k+1:m-(i-1)*k+k), T(:, endr+1:m)] =  ...
       applyg(GI(k+1:2*k, 2*k+1:m-(i-1)*k+k), T(:, endr+1:m), T(:, startr:endr), rho, q);
    if cr,
      R(startr:endr, startr:m) = GI(k+1:2*k, k+1:m-(i-1)*k+k);
    end
    
    % Now, do the inverse stuff.
    GI(k+1:2*k, m-i*k+1+k:m-(i-1)*k+k) = zeros(k);
    % The block we swapped to the beginning,
    [GI(k+1:2*k, 1:k), GI(1:k, (i-1)*k+1:i*k)] = ...
      applyg(GI(k+1:2*k, 1:k), GI(1:k, (i-1)*k+1:i*k), T(:, startr:endr), rho, q);
    % and the rest.
    [GI(k+1:2*k, m-i*k+1+k:m), GI(1:k, 1:(i-1)*k)] = ...
      applyg(GI(k+1:2*k, m-i*k+1+k:m), GI(1:k, 1:(i-1)*k), T(:, startr:endr), rho, q);
    if cli,
      LI(startr:endr, 1:i*k) = [GI(k+1:2*k, m-i*k+1+k:m), GI(k+1:2*k, 1:k)];
    end   
  end
  
  GI(k+1:2*k, 1:k) = zeros(k);
  Out1 = GI;
  if (job == 'R')
    Out2 = R;
  end
  if (job == 'L')
    Out2 = LI;
  end
  if (job == 'A')
    Out2 = LI;
    Out3 = R;
  end
else
  % We have to use the R-factor as working array.
  R(k+1:2*k, k+1:m) = T(:, 1:(n-1)*k);
  
  for i = 2:n,
    startr = (i-1)*k+1;
    endr = i*k;
    [R(startr:endr, startr:endr), T(:, startr:endr), rho, q] =  ...
       createg(R(startr:endr, startr:endr), T(:, startr:endr));
    [R(startr:endr, endr+1:m), T(:, endr+1:m)] =  ...
       applyg(R(startr:endr, endr+1:m), T(:, endr+1:m), T(:, startr:endr), rho, q);
    if i<n,
      R(i*k+1:(i+1)*k, i*k+1:m) = R(startr:endr, startr:m-k);
    end
  end
  Out1 = R;
end
