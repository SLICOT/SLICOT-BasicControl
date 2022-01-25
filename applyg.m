function [A1, B1] = applyg(A1, B1, B, rho, qt)
% [A1, B1] = applyg(A1, B1, B, rho, qt);
%
% input:
%
%    A1 - real p x n matrix
%    B1 - real q x n matrix
%     B - stored Householder vectors from createg, q x k
%   rho - stored hyperbolic rotations from createg, p x 1
%    qt - orthogonal transformation matrix from createg
%
% output:
%
%    A1 - altered A1
%    B1 - altered B1
%
% purpose:
%
%    To apply the transformations obtained from createg
%    to further columns of the generator matrix, denoted
%    by A1 and B1.

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Daniel Kressner, June 2000.
%
% Revised: V. Sima, Nov. 2000, Jan. 2001, Mar. 2009. 

% Input processing.
if nargin < 5,
  error('Not enough input arguments.');
end

[p,n] = size(A1);
[q,k] = size(B);

if n <= 0,
  return;
end

if size(B1,2) ~= n,
  error('A1 and B1 must have the same number of columns.');
end
if length(rho) ~= p,
  error('rho must have the length equal to the number of rows in A1.');
end
if p < k,
  warning('SLICOT:applyg', ...
          'The number of rows in A1 must be greater than or equal to k.');
end

% Applying the transformations.
B1 = qt' * B1;
for i = 1:k,
  if i > 1,
    v = B(1:min(i,q), i);
    x = v' * B1(1:min(i,q),:);
    B1(1:min(i,q),:) = B1(1:min(i,q),:) - v*x;
  end
  
  if i <= p,
    r = rho(i);
    s = sqrt(1-abs(r)^2); c = r;
    A1(i,:) = ([ 1 -conj(c) ]/s)*[A1(i,:); B1(1,:)];
    B1(1,:) = ([-c s])*[A1(i,:); B1(1,:)];
  end
end

return
