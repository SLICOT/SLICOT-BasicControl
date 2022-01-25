function [A, B, rho, qt] = createg(A, B)
% [A, B, rho] = createg(A, B)
%
% input:
%
%    A - real upper triangular p x k matrix
%    B - real q x k matrix
%
% output:
%   
%    A - altered, but still upper triangular
%    B - contains the Householder vectors
%  rho - contains the information for the hyperbolic rotations
%   qt - orthogonal transformation matrix
%
% purpose:
%
%    To reduce A and B as the leading blocks of a generator via symplectic
%    transformations to proper form. Information about the transformation
%    will be stored in B and rho.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Daniel Kressner, June 2000.
%
% Revised: V. Sima, Nov. 2000, Jan. 2001, Mar. 2009. 

% Input processing.

if nargin < 2,
  error('Not enough input arguments.');
end

[p,k] = size(A);
q     = size(B, 1);
rho   = [];

if size(B,2) ~= k,
  error('A and B must have the same number of columns.');
end
if p < k,
  warning('SLICOT:createg', ...
      'The number of rows in A should be greater than or equal to the number of columns.');
end
if (k==0)||(q==0)
  return
end

% Reduce B to upper triangular form.
[qt, B] = qr(B);

rho = zeros(1,min(k,p));

% Further reduction.
for i = 1:k,
  
  % Step 1: annihilate the upper rest of the i-th column of B
  %         with a Householder transformation.
  if i > 1,
    [v, beta] = house(B(1:min(i,q),i));
    v = v * sqrt(beta);
    x = v' * B(1:min(i,q), i:k);
    B(1:min(i,q), i:k) = B(1:min(i,q), i:k) - v*x;
  end
  
  % Step 2: annihilate the top entry in the i-th column of B
  %         with a hyperbolic rotation.
  if i <= p,
    if abs(B(1,i)) >= abs(A(i,i)),
      error('|Rho| >= 1, (numerically semi/indefiniteness)');
    end;
    r =  B(1,i) / A(i,i);
    rho(i) = r;	
    s = sqrt(1-abs(r)^2); c = r;
    A(i, i:k) = ([ 1 -conj(c) ]/s)*[A(i, i:k); B(1, i:k)];
    B(1, i:k) = ([-c s])*[A(i, i:k); B(1, i:k)];
  end  
  if i > 1,
    B(1:i,i) = v;
  end
end % for i

return
