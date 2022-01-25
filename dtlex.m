function [E,A,Y,B,X,U] = dtlex(nr,parin)
%DTLEX
%
% Usage:  [E,A,Y,B,X,U] = dtlex(nr,parin)
%         [E,A,Y,B,X,U] = dtlex(nr)
%
% Main routine of the benchmark library DTLEX (Version 1.0) described 
% in [1]. It generates benchmark examples of (generalized) discrete-time 
% Lyapunov equations
%
%        T           T 
%       A  X  A  -  E  X E  =  Y .                                   (1)
%
% In some examples, the right hand side has the form
%
%                T
%       Y  =  - B  B
%
% and the solution can be represented as a product of Cholesky factors
%
%              T 
%       X  =  U  U .
%
% E, A, Y, X, and U are real n-by-n matrices, and B is m-by-n. Note 
% that E can be the identity matrix. For some examples, B, X, or U are 
% not provided.
%
% Input:
%  - nr    : index of the desired example according to [1]; 
%            nr is a 1-by-2 matrix;
%            nr(1) defines the group:
%             = 1 : parameter-free problems of fixed size
%             = 2 : parameter-dependent problems of fixed size
%             = 3 : parameter-free problems of scalable size
%             = 4 : parameter-dependent problems of scalable size
%            nr(2) defines the number of the benchmark example within
%            a certain group.
%  - parin : parameters of the chosen example; 
%            referring to [1], the entries in parin have the following 
%            meaning:
%            Ex. 4.1 : parin(1:3) = [n r s]
%            Ex. 4.2 : parin(1:3) = [n lambda s]
%            Ex. 4.3 : parin(1:2) = [n t]
%            Ex. 4.4 : parin(1:2) = [q t]
%            parin is optional; default values as defined in [1] are
%            used as example parameters if parin is ommited. Note that
%            parin is not referenced if nr(1) = 1.
%
% Output:
%  - E, A, Y, B, X, U :  matrices of the Lyapunov equation (1).
%
% References: 
%
% [1]   D. Kressner, V. Mehrmann, and T. Penzl.
%       DTLEX - a Collection of Benchmark Examples for Discrete-Time 
%       Lyapunov Equations.
%       SLICOT Working Note 1999-7, 1999.

%  RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%  Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%  D. Kressner, V. Mehrmann, and T. Penzl (TU Chemnitz).
%  Feb 1, 1999.
%
%  Revisions:
%  V. Sima, Mar. 2, 2009.

E = []; A = []; Y = []; B = []; X = []; U = [];

if nargin < 1,
  error('Not enough input arguments.');
end;

if length(nr) < 2,
  error('Please use the nr = [group, example] notation.');
end;

if nr(1) == 4,
  if nr(2) == 1,
    %  Example 4.1:
    if nargin < 2,
      parin = [10 1.5 1.5];
    else
      if length(parin) < 3, error('Not enough input parameters.'); end
    end
    n = parin(1);
    r = parin(2);
    s = parin(3);
    if n ~= round(n) || n < 2,  error('Invalid value of parameter n.'); end
    if r <= 1,  error('Invalid value of parameter r.'); end
    if s <= 1,  error('Invalid value of parameter s.'); end

    E  = eye(n);  A = zeros( 1,n );  S = A;  SI = A;  f = A;
    for i = 1:n,
      A(i)  = (r^(i-1) - 1) / (r^(i-1) + 1);
      S(i)  =  s^(i-1);
      SI(i) =  s^(1-i);
      f(i)  = (-1)^(i-1);
    end;
    H1 = eye(n) - 2/n * ones(n);
    H2 = eye(n) - 2/n * f'*f;
    A  = H2 * diag(S) * H1 * diag(A) * H1 * diag(SI) * H2;
    B  = eye(1,n) * H1 * diag(SI) * H2;
    Y  = -B' * B;
    X  = -Y;

  elseif nr(2) == 2
    %  Example 4.2:
    if nargin < 2,
      parin = [10 -.5 1.5];
    else
      if length(parin) < 3, error('Not enough input parameters.'); end
    end
    n      = parin(1);
    lambda = parin(2);
    s      = parin(3);
    if n ~= round(n) || n < 2,  error('Invalid value of parameter n.'); end
    if (abs(lambda) >= 1),  error('Invalid value of parameter lambda.'); end
    if s <= 1,  error('Invalid value of parameter s.'); end

    E = eye(n);  S = zeros( 1,n );  SI = S;  f = S;
    A = lambda * eye(n) + diag(ones(n-1,1), 1);
    for i = 1:n,
      S(i)  =  s^(i-1);
      SI(i) =  s^(1-i);
      f(i)  = (-1)^(i-1);
    end;
    H1 = eye(n) - 2/n * ones(n);
    H2 = eye(n) - 2/n * f'*f;
    A  = H2 * diag(S) * H1 * A * H1 * diag(SI) * H2;
    B  = eye(1,n) * H1 * diag(SI) * H2;
    Y  = -B'*B;

  elseif nr(2) == 3
    %  Example 4.3:
    if nargin < 2,
      parin = [10 10];
    else
      if length(parin) < 2, error('Not enough input parameters.'); end
    end
    n = parin(1);
    t = parin(2);
    if n ~= round(n) || n < 2,  error('Invalid value of parameter n.'); end
    if t < 0,  error('Invalid value of parameter t.'); end

    E = eye(n) + 2^(-t) * tril(ones(n), -1);
    A = diag(2^(-t) * ones(1,n) +  (1:n)) + triu(ones(n), 1);
    X = ones(n);
    Y = A'*X*A - E'*X*E;

  elseif nr(2) == 4
    %  Example 4.4:
    if nargin < 2,
      parin = [10 1.5];
    else
      if length(parin) < 2, error('Not enough input parameters.'); end
    end
    q = parin(1);
    t = parin(2);
    n = q*3;
    if q ~= round(q) || q < 1,  error('Invalid value of parameter q.'); end
    if t < 1,  error('Invalid value of parameter t.'); end

    E = fliplr(tril(ones(n))) * tril(ones(n));
    A = zeros(n,n);
    for i = 1:q,
     j = 3*i;
     A(j-2,j-2) = 1 - 1/t^i;
     A(j-1:j,j-1:j) = -A(j-2,j-2) / sqrt(2) * [1 1; -1 1];
    end;
    A = fliplr(tril(ones(n))) * A * tril(ones(n));
    B = 1:n;
    Y = -B' * B;

  else
    error(['Example #%i is not available in Group #%i !',nr(2),nr(1)]);
  end;

else
  error(['Group #%i is not available !',nr(1)]);
end;
