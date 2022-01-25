function [E,A,B,C,D] = dtdsx(nr,parin)
%DTDSX
%
% Usage:  [E,A,B,C,D] = dtdsx(nr,parin)
%         [E,A,B,C,D] = dtdsx(nr)
%
% Main routine of the benchmark library DTDSX (Version 1.0) described 
% in [1]. It generates benchmark examples for time-invariant, 
% discrete-time, dynamical systems
%
%       E x_k+1 = A x_k + B u_k
%                                                                   (1)
%           y_k = C x_k + D u_k
%
% E, A are real n-by-n matrices, B is n-by-m, C is p-by-n, and 
% D is p-by-m. 
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
%            Ex. 2.1 : parin(1:3) = [tau delta K]
%            Ex. 3.1 : parin(1)   = n
%            parin is optional; default values as defined in [1] are
%            used as example parameters if 'parin' is omitted. Note that
%            parin is not referenced if nr(1) = 1.
%
% Output:
%  - E, A, B, C, D :  matrices of the dynamical system (1).
%
% References: 
%
% [1] D. Kressner, V. Mehrmann, and T. Penzl.
%     DTDSX - a Collection of Benchmark Examples for State-Space 
%     Realizations of Discrete-Time Dynamical Systems.
%     SLICOT working note 1998-10. 1998.
%
%     For questions concerning the collection or for the submission of
%     test examples, please contact Volker Mehrmann
%     (Email: volker.mehrmann@mathematik.tu-chemnitz.de).

%  RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%  Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%  D. Kressner, V. Mehrmann, and T. Penzl (TU Chemnitz).
%  Dec 1, 1998. 
%
%  Revisions:
%  V. Sima, Research Institute for Informatics, Bucharest, Oct. 2004,
%           Mar. 2009.

if nargin < 1,
  error('Not enough input arguments.');
end;

if length(nr) < 2,
  error('Please use the nr = [group, example] notation.');
end;

	 
if nr(1) == 1,
  if nr(2) == 1,
    %  Example 1.1:  Laub 1979, Ex. 2: uncontrollable-unobservable data
    E = eye(2);
    A = [4, 3; -4.5, -3.5];
    B = [1; -1];  
    C = [3, 2];                    
    D = zeros(size(C,1),size(B,2));

  elseif nr(2) == 2,
    %  Example 1.2:  Laub 1979, Ex. 3
    E = eye(2);
    A = [0.9512, 0; 0, 0.9048];
    B = [4.877, 4.877; -1.1895, 3.569];
    C = eye(2);
    D = zeros(size(C,1),size(B,2));

  elseif nr(2) == 3,
    %  Example 1.3:  Van Dooren 1981, Ex. II
    E = eye(2);
    A = [2, -1; 1, 0];
    B = [1; 0];
    C = [0 1];
    D = zeros(size(C,1),size(B,2));

  elseif nr(2) == 4,
    %  Example 1.4:  Ionescu/Weiss 1992
    E = eye(2);
    A = [0, 1; 0, -1];
    B = [1, 0; 2, 1];
    C = eye(2);
    D = zeros(size(C,1),size(B,2));

  elseif nr(2) == 5,
    %  Example 1.5:  Jonckheere 1981
    E = eye(2);
    A = [0, 1; 0, 0];
    B = [0; 1];
    C = eye(2);
    D = zeros(size(C,1),size(B,2));

  elseif nr(2) == 6,
    %  Example 1.6:  Ackerson/Fu 1970: satellite control problem
    E = eye(4);
    A = 0.998*eye(4);
    A(1,2) = 0.067;   A(2,1) = -A(1,2);
    A(3,4) = 0.153;   A(4,3) = -A(3,4);
    B = [0.0033, 0.02; 0.1, -0.0007; 0.04, 0.0073; -0.0028, 0.1];
    C = eye(4);
    D = zeros(size(C,1),size(B,2));

  elseif nr(2) == 7,
    %  Example 1.7:  Litkouhi 1983: system with slow and fast modes
    E = eye(4);
    A = [ 0.98475,  -0.079903,  0.0009054, -0.0010765;... 
          0.041588,  0.99899,  -0.035855,   0.012684;...
         -0.54662,   0.044916, -0.32991,    0.19318;...
          2.6624,   -0.10045,  -0.92455,   -0.26325];
    B = [0.0037112, 0.0007361; -0.087051, 9.3411e-6; -1.19844, -4.1378e-4;...
         -3.1927, 9.2535e-4];
    C = eye(4);
    D = zeros(size(C,1),size(B,2));

  elseif nr(2) == 8,
    %  Example 1.8:  Lu/Lin 1993, Ex. 4.3
    E = eye(4);
    B  = -triu(ones(4)) + diag(2*ones(4,1));
    C  = triu(ones(4)) + diag(ones(2,1),2) + diag(3,3);
    A  = B*[0.4, 0, 0, 0; 1, 0.6, 0, 0; 0, 1, 0.8, 0; 0, 0, 0, -0.999982]*C;
    D = zeros(size(C,1),size(B,2));

  elseif nr(2) == 9,
    %  Example 1.9:  Gajic/Shen 1993, Section 2.7.4: chemical plant
    E = eye(5);
    A = 0.01*[ 95.407,   1.9643,  0.3597,  0.0673,  0.019;... 
               40.849,  41.317,  16.084,   4.4679,  1.1971;...
               12.217,  26.326,  36.149,  15.93,   12.383;...
                4.1118, 12.858,  27.209,  21.442,  40.976;...
                0.1305,  0.5808,  1.875,   3.6162, 94.28];
    B = 0.01*[ 0.0434,  2.6606,  3.753,  3.6076,  0.4617;...
              -0.0122, -1.0453, -5.51,  -6.6,    -0.9148]';
    C = eye(5);
    D = zeros(size(C,1),size(B,2));

  elseif nr(2) == 10,
    %  Example 1.10:  Davison/Wang 1974
    E = eye(6);
    A = kron(eye(2), diag([1,1], 1));
    B = [0, 0, 1, 0, 0, 0;  0, 0, 0, 0, 0, 1]'; 
    C = [1, 1, 0, 0, 0, 0;  0, 0, 0, 1, -1, 0];
    D = [1, 0;  1, 0];

  elseif nr(2) == 11,
    %  Example 1.11:  Patnaik et al. 1980: tubular ammonia reactor
    E = eye(9);
    A  = [870.1,  135.0   11.59   0.5014 -37.22  0.3484 0.0     4.242  7.249;...
           76.55 897.4   12.72   0.5504 -40.16  0.3743 0.0     4.53   7.499;...
         -127.2  357.5  817.0    1.455 -102.8   0.987  0.0    11.85  18.72;...
         -363.5  633.9   74.91 796.6   -273.5   2.653  0.0    31.72  48.82;...
         -960.0 1645.9 -128.9   -5.597   71.42  7.108  0.0    84.52 125.9;...
         -664.4  112.96 -88.89  -3.854   84.47 13.6    0.0   144.3  101.6;...
         -410.2  693.0  -54.71  -2.371   66.49 12.49   0.1063 99.97  69.67;...
         -179.9  301.7  -23.93  -1.035   60.59 22.16   0.0   213.9   35.54;...
         -345.1  580.4  -45.96  -1.989  105.6  19.86   0.0   219.1  215.2 ];
    A  = 0.001*A;
    B  = [  4.76    0.879   1.482  3.892 10.34   7.203  4.454  1.971  3.773;...
           -0.5701 -4.773 -13.12 -35.13 -92.75 -61.59 -36.83 -15.54 -30.28;...
          -83.68   -2.73    8.876 24.8   66.8   38.34  20.29   6.937 14.69 ]';
    B  = 0.0001*B;
    C  = [1, zeros(1,8); zeros(1,4), 1, zeros(1,4)];
    D = zeros(size(C,1),size(B,2));

  elseif nr(2) == 12,
    %  Example 1.12:  Smith 1969: two-stand cold rolling mill
    E = eye(10);
    A = [zeros(1,9) 0.1120; eye(9) zeros(9,1)];
    B = [ 2.7600  -1.3500  -0.4600; zeros(9,3) ];
    C = [eye(5,1) zeros(5,8) ...
        [0 0.8940 -16.9300 0.0700 0.3980 ]'];
    D = [zeros(1,3);
         -0.2230   1.8500 -0.5420;...
         28.3000 204.0000 68.7000;...
         -5.2100  -0.8430 -0.2850;...
         -0.1010  -6.7500 -0.2460];

  else
    error(['Example #%i is not available in Group #%i !',nr(2),nr(1)]);
  end;


elseif nr(1) == 2,
  if nr(2) == 1,
    %  Example 2.1:  Pappas et al. 1980: process control of paper machine
    if nargin < 2,
      parin = [10^8 1 1];
    elseif length(parin) < 3,
      error('Not enough parameters.');
    end;
    tau   = parin(1);
    delta = parin(2);
    K     = parin(3);
    if tau == 0, error('Parameter tau must not be zero.'); end
    alpha = 1 - delta/tau;
    beta  = K*delta/tau;

    E = eye(4);
    A = diag(ones(3,1),-1);  A(1,1) = alpha;
    B = [beta; 0; 0; 0];
    C = [0, 0, 0, 1];
    D = zeros(size(C,1),size(B,2));

  else
    error(['Example #%i is not available in Group #%i !',nr(2),nr(1)]);
  end;


elseif nr(1) == 3,
  if nr(2) == 1,
    %  Example 3.1:  Pappas et al. 1980, Ex. 3
    if nargin < 2, n = 100;  else n = parin(1); end 
    if n ~= round(n) || n < 2, error('Invalid value of parameter n.'); end

    E = eye(n);
    A = diag(ones(n-1,1),1);
    B = flipud(eye(n,1));
    C = eye(n);
    D = zeros(size(C,1),size(B,2));

  else
    error(['Example #%i is not available in Group #%i !',nr(2),nr(1)]);
  end;


elseif nr(1) == 4,
  error('There are no examples available in Group 4 !');

else
  error(['Group #%i is not available !',nr(1)]);
end;