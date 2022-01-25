% CONDIS.F   - MEX-function to perform a transformation on the
%              parameters (A,B,C,D) of a system (equivalent to a
%              bilinear transformation of the corresponding transfer
%              function matrix), using SLICOT routine AB04MD.
%
%   [Ao,Bo,Co,Do] = CONDIS(task,A,B,C,D(,alpha,beta))
%
%   CONDIS performs a transformation on the parameters (A,B,C,D) of a
%   system, which is equivalent to a bilinear transformation of the
%   corresponding transfer function matrix.
%
%   Description of input parameters: 
%   task   - integer specifying the type of the transformation to be
%            performed.
%            task = 1 :  discrete-time   -> continuous-time;
%            task = 2 :  continuous-time -> discrete-time.
%   A      - the n-by-n state dynamics matrix A.
%   B      - the n-by-m input/state matrix B.
%   C      - the p-by-n state/output matrix C.
%   D      - the p-by-m input/output matrix D.
%   alpha, - nonzero parameters specifying the bilinear transformation.
%   beta     Recommended values for stable systems: alpha = 1, beta = 1.
%            Default:  alpha = 1, beta = 1.
%
%   Description of output parameters:
%                                                         _
%   Ao     - the n-by-n transformed state dynamics matrix A.
%                                                       _
%   Bo     - the n-by-m transformed input/state  matrix B.
%                                                       _
%   Co     - the p-by-n transformed state/output matrix C.
%                                                       _
%   Do     - the p-by-m transformed input/output matrix D.
%
%   Method
%   The parameters of the discrete-time system are transformed into
%   the parameters of the continuous-time system (task = 1), or
%   vice-versa (task = 2) by the transformation:
%
%   1.  Discrete -> continuous
%       _                     -1
%       A = beta*(alpha*I + A)  * (A - alpha*I)
%       _                                     -1
%       B = sqrt(2*alpha*beta) * (alpha*I + A)  * B
%       _                                         -1
%       C = sqrt(2*alpha*beta) * C * (alpha*I + A)
%       _                        -1
%       D = D - C * (alpha*I + A)  * B
%
%   which is equivalent to the bilinear transformation
%
%                     z - alpha
%       z -> s = beta --------- 
%                     z + alpha
%
%   of one transfer matrix onto the other.
%
%   2.  Continuous -> discrete
%       _                     -1
%       A = alpha*(beta*I - A)  * (beta*I + A)
%       _                                    -1
%       B = sqrt(2*alpha*beta) * (beta*I - A)  * B
%       _                                        -1
%       C = sqrt(2*alpha*beta) * C * (beta*I - A)
%       _                       -1
%       D = D + C * (beta*I - A)  * B
%
%   which is equivalent to the bilinear transformation
%
%                      beta + s
%       s -> z = alpha -------- 
%                      beta - s
%
%   of one transfer matrix onto the other.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, July 2003.
%
% Revisions:
%   -
%
