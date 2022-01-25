% MUCOMP.F - Gateway function for computing the structured
%            singular value using the SLICOT routine AB13MD.
%
% Matlab call:  
%   [BOUND, D, G] = mucomp(Z, K, T)
%      [BOUND, D] = mucomp(Z, K, T)
%           BOUND = mucomp(Z, K, T)
%
% Purpose:
%   To compute an upper bound on the structured singular value for a 
%   given square complex matrix and given block structure of the
%   uncertainty.
%
% Input parameters:
%   Z      - the complex n-by-n matrix for which the structured 
%            singular value is to be computed.
%   K      - the vector of length m containing the block structure
%            of the uncertainty; K(i) is the size of block i, i = 1:m. 
%   T      - the vector of length m indicating the type of each block.
%            T(i) = 1 if the corresponding block is real;
%            T(i) = 2 if the corresponding block is complex.
%
% Output parameters:
%   BOUND  - the upper bound on the structured singular value.
%   D, G   - vectors of length n containing the diagonal entries
%            of the diagonal matrices D and G, respectively, such that
%            the matrix Z'*D^2*Z + sqrt(-1)*(G*Z-Z'*G) - bound^2*D^2
%            is negative semidefinite.
%
% Comments:
%   Currently, there are two limitations:
%   1. The size of a real block should be 1 (K(i) = 1 if T(i) = 1).
%   2. The sum of block sizes must be equal to n.
% 
% References
%   [1] Fan, M.K.H., Tits, A.L., and Doyle, J.C.
%       Robustness in the presence of mixed parametric uncertainty and
%       unmodeled dynamics.                                           
%       IEEE Trans. Automatic Control, vol. AC-36, 1991, pp. 25-38.   

% RELEASE 2.0 of SLICOT Robust Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   P.Hr. Petkov, TU Sofia, Bulgaria, Oct. 2000.
%
% Revisions:
%   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2005.
