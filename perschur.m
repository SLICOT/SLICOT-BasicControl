% PERSCHUR.F - MEX-function to compute the periodic Schur decomposition
%              using SLICOT routines MB03VD, MB03VY, and MB03WD.
%
%   [Ao(,WR,WI,Z)] = perschur(A(,job,compz,index))
%
%   [Ao(,Z)]       = perschur(A, 0 (,compz,index))
%   [WR,WI(,Z)]    = perschur(A(,1  ,compz,index))
%   [Ao(,WR,WI,Z)] = perschur(A, 2 (,compz,index))
%
%   PERSCHUR computes either the upper Hessenberg form or the Schur 
%   decomposition and the eigenvalues of a product of n-by-n matrices,
%   A = A_1*A_2*...*A_p, p >= 1, using orthogonal transformations, and 
%   without evaluating the product.
%
%   Description of input parameters:
%   A     - the n-by-n-by-p matrices A_1, ..., A_p.
%   job   - (optional) scalar indicating how A should be transformed,
%           as follows:
%           = 0 :  compute the factors H_1, ..., H_p of the upper
%                  Hessenberg form H of A, i.e., H = H_1*H_2*...*H_p,
%                  where H_1 is upper Hessenberg and H_i are upper
%                  triangular, i = 2, ..., p; specifically, the
%                  transformations are
%                     Z_i'*A_i*Z_(i+1) = H_i, i = 1, ..., p-1, and
%                     Z_p'*A_p*Z_1     = H_p;
%           = 1 :  compute the eigenvalues only (default);
%           = 2 :  compute the factors T_1, ..., T_p of the full Schur
%                  form, T = T_1*T_2*...*T_p, where T_1 is in real Schur
%                  form and T_i are upper triangular, i = 2, ..., p;
%                  specifically, the transformations are
%                     Z_i'*A_i*Z_(i+1) = T_i, i = 1, ..., p-1, and
%                     Z_p'*A_p*Z_1     = T_p.
%   compz - (optional) scalar indicating whether or not to compute
%           the orthogonal matrices Z_1, ..., Z_p, as follows:
%           = 0 :  the matrices Z_1,..., Z_p are not computed (default);
%           = 1 :  the matrices Z_1,..., Z_p are computed, i = 1,..., p.
%   index - (optional) vector of length at most 4 containing
%           [ILO; IHI; ILOZ; IHIZ]. If index is specified, it is assumed
%           that all matrices A_j, j = 2, ..., p, are already upper
%           triangular in rows and columns 1:ILO-1 and IHI+1:n, and A_1
%           is upper Hessenberg in rows and columns 1:ILO-1 and IHI+1:n,
%           with A_1(ILO,ILO-1) = 0 (unless ILO = 1), and
%           A_1(IHI+1,IHI) = 0 (unless IHI = n). 
%           Otherwise, ILO = 1, IHI = n.
%           1 <= ILO <= max(1,n); min(ILO,n) <= IHI <= n.
%           If compz = 1, ILOZ and IHIZ specify the rows of Z to which
%           the transformations must be applied.
%           1 <= ILOZ <= ILO; IHI <= IHIZ <= n.
%           Default: index = [1; n; 1; n].
%
%   Description of output parameters:
%   Ao    - if job = 0, the computed n-by-n-by-p matrix H, containing
%           the factors H_1, ..., H_p of the periodic Hessenberg form;
%           if job = 2, the computed n-by-n-by-p matrix T, containing
%           the factors T_1, ..., T_p of the periodic Schur form;
%           if job = 1, Ao is not returned, since no useful information
%           is available.
%   WR    - if job >= 1, the n-vector of real parts of the computed
%           eigenvalues.
%   WI    - if job >= 1, the n-vector of imaginary parts of the computed
%           eigenvalues. Complex conjugate pairs appear consecutively.
%           If index is specified, and ILO > 1 and/or IHI < n, then only
%           the locations ILO to IHI are set in WR and WI.
%   Z     - (optional) if compz = 1, the n-by-n-by-p matrix Z,
%           containing the factors Z_i, i = 1 : p, which reduce the 
%           matrix product to the upper Hessenberg form, or periodic
%           Schur form, according to job. 
%           The transformations are applied only to the submatrices
%           Z_i(ILOZ:IHIZ,ILO:IHI), i = 1, ..., p.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, June 2002.
%
% Revisions:
%   -
%
