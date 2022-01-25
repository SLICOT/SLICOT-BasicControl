% CFSYS.F    - MEX-function for constructing the state-space
%              representation for the system G = (A,B,C,D) from the
%              factors Q and R of its left or right coprime
%              factorization, where G, Q and R are the corresponding
%              transfer-function matrices, using SLICOT routines SB08GD
%              and SB08HD.
%
%   [Ao,Bo,Co,Do(,rcnd)] = cfsys(task,A,B,C,D(,BR)(,CR),DR)
%   [Ao,Bo,Co,Do(,rcnd)] = cfsys(1,A,B,C,D,BR,DR)
%   [Ao,Bo,Co,Do(,rcnd)] = cfsys(2,A,B,C,D,CR,DR)
%
%   CFSYS constructs the state-space representation for the system
%   G = (A,B,C,D) from the factors Q = (AQR,BQ,CQR,DQ) and
%   R = (AQR,BR,CQR,DR) of its left coprime factorization,
%                   -1
%              G = R  * Q,
%
%   or from the factors Q = (AQR,BQR,CQ,DQ) and R = (AQR,BQR,CR,DR)
%   of its right coprime factorization,
%                       -1
%              G = Q * R  ,
%
%   where G, Q and R are the corresponding transfer-function matrices.
%
%   Description of input parameters:
%   task   - integer specifying the computations to be performed.
%            = 1 :  compute the state-space representation from the
%                   factors Q and R of its left coprime factorization;
%            = 2 :  compute the state-space representation from the
%                   factors Q and R of its right coprime factorization.
%   A      - the n-by-n state dynamics matrix AQR of the systems
%            Q and R.
%   B      - if task = 1, the n-by-m input/state matrix BQ of the
%            system Q;
%          - if task = 2, the n-by-m input/state matrix BQR of the
%            systems Q and R.
%   C      - if task = 1, the p-by-n state/output matrix CQR of the
%            systems Q and R;
%          - if task = 2, the p-by-n state/output matrix CQ of the
%            system Q.
%   D      - the p-by-m input/output matrix DQ of the system Q.
%   BR     - if task = 1, the n-by-p input/state matrix BR of the
%            system R.
%   CR     - if task = 2, the m-by-n state/output matrix CR of the
%            system R.
%   DR     - if task = 1, the p-by-p input/output matrix DR of the
%            system R.
%          - if task = 2, the m-by-m input/output matrix DR of the
%            system R.
%
%   Description of output parameters:
%   Ao     - the n-by-n state dynamics matrix of the system G.
%   Bo     - the n-by-m input/state matrix of the system G.
%   Co     - the p-by-n state/output matrix of the system G.
%   Do     - the p-by-m input/output matrix of the system G.
%   rcnd   - (optional) an estimate of the reciprocal condition number
%            of the matrix DR.
%
%   Method:
%   If task = 1, the matrices of the state-space representation
%   G = (A,B,C,D) are computed by using the formulas:
%
%                      -1              -1
%     A = AQR - BR * DR  * CQR,  C = DR  * CQR,
%                      -1              -1
%     B = BQ  - BR * DR  * DQ,   D = DR  * DQ.
%
%   If task = 2, the matrices of the state-space representation
%   G = (A,B,C,D) are computed by using the formulas:
%
%                       -1                   -1
%     A = AQR - BQR * DR  * CR,  B = BQR * DR  ,
%                      -1                   -1
%     C = CQ  - DQ * DR  * CR,   D = DQ * DR  .
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2003.
%
% Revisions:
%   -
%
