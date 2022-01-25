% HESSOL.F - MEX-function for analysing and solving a system of
%            linear equations with an upper Hessenberg coefficient
%            matrix using SLICOT routines MB02RD, MB02RZ, MB02SD,
%            MB02SZ, MB02TD, and MB02TZ.
%
%   [(X)(,rcnd,LU,ipiv,Hn)] = Hessol(job,H(,ipiv,Hn,mtype,normh,B,tran))
%
%   [rcnd]                = Hessol(-2,LU,ipiv,Hn(,mtype,normh))
%   [X]                   = Hessol(-1,LU,ipiv,mtype,B(,tran))
%   [LU,ipiv(,Hn)]        = Hessol( 0,H(,mtype,normh))
%   [X(,LU,ipiv,Hn)]      = Hessol( 1,H,mtype,normh,B(,tran))
%   [rcnd(,LU,ipiv,Hn)]   = Hessol( 2,H(,mtype,normh))
%   [X,rcnd(,LU,ipiv,Hn)] = Hessol( 3,H,mtype,normh,B(,tran))
%
%  HESSOL performs analysis and solution of a system of linear equations
%     H * X = B,  H' * X = B  or  H**H * X = B, 
%  with a real or complex upper Hessenberg coefficient matrix H.
%  An LU factorization with row pivoting, H = P*L*U, is used, where
%  P is a permutation matrix, L is lower triangular with unit diagonal
%  elements (and one nonzero subdiagonal), and U is upper triangular.
%
%   Description of input parameters:
%   job    - option parameter indicating the task to be performed.
%            =-2 :  condition estimation using LU factorization of H;
%            =-1 :  solution of the system using LU factorization of H;
%            = 0 :  LU factorization of H;
%            = 1 :  LU factorization of H and solution of the system;
%            = 2 :  LU factorization of H and condition estimation;
%            = 3 :  LU factorization of H, condition estimation, and
%                   solution of the system.
%   H      - the n-by-n Hessenberg matrix H or its LU factorization
%            (if job < 0).
%   ipiv   - (optional) if job < 0, the pivot indices defining P in the
%            LU factorization; for 1 <= i <= n, row i of the matrix H
%            was interchanged with row ipiv(i) of H.
%   Hn     - if job = -2, the 1-norm or the infinity-norm of the initial
%            matrix H, as specified by normh.
%   mtype  - (optional) option parameter indicating the type of the
%            matrix H.
%            = 0 :  real matrix;
%            = 1 :  complex matrix.
%            Default:  mtype = 0.
%   normh  - (optional) if job = -2, or job >= 0, specifies whether the
%            1-norm or the infinity-norm reciprocal condition number is
%            required.
%            = 1 :  1-norm;
%            = 2 :  infinity-norm.
%            Default:  normh = 1.
%   B      - (optional) if |job| = 1, or job = 3, the n-by-m matrix B.
%   tran   - (optional) if |job| = 1, or job = 3,  option parameter
%            indicating whether the matrix H or its (conjugate) transpose
%            should be used.
%            = 0 :  use the matrix H;
%            = 1 :  use the matrix H';
%            = 2 :  use the matrix H**H (if mtype = 1).
%            Default:  tran = 0.
%
%   Description of output parameters:
%   X      - if |job| = 1, or job = 3, the n-by-m solution matrix.
%   rcnd   - if |job| = 2, or job = 3, the reciprocal of the condition
%            number of the matrix H, approximating (in the norm defined
%            by normh) 1/(norm(H) * norm(inv(H))).
%   LU     - if job >= 0, the n-by-n matrix containing the factors
%            L and U from the factorization H = P*L*U.
%   ipiv   - (optional) if job >= 0, the pivot indices defining P in the
%            LU factorization; for 1 <= i <= n, row i of the matrix H
%            was interchanged with row ipiv(i) of H.
%   Hn     - if job >= 0, the 1-norm or the infinity-norm of the initial
%            matrix H, as specified by normh.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2002.
%
% Revisions:
%   -
%
