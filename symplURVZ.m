% SYMPLURVZ.F - Gateway function for computing the eigenvalues of a
%               complex skew-Hamiltonian/Hamiltonian pencil in factored
%               form, using SLICOT routine MB04AZ.
%               The gateway accepts real or complex input matrices.
%
%   [ALPHA,BETA(,Ao,Do,Co,Bo,Fo)(,Q,U)] = symplURVZ(Z,B,FG
%                                                   (,job,compq,compu))
%   [ALPHA,BETA(,Q,U)] =              symplURVZ(Z,B,FG(,0,compq,compu))
%
%   SYMPLURVZ computes the eigenvalues of a complex N-by-N skew-Hamiltonian/
%   Hamiltonian pencil aS - bH, with
%
%            H  T           (  B  F  )       (  Z11  Z12  )
%     S = J Z  J  Z and H = (      H ), Z =: (            ).       (1)
%                           (  G -B  )       (  Z21  Z22  )
%
%   The structured Schur form of the embedded real skew-Hamiltonian/
%                                                         H  T
%   skew-Hamiltonian pencil, aB_S - bB_T, with B_S = J B_Z  J  B_Z,
%
%           (  Re(Z11)  -Im(Z11)  |  Re(Z12)  -Im(Z12)  )
%           (                     |                     )
%           (  Im(Z11)   Re(Z11)  |  Im(Z12)   Re(Z12)  )
%           (                     |                     )
%     B_Z = (---------------------+---------------------) ,
%           (                     |                     )
%           (  Re(Z21)  -Im(Z21)  |  Re(Z22)  -Im(Z22)  )
%           (                     |                     )
%           (  Im(Z21)   Re(Z21)  |  Im(Z22)   Re(Z22)  )
%                                                                  (2)
%           ( -Im(B)  -Re(B)  | -Im(F)  -Re(F)  )
%           (                 |                 )
%           (  Re(B)  -Im(B)  |  Re(F)  -Im(F)  )
%           (                 |                 )
%     B_T = (-----------------+-----------------) ,  T = i*H,
%           (                 |      T       T  )
%           ( -Im(G)  -Re(G)  | -Im(B )  Re(B ) )
%           (                 |      T       T  )
%           (  Re(G)  -Im(G)  | -Re(B ) -Im(B ) )
%
%   is determined and used to compute the eigenvalues. Optionally, if
%   job = 1, the pencil aB_S - bB_H is transformed by a unitary matrix Q
%   and a unitary symplectic matrix U to the structured Schur form
%                                            H  T
%   aB_Sout - bB_Hout, with B_Sout = J B_Zout  J  B_Zout,
%
%              ( Ao  Do  )              ( Bo  Fo  )
%     B_Zout = (         ) and B_Hout = (       H ),               (3)
%              (  0  Co  )              (  0 -Bo  )
%
%   where Ao and Bo are upper triangular, Co is lower triangular, and Fo
%   is Hermitian. The embedding doubles the multiplicities of the
%   eigenvalues of the pencil aS - bH.
%   Optionally, if compq > 0, the unitary matrix Q is computed.
%   Optionally, if compu > 0, the unitary symplectic matrix U is
%   computed.
%
% Description of input parameters:
%   Z      - the n-by-n matrix Z, n = 2*m.
%   B      - the m-by-m matrix B.
%   FG     - an  m-by-(m+1) matrix containing the triangles of the
%            Hermitian matrices F and G, as follows:
%            the leading m-by-m lower triangular part contains the lower
%            triangle of the matrix G, and the m-by-m upper triangular
%            part of the submatrix in the columns 2 to m+1 contains the
%            upper triangle of the matrix F of H in (1).
%            So, if i >= j, then G(i,j) = G(j,i) is stored in FG(i,j)
%            and F(j,i) = F(i,j) is stored in FG(j,i+1).
%            FG is an empty matrix if m = 0.
%   job    - (optional) scalar indicating the computation to be
%            performed, as follows:
%            = 0 :  compute the eigenvalues only (default);
%            = 1 :  compute the eigenvalues and the matrices of the
%                   transformed pencil in (3).
%   compq  - (optional) scalar indicating whether the unitary
%            transformation matrix Q is returned, or if Q is not
%            required, as follows:
%            = 0 :  Q is not required (default);
%            = 1 :  on exit, Q contains the unitary matrix Q.
%   compu  - (optional) scalar indicating whether the unitary symplectic
%            transformation matrix U is returned, or if U is not
%            required, as follows:
%            = 0 :  U is not required (default);
%            = 1 :  on exit, U contains the relevant part of the unitary 
%                   symplectic matrix U (see below).
%
% Description of output parameters:
%   ALPHA  - the n-vector of the numerators alpha defining the
%            eigenvalues of the pencil aS - bH.
%   BETA     the n-vector of the denominators beta defining the
%            eigenvalues of the pencil aS - bH.
%            Together, the quantities alpha = ALPHA(j) and
%            beta = BETA(j) represent the j-th eigenvalue of the pencil
%            aS - bH, in the form lambda = alpha/beta. Since lambda may
%            overflow, the ratios should not, in general, be computed.
%   Ao     - (optional) if job = 1, the n-by-n matrix Ao in (3).
%            The strictly lower triangular part is not zeroed.
%   Do     - (optional) if job = 1, the n-by-n matrix Do in (3).
%   Co     - (optional) if job = 1, the n-by-n matrix Co in (3).
%   Bo     - (optional) if job = 1, the n-by-n matrix Bo in (3).
%            The strictly lower triangular part is not zeroed.
%   Fo     - (optional) if job = 1, the n-by-n matrix Fo in (3).
%            The strictly lower triangular part is not zeroed.
%   Q      - (optional) if compq > 0, a 2*n-by-2*n matrix containing the
%            unitary matrix Q.
%   U      - (optional) if compu > 0, the n-by-2*n leading rows of the
%            unitary symplectic matrix U.
%            If job = 0 and compq > 0 or compu > 0, Q or U contain the
%            corresponding matrices which transformed B_Z and B_T
%            in (2). Their imaginary part is zero.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, July 2012.
%
% Revisions:
%   -
%
