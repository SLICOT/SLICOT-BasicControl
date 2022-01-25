% HAMILEIG.F - MEX-function to compute the eigenvalues of a
%              Hamiltonian matrix using the square-reduced approach
%              implemented in SLICOT routines MB03SD and MB04ZD.
%
%   [Ao,QGo(,U)]       = Hamileig(A,QG,job(,compu,S)),        job = -1;
%   [WR,WI(,Ao,QGo,U)] = Hamileig(A,QG(,job,jobscl,compu,S)), job =  0,
%                                                          or job =  1.
%
%   [Ao,QGo(,U)]       = Hamileig(A,QG,-1(,compu,S))
%   [WR,WI(,Ao,QGo,U)] = Hamileig(A,QG(,0, jobscl,compu,S))
%   [WR,WI]            = Hamileig(A,QG, 1(,jobscl))
%
%   HAMILEIG transforms a Hamiltonian matrix
% 
%             ( A   G  )
%         H = (      T )                                           (1)
%             ( Q  -A  )
% 
%   (with G and Q symmetric matrices) into a square-reduced Hamiltonian
%   matrix
% 
%              ( A'   G'  )
%         H' = (        T ),                                       (2)
%              ( Q'  -A'  )
%                                                               T
%   by an orthogonal symplectic similarity transformation H' = U H U,
%   where
%
%               (  U1   U2 )
%           U = (          ),                                      (3)
%               ( -U2   U1 )
%
%   and computes the eigenvalues of H.  Therefore, H' is such that
% 
%         2    ( A''   G'' )
%       H'  =  (         T ),                                      (4)
%              ( 0    A''  )
%
%   with A'' upper Hessenberg and G'' skew symmetric.  The square roots
%   of the eigenvalues of A'' = A'*A' + G'*Q' are the eigenvalues of H.
%
% Description of input parameters: 
%   A      - the n-by-n matrix A.
%   QG     - an  n-by-(n+1) matrix containing the triangles of the 
%            symmetric matrices Q and G, as follows:
%            the leading n-by-n lower triangular part contains the lower
%            triangle of the matrix Q, and the n-by-n upper triangular
%            part of the submatrix in the columns 2 to n+1 contains the
%            upper triangle of the matrix G of H in (1).
%            So, if i >= j, then Q(i,j) = Q(j,i) is stored in QG(i,j)
%            and G(i,j) = G(j,i) is stored in QG(j,i+1).
%            QG is an empty matrix if n = 0.
%   job    - (optional) scalar indicating the computation to be
%            performed, as follows:
%            =-1 :  compute the square-reduced matrix H' in (2);
%            = 0 :  compute the eigenvalues of H (default);
%            = 1 :  compute the eigenvalues of H, assuming that the
%                   given Hamiltonian matrix is already in the reduced
%                   form (2).
%   jobscl - (optional) if job >= 0, scalar specifying whether or not 
%            balancing operations should be performed when computing the
%            eigenvalues of A'', as follows:
%            = 0 :  do not use balancing;
%            = 1 :  do scaling in order to equilibrate the rows
%                   and columns of A'' (default).
%            If job = -1, jobscl is not used.
%   compu  - (optional) if job <= 0, scalar indicating whether the
%            orthogonal symplectic similarity transformation matrix U is
%            returned or accumulated into an orthogonal symplectic
%            matrix, or if the transformation matrix is not required,
%            as follows:
%            = 0 :  U is not required (default);
%            = 1 :  on entry, U need not be set;
%                   on exit, U contains the orthogonal symplectic
%                   matrix U;
%            = 2 :  the orthogonal symplectic similarity transformations
%                   are accumulated into U;
%                   on input, U must contain an orthogonal symplectic
%                   matrix S;
%                   on exit, U contains S*U.
%            If job = 1, compu is not used.
%   S      - (optional) if job <= 0 and compu = 2, an n-by-2*n matrix 
%            containing the first n rows of the given orthogonal
%            symplectic matrix S.
%
% Description of output parameters:
%   WR,WI  - if job >= 0, the n-vectors of real parts and imaginary
%            parts, respectively, of the n computed eigenvalues of H'
%            with non-negative real part. The remaining n eigenvalues
%            are the negatives of these eigenvalues.
%            Eigenvalues are stored in WR and WI in decreasing order of
%            magnitude of the real parts, i.e., WR(I) >= WR(I+1).
%            (In particular, an eigenvalue closest to the imaginary
%             axis is WR(N)+WI(N)i.)
%            In addition, eigenvalues with zero real part are sorted in
%            decreasing order of magnitude of imaginary parts.  Note
%            that non-real eigenvalues with non-zero real part appear
%            in complex conjugate pairs, but eigenvalues with zero real
%            part do not, in general, appear in complex conjugate
%            pairs. 
%   Ao     - if job <= 0, the computed n-by-n submatrix A' of H' in (2).
%   QGo    - if job <= 0, the computed n-by-(n+1) matrix containing the 
%            triangles of the symmetric matrices Q' and G' of H' in (2),
%            stored in the same way as the initial matrices Q and G.
%   U      - if job <= 0 and compu > 0, an n-by-2*n matrix containing
%            the first n rows of the computed orthogonal symplectic 
%            matrix U.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2002.
%
% Revisions:
%   -
%
