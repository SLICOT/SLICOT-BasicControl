function W = Hameig( A, QG, scale )
%HAMEIG  computes the eigenvalues of a Hamiltonian matrix given in a
%        compressed form.
% 
%        W = HAMEIG(A,QG,SCALE)  computes the eigenvalues of a Hamiltonian 
%        matrix 
% 
%            ( A   G  )
%        H = (      T ),  with A, Q, and G n-by-n.                      (1)
%            ( Q  -A  )
%        
%        The symmetric matrices Q and G are stored in the n-by-(n+1) array QG, 
%        as follows: the leading n-by-n lower triangular part contains the
%        lower triangle of the matrix Q, and the n-by-n upper triangular part
%        of the submatrix in the columns 2 to n+1 contains the upper triangle
%        of the matrix G of H in (1). So, if i >= j, then Q(i,j) = Q(j,i) 
%        is stored in QG(i,j) and G(i,j) = G(j,i) is stored in QG(j,i+1).
%
%        SCALE is a scalar specifying whether or not balancing operations
%        should be performed when computing the eigenvalues, as follows:
%        SCALE = 0 :  do not use balancing;
%        SCALE = 1 :  do scaling in order to equilibrate the rows and columns.
%        Default: SCALE = 1.
%
%        W is an 2*n-vector containing the eigenvalues of H. 
%
%        Comments
%        Eigenvalues are stored in W as follows: The first n entries of W are
%        the computed eigenvalues of H with non-negative real part. They are
%        stored in decreasing order of magnitude of the real parts, i.e.,
%        real( W(I) ) >= real( W(I+1) ). (In particular, an eigenvalue closest
%        to the imaginary axis is W(N).) In addition, eigenvalues with zero real
%        part are sorted in decreasing order of magnitude of imaginary parts.
%        Note that non-real eigenvalues with non-zero real part appear in
%        complex conjugate pairs, but eigenvalues with zero real part (in the 
%        first n entries of W) do not, in general, appear in complex conjugate
%        pairs. The remaining n eigenvalues are the negatives of the first
%        n eigenvalues. 
%
%        See also HAMILEIG
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, Nov. 2002.
%        Revised: Mar. 2009.
%

ni = nargin;
%
if ni < 2 || nargout < 1,
    error( [ 'Usage: W = Hameig(A,QG)',  sprintf('\n'),...
             '       W = Hameig(A,QG,0)' ] )
end
%
[n, m] = size( A );  [n1, m1] = size( QG );
if ~( m == n && n1 == n && m1 == n + min( n, 1 ) ),
    error( 'Arrays A and/or QG have wrong or incompatible sizes' )
end
%
if ni < 3,  scale = 1;  end
%
job = 0;
[ Wr, Wi ] = Hamileig( A, QG, job, scale );
W = Wr + 1i*Wi;  W(n+1:2*n) = -W;  
for ii = n + 1 : 2*n,
    if ~( Wr(ii-n) == 0 ),  W(ii) = conj( W(ii) );  end
end
W = reshape( W, 2*n, min( n, 1 ) );
%
% end Hameig
