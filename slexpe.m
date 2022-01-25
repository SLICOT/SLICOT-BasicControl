function [F,V,Y,W] = slexpe(A,delta,scale)
% SLEXPE Computes the exponential of a matrix using an
%        eigenvalue/eigenvector decomposition technique, or
%        a diagonal Pade approximant with scaling and squaring,
%        if the matrix appears to be defective.
%
%        F = SLEXPE(A)  computes the matrix exponential, F = exp(A).
%
%        F = SLEXPE(A,DELTA)  computes F = exp(A*DELTA), where DELTA
%        is a real scalar.
%
%        [F,V,Y,W] = SLEXPE(A,DELTA,SCALE)  has additional input and
%        output arguments.
%
%        SCALE is an integer scalar indicating whether or not the matrix
%        should be diagonally scaled:
%        SCALE = 0 :  do not scale the matrix;
%        SCALE = 1 :  diagonally scale the matrix, i.e., replace A by
%                     D*A*D**(-1), where D is a diagonal matrix chosen to
%                     make the rows and columns of A more equal in norm.
%        Default:  SCALE = 1.
%
%        V is the eigenvector matrix for A, but if prod(size(V)) is 1,
%        and W equals sqrt(-1), then V stores MDIG (see below).
%
%        Y is a real matrix satisfying the equation Vr*Y = F, where Vr is
%        the real, compact representation of V. If the k-th eigenvalue is 
%        real, the k-th column of Vr holds the eigenvector corresponding
%        to that eigenvalue. Otherwise, the k-th and (k+1)-th columns of
%        Vr hold the real and imaginary parts of the (right) eigenvector 
%        corresponding to the complex eigenvalue with positive imaginary
%        part. If all eigenvalues are real, then Y = exp(W*DELTA)*inv(V).
%        If prod(size(Y)) is 1, and W equals sqrt(-1), then Y stores IDIG
%        (see below).
%
%        W is a column vector of eigenvalues. The eigenvalues are
%        unordered except that complex conjugate pairs of values
%        appear consecutively, with the eigenvalue having positive
%        imaginary part first.
%
%        [F,MDIG,IDIG,W] = SLEXPE(A,DELTA,SCALE)  with W = sqrt(-1)
%        means that the eigenvalue/eigenvector decomposition technique
%        failed, and the Pade approximation was used instead.
%
%        MDIG is the minimal number of accurate digits in the 1-norm
%        of exp(A*DELTA).
%
%        IDIG is the number of accurate digits in the 1-norm of
%        exp(A*DELTA) at 95% confidence level.
%
%        See also SLMEXP, SLEXPM, SLEXPI
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, April 2003.
%
%        Revisions:
%        March 2009.
%

ni   = nargin;
nout = nargout;
%
if min( ni, nout ) < 1,
    error(['Usage: F               = SLEXPE(A)', sprintf('\n'),...
           '       [F,V,Y,W]       = SLEXPE(A,DELTA,SCALE)'])
end
% 
if ni == 1,
    delta = 1;
    scale = 1;
elseif ni == 2,
    scale = 1;
end
%
[ F, V, Y, Wr, Wi ] = slmexp( 0, A, delta, scale );
if min( length( V ), length( Y ) ) == 1 && Wi ~= 0,  
    W = 1i;  
else  
    W = Wr + 1i*Wi;  
end
%
if nout > 1 && length( V ) > 1,
    ind = find( Wi > 0 );
    for ii = ind,
        V(:,ii)   = V(:,ii) + 1i*V(:,ii+1);
        V(:,ii+1) = conj( V(:,ii) );
    end
end
%
% end slexpe
