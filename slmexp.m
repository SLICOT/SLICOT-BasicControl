% SLMEXP.F   - MEX-function to compute the exponential of a matrix
%              and, optionally, its integral, using SLICOT routines
%              MB05MD, MB05ND, and MB05OD.
%
%   [F(,V,Y,Wr,Wi)(,mdig,idig)(,H)] = SLMEXP(meth,A(,delta)(,scale)
%                                                  (,ordPad)(,tol))
%
%   [F(,V,Y,Wr,Wi)] = SLMEXP(0,A(,delta,scale))
%   [F(,V,Y,sw,sw)] = SLMEXP(0,A(,delta,scale))
%   [F(,mdig,idig)] = SLMEXP(1,A(,delta,scale,ordPad))
%   [F(,H)]         = SLMEXP(2,A(,delta)(,tol))
%
%   SLMEXP computes F(delta) = exp(A*delta), where A is a real n-by-n
%   matrix and delta is a scalar value. Optionally (for meth = 2), the
%   matrix integral H(delta), defined by
%
%      H(delta) = Int[F(s) ds] from s = 0 to s = delta,
%
%   is returned. Other useful results can also be obtained if required.
%   When meth = 0, matrix A is assumed to be non-defective.
%
%   Description of input parameters:
%   meth   - scalar indicating the method to be used for computing the
%            matrix exponential, as follows:
%            = 0 :  use an eigenvalue/eigenvector decomposition
%                   technique;
%            = 1 :  use a diagonal Pade approximant with scaling and
%                   squaring;
%            = 2 :  use a Pade approximation.
%   A      - the n-by-n matrix A. If delta = 0, the given A is not used.
%   delta  - (optional) the scalar value delta of the problem.
%            Default:  delta = 1.
%   scale  - (optional) if meth = 0 or 1, scalar indicating whether or
%            not the matrix should be diagonally scaled, as follows:
%            = 0 :  do not scale the matrix;
%            = 1 :  diagonally scale the matrix, i.e., replace A by
%                   D*A*D**(-1), where D is a diagonal matrix chosen to
%                   make the rows and columns of A more equal in norm.
%            Default:  scale = 1.
%   ordPad - (optional) if meth = 1, scalar specifying the order of the
%            diagonal Pade approximant. 1 <= ordPad <= 15. In the
%            absence of further information, ordPad should be set to 9.
%            Default:  ordPad = 9.
%   tol    - (optional) if meth = 2, real scalar indicating the
%            tolerance to be used in determining the order of the Pade
%            approximation to H(t), where t is a scale factor determined
%            internally. A reasonable value for tol may be sqrt(EPS),
%            where EPS is the machine precision (see LAPACK Library
%            routine DLAMCH).
%            Default:  tol = sqrt(EPS).
%
%   Description of output parameters:
%   F      - the n-by-n solution matrix exp(A*delta).
%   V      - (optional) if meth = 0, the n-by-n eigenvector matrix
%            for A. If the k-th eigenvalue is real, the k-th column of
%            the eigenvector matrix holds the eigenvector corresponding
%            to the k-th eigenvalue.
%            Otherwise, the k-th and (k+1)-th eigenvalues form a
%            complex conjugate pair and the k-th and (k+1)-th columns
%            of the eigenvector matrix hold the real and imaginary
%            parts of the eigenvectors corresponding to these
%            eigenvalues as follows.
%            If p and q denote the k-th and (k+1)-th columns of the
%            eigenvector matrix, respectively, then the eigenvector
%            corresponding to the complex eigenvalue with positive
%            (negative) imaginary value is given by
%                                      2
%            p + q*j (p - q*j), where j  = -1.
%   Y      - (optional) if meth = 0, the n-by-n real matrix satisfying
%            the equation V*Y = F. If all eigenvalues are real, then Y
%            is the matrix product exp(Lambda*delta) times the inverse
%            of the (right) eigenvector matrix V of A, where Lambda is
%            the diagonal matrix of eigenvalues.
%   Wr,Wi  - (optional) if meth = 0, n-vectors of real and imaginary
%            parts, respectively, of the eigenvalues of the matrix A.
%            The eigenvalues are unordered except that complex conjugate
%            pairs of values appear consecutively, with the eigenvalue
%            having positive imaginary part first.
%   sw     - (optional) if meth = 0, but the matrix A is found to be
%            defective, sw is set to 1, and the calculations are redone
%            automatically for meth = 1 with ordPad = 9. Then, the
%            parameters V and Y (if specified) will contain mdig and
%            idig, respectively. Switching the methods also occurs when
%            sw is not present (the output list has less than four
%            parameters).
%   mdig   - (optional) if meth = 1, the minimal number of accurate
%            digits in the 1-norm of exp(A*delta).
%   idig   - (optional) if meth = 1, the number of accurate digits in
%            the 1-norm of exp(A*delta) at 95% confidence level.
%   H      - (optional) if meth = 2, the n-by-n matrix containing an
%            approximation to H(delta).
%
%   Method
%   If meth = 0, an eigenvalue/eigenvector decomposition technique is
%   used, based on a modification of LAPACK Library routine DGEEV for
%   obtaining the right eigenvector matrix. A condition estimate is
%   then employed to determine if the matrix A is near defective and
%   hence the exponential solution is inaccurate. If A is defective, a
%   warning is returned, and meth = 1 is then used. If the output
%   parameter list in that case had 4 items, the parameters are then
%   F, mdig, idig, and sw (set to 1); for 5 items, sw is duplicated.
%
%   If meth = 1, the exponential of the matrix A is evaluated from a
%   diagonal Pade approximant. The algorithm exploits the identity
%
%       (exp[(2**-m)*A]) ** (2**m) = exp(A),
%
%   where m is an integer determined internally, to improve the accuracy
%   for matrices with large norms.
%
%   If meth = 2, a Pade approximation to H(t) for some small value of t
%   (where 0 < t <= delta) is used, and then F(t) is calculated from
%   H(t). Finally, the results are re-scaled to give F(delta) and
%   H(delta).
%
%   Comments
%   The general recommended method for computing the matrix exponential
%   is meth = 1. Setting meth = 0 could be more efficient for
%   non-defective matrices.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003.
%
% Revisions:
%   -
%
