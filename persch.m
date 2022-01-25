function [Ao,W,Z] = persch(A,job,compz,index)
%PERSCH  computes either the upper Hessenberg form or the Schur 
%        decomposition and the eigenvalues of a product of matrices.
% 
%        W = PERSCH(A)  computes the eigenvalues of a product of n-by-n
%        matrices, A = A_1*A_2*...*A_p, p >= 1, without evaluating 
%        the product.
%
%        [AO,Z]   = PERSCH(A,JOB,COMPZ,INDEX)  with JOB = 0, or
%        [AO,W,Z] = PERSCH(A,JOB,COMPZ,INDEX)  also computes either the upper 
%        Hessenberg form or the Schur form and/or the transformation matrices.
%
%        JOB is a scalar indicating how A should be transformed, as follows:
%        JOB = 0 :  compute the factors H_1, ..., H_p of the upper Hessenberg
%                   form H of A, i.e., H = H_1*H_2*...*H_p, where H_1 is upper
%                   Hessenberg and H_i are upper triangular, i = 2, ..., p;
%                   specifically, the transformations are
%                      Z_i'*A_i*Z_(i+1) = H_i, i = 1, ..., p-1, and
%                      Z_p'*A_p*Z_1     = H_p;
%        JOB = 1 :  compute the eigenvalues only (default);
%        JOB = 2 :  compute the factors T_1, ..., T_p of the full Schur form, 
%                   T = T_1*T_2*...*T_p, where T_1 is in real Schur form and
%                   T_i are upper triangular, i = 2, ..., p; specifically,
%                      Z_i'*A_i*Z_(i+1) = T_i, i = 1, ..., p-1, and
%                      Z_p'*A_p*Z_1     = T_p.
%
%        COMPZ is a scalar indicating whether or not the orthogonal matrices
%        Z_1, ..., Z_p should be computed, as follows:
%        COMPZ = 0 :  the matrices Z_1,..., Z_p are not computed (default);
%        COMPZ = 1 :  the matrices Z_1,..., Z_p are computed, i = 1,..., p.
%
%        INDEX is a vector of length at most 4 containing [ILO; IHI; ILOZ; IHIZ].
%        If INDEX is specified, it is assumed that all matrices A_j, j = 2, ..., p,
%        are already upper triangular in rows and columns 1:ILO-1 and IHI+1:n,
%        and A_1 is upper Hessenberg in rows and columns 1:ILO-1 and IHI+1:n,
%        with A_1(ILO,ILO-1) = 0 (unless ILO = 1), and A_1(IHI+1,IHI) = 0
%        (unless IHI = n).  1 <= ILO <= max(1,n); min(ILO,n) <= IHI <= n.
%        If COMPZ = 1, ILOZ and IHIZ specify the rows of Z to which the
%        transformations must be applied.  1 <= ILOZ <= ILO; IHI <= IHIZ <= n.
%        Default: INDEX = [1; n; 1; n].
%
%        If JOB = 0, AO is the computed n-by-n-by-p matrix H, containing
%        the factors H_1, ..., H_p of the periodic Hessenberg form.
%        If JOB = 2, AO is the computed n-by-n-by-p matrix T, containing
%        the factors T_1, ..., T_p of the periodic Schur form.
%        If JOB = 1, AO is not returned, and the call of the function is
%        [W(,Z)] = PERSCH(A,1(,COMPZ,INDEX))  with possible optional arguments
%        appearing in the brackets.
%
%        W is an n-vector containing the computed eigenvalues. Complex conjugate
%        pairs appear consecutively. If INDEX is specified, and ILO > 1 and/or
%        IHI < n, then only the locations ILO to IHI are set in W.
%
%        If COMPZ = 1, Z is the computed n-by-n-by-p matrix Z, containing the
%        factors Z_i, i = 1, ..., p, which reduce the matrix product to either
%        the upper Hessenberg form, or periodic Schur form, according to JOB.
%        The transformations are applied only to the submatrices
%        Z_i(ILOZ:IHIZ,ILO:IHI), i = 1, ..., p.
%
%        [AO(,Z)] = PERSCH(A,0(,COMPZ,INDEX))  returns the periodic Hessenberg
%        form AO, and, optionally, the corresponding transformation matrix.
%
%        Comments
%        A, AO, and Z (if specified) are stored as tri-dimensional arrays, with
%        M(1:n,1:n,i) containing the matrix M_i, i = 1 : p.
%
%        See also PERSCHUR
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, June 2002.
%        Revised: V. Sima, July 2002, March 2009.
%

ni   = nargin;
nout = nargout;
%
if min( ni, nout ) < 1,
    error(['Usage: W  = persch(A)',  sprintf('\n'),...
           '       AO = persch(A,0)',sprintf('\n'),...
           '       AO = persch(A,2)'])
end
[n, p] = size(A);
if n == 0,
    p = 1;
else
    if ( fix( p/n ) )*n == p,
        p = p/n;
    else
        error(['The length of A must be a multiple of ',num2str( n*n )])
    end
end
if ni == 1,
    job   = 1;
    compz = 0;
    index = [1; n; 1; n];
elseif ni == 2,
    compz = 0;
    index = [1; n; 1; n];
elseif ni == 3,
    index = [1; n; 1; n];
end
%
if nout == 1,
    if job == 1,
        [Wr,Wi] = perschur(A,job,compz,index);  Ao = Wr + 1i*Wi;
    else
        Z = perschur(A,job,compz,index);  Ao = zeros(n,n,p);
        for ii = 1 : p, 
            Ao(:,:,ii) = Z(:,(ii-1)*n+1:ii*n);  
        end
    end
elseif nout == 2,
    if job == 0,
        [Z,Ww] = perschur(A,job,compz,index);  Ao = zeros(n,n,p);  W = zeros(n,n,p);
        for ii = 1 : p, 
            Ao(:,:,ii) =  Z(:,(ii-1)*n+1:ii*n);  
            W(:,:,ii)  = Ww(:,(ii-1)*n+1:ii*n);  
        end
    elseif job == 1,
        [Wr,Wi,Z] = perschur(A,job,compz,index);  Ao = Wr + 1i*Wi;  W = zeros(n,n,p);
        for ii = 1 : p, 
            W(:,:,ii) = Z(:,(ii-1)*n+1:ii*n);  
        end
    elseif job == 2,
        [Z,Wr,Wi] = perschur(A,job,compz,index);  Ao = zeros(n,n,p);  W = Wr + 1i*Wi;
        for ii = 1 : p, 
            Ao(:,:,ii) = Z(:,(ii-1)*n+1:ii*n);  
        end
    end
elseif nout == 3,
    [Aw,Wr,Wi,Zw] = perschur(A,job,compz,index);  Ao = zeros(n,n,p);  Z = zeros(n,n,p);
    for ii = 1 : p, 
        Ao(:,:,ii) = Aw(:,(ii-1)*n+1:ii*n);  
        Z(:,:,ii)  = Zw(:,(ii-1)*n+1:ii*n);  
    end
    W = Wr + 1i*Wi;
end
%
% end persch
