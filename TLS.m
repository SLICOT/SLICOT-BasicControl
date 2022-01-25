function [ X, ro, rcnd, to ] = TLS( A, B, meth, tol, r, printw, t )
%TLS     Solves the Total Least Squares (TLS) problem using a 
%        singular value decomposition (SVD) approach or a 
%        Partial SVD (PSVD) approach.
%
%        X = TLS( A, B )  solves the TLS problem A*X = B using
%        the PSVD approach and default options.
%        
%        X = TLS( A, B, METH )  solves the TLS problem A*X = B using
%        default options and the approach specified by METH:
%        METH = 1 :  SVD  approach;
%        METH = 2 :  PSVD approach.
%        Default:  METH = 2.
%        
%        [ X, Ro, RCND ]     = TLS( A, B, 1, TOL, R, PRINTW )  or
%        [ X, Ro, RCND, To ] = TLS( A, B, 2, TOL, R, PRINTW, T )  have 
%        additional input and output arguments.
%        
%        TOL is a real scalar defining the tolerance to be used.
%        If METH = 1, the tolerance is used to determine the
%        rank of the TLS approximation [A+DA|B+DB] and to check
%        the multiplicity of the singular values of matrix C = [A|B].
%        Specifically, S(i) and S(j) (i < j) are considered to be
%        equal if SQRT(S(i)**2 - S(j)**2) <= TOL, and the TLS
%        approximation [A+DA|B+DB] has rank R if S(i) > TOL*S(1),
%        for i = 1,2,...,R. The value TOL is also used to check
%        the singularity of an upper triangular matrix F (see TOTALLS).
%        If TOL <= 0, the tolerance is taken as eps, where eps is 
%        the machine precision.
%        If TOL is missing, the tolerance is computed internally 
%        using sdev = eps, where sdev is the estimated standard
%        deviation of the error on each element of the matrix C.
%        (See TOTALLS for job = 3).
%        If METH = 2, TOL defines the multiplicity of singular
%        values by considering all singular values within an interval
%        of length TOL as coinciding. The value TOL is used in
%        checking how many singular values are less than or equal
%        to T (see T below). Also in computing an appropriate upper
%        bound T by a bisection method, TOL is used as a stopping
%        criterion defining the minimum (absolute) subinterval width.
%        The value TOL is also taken as an absolute tolerance for
%        negligible elements in the QR/QL iterations.
%        If TOL is missing, or TOL <= 0, then the tolerance is
%        computed internally.
%
%        R is an integer scalar. If R >= 0, R specifies the rank of
%        the TLS approximation [A+DA|B+DB];  R <= min(size(A)).
%        Otherwise, or if R is missing, the rank is computed internally.
%        Default:  R = -1.
%
%        PRINTW is a switch for printing the warning messages.
%        PRINTW = 1 :  print warning messages;
%               = 0 :  do not print warning messages.
%        Default:  PRINTW = 0.
%
%        T is a real scalar. If METH = 2, and R < 0, then the rank of
%        the TLS approximation [A+DA|B+DB] is computed using the real
%        number T as (min(size(C)) - d), where d is the number of
%        singular values of [A|B] less than or equal to T;  T >= 0.
%        If METH = 2 and R >= 0, T is an initial estimate for
%        computing a lower bound on the R largest singular values
%        of [A|B]. If T < 0 on entry however, then T is computed
%        internally.
%        Default:  T = eps, if R <  0;
%                  T = -1,  if R >= 0.
%
%        X is the solution matrix to the TLS problem specified by
%        A and B.
%
%        Ro is the computed rank of the TLS approximation [A+DA|B+DB].
%        The input value R (if given) may be changed if the R-th
%        and the (R+1)-th singular values of C are considered to be
%        equal, or if the upper triangular matrix F is (numerically)
%        singular.
%
%        RCND is a real scalar containing the reciprocal of the
%        condition number of the matrix F.
%
%        To is an integer scalar. If METH = 2 and R >= 0 on entry, then
%        To contains the computed bound such that precisely Ro singular
%        values of C are greater than To + TOL.
%        If METH = 2, and R < 0, then To = T.
%
%        See also TOTALLS
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Sep-16-2003.
%        Revised Mar-03-2009.
%

ni = nargin;  nout = nargout;
if ni < 2 || nout < 1,  
    error( ['Usage: X = TLS( A, B )', sprintf('\n'), ...
            '       X = TLS( A, B, METH )', sprintf('\n'), ...
            '       [ X, Ro, RCND ]     = TLS( A, B, 1, TOL, R, PRINTW )', sprintf('\n'), ...
            '       [ X, Ro, RCND, To ] = TLS( A, B, 2, TOL, R, PRINTW, T )' ] );
end
%
if ni == 2,
    meth = 2;  tol = 0;  r = -1;  printw = 0;  t = eps;
elseif ni == 3,
    if meth == 1,  tol = eps;  else  tol = -1;  t = eps;  end
    r = -1;  printw = 0;
elseif ni == 4,
    r = -1;  printw = 0;  t = eps;
elseif ni == 5,
    printw = 0;  
    if meth == 2,
        if r < 0,  t = eps;  else  t = -1;  end
    end
elseif ni == 6 && meth == 2,
    if r < 0,  t = eps;  else  t = -1;  end
end
if ( meth == 1 && ni > 6 ) || ( meth == 2 && ni > 7 )
    error( 'Too many input parameters' );
end
%
if meth == 1,
    if ni <= 3,
        job = 3;
        [ X, V, S, ro, rcnd ] = TotalLS( A, B, meth, job, tol, printw );
    elseif ni <= 4 || r < 0,
        job = 1;
        [ X, V, S, ro, rcnd ] = TotalLS( A, B, meth, job, tol, printw );
    else
        job = 0;
        [ X, V, S, ro, rcnd ] = TotalLS( A, B, meth, job, r, tol, printw );
    end
else
    reltol = 0;
    [ X, V, Q, NL, ro, to, rcnd ] = TotalLS( A, B, meth, r, t, tol, reltol, printw );
end
%
% end TLS
