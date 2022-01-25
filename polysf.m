function [ E, res, B ] = polysf( A, discr, form )
%POLYSF  Computes the spectral factorization of a real polynomial,
%        arising from optimality problems.
%
%        E = POLYSF(A)  computes a real polynomial E(s) such that
%            (a)  E(-s) * E(s) = A(-s) * A(s) =: B(s), and
%            (b)  E(s) is stable - that is, all the zeros of E(s) have
%                 non-positive real parts;
%        this corresponds to computing the spectral factorization of the
%        real polynomial A(s) arising from continuous optimality problems.
%        
%        [E,RES,B] = POLYSF(A,DISCR,FORM)  has additional input and output
%        parameters.
%        
%        DISCR specifies the computations to be performed:
%        DISCR = 0 :  compute the spectral factorization for the
%                     continuous optimality problems;
%        DISCR = 1 :  compute the spectral factorization for the
%                     discrete optimality problems.
%        For DISCR = 1, POLYSF computes a real polynomial E(z) such that
%            (a)  E(1/z) * E(z) = A(1/z) * A(z) =: B(z), and
%            (b)  E(z) is stable - that is, E(z) has no zeros with modulus
%                 greater than 1;
%        this corresponds to computing the spectral factorization of the
%        real polynomial A(z) arising from discrete optimality problems.
%        Default: DISCR = 0.
%
%        FORM indicates whether the coefficients of A or B are
%        supplied, as follows:
%        FORM = 1 :  the coefficients of A are supplied;
%        FORM = 2 :  the coefficients of B are supplied.
%        Default: FORM = 1.
%        
%        E is a vector containing the coefficients of the spectral
%        factor E in increasing powers of s or z.
%        
%        RES is an estimate of the accuracy with which the coefficients
%        of the polynomial E have been computed.
%        
%        B is a vector containing the coefficients of the polynomial B,
%        B(s) = b(0) + b(1) * s**2  + ... + b(DA) * s**(2*DA), 
%        if DISCR = 0, or
%        B(z) = b(0) + b(1) * (z + 1/z) + ... + b(DA) * (z**DA + 1/z**DA),
%        if DISCR = 1, where DA is the degree of the polynomial A.
%
%        See also SPECFACT
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Dec-22-2003.
%        Revised Mar-03-2009.
%

%
ni = nargin;
if ni < 1 || nargout < 1,  
    error( ['Usage: E = POLYSF(A)',    sprintf('\n'), ...
            '       E = POLYSF(A,DISCR)', sprintf('\n'), ...
            '       [E,RES,B] = POLYSF(A,DISCR,FORM)' ] );
end
%
if ni == 1,  
    discr = 0;  form = 1;
elseif ni == 2,  
    form = 1;
end
task = discr + 1;
%
[ E, res, B ] = specfact( task, A, form );
%
% end polysf
