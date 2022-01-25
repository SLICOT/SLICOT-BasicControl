function S = sncs(A,T,job)
%SNCS    Computes the sine or cosine transform of a real signal using
%        the formulas implemented in SLICOT routine DF01MD. 
% 
%        S = SNCS(A,T)  computes the sine transform of the real 
%        signal A, with the sampling time T. If T is omitted,
%        T = 1 is assumed.
% 
%        S = SNCS(A,T,JOB)  computes the transform specified by JOB:
%        sine transform,   if JOB = 0 (default), or 
%        cosine transform, if JOB = 1.
% 
%        The length N of A should satisfy N = 2^m + 1, m >= 2.

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, Jan. 2003.
%
%        Revisions:
%        Mar. 2009.

ni = nargin;
%
if ni < 1,
    error('Usage: S = sncs(A,T,JOB)')
elseif ni == 1,
    T   = 1;
    job = 0;
elseif ni == 2,
    job = 0;
end
%
n = length( A );
if n < 5,
    error('length(A) must be at least 5')
end
%
if ~( n == 2^log2( n - 1 ) + 1 ),
    error('length(A) must be of the form 2^m+1')
end
%
if job == 0,
    % 
    % Sine transform.  Transform to a complex B. 
    % 
    Br = [ -2*A(2); A(2:2:n-3)-A(4:2:n); 2*A(n-1) ];  Bi = [ 0; -A(3:2:n-1); 0 ]; 
    %
    % Inverse Fourier transform of B.
    %
    [Yr,Yi] = datana( -5,Br,Bi );
    %
    D(1:2:n-1) = (n-1)*Yr(1:(n-1)/2);  D(2:2:n) = (n-1)*Yi(1:(n-1)/2);
    %
    % Sine transform coefficients.
    %
    S(1) = 0;  t = pi*(1:n-2)/(n-1);
    S(2:n-1) = T*( ( D(2:n-1) - D(n-1:-1:2) ) - ...
                   ( D(2:n-1) + D(n-1:-1:2) )./( 2*sin( t ) ) );
    S(n) = 0;
    %
else
    %
    % Cosine transform.  Transform to a complex B. 
    % 
    Br = 2*[ A(1:2:n-1); A(n) ];  Bi = 2*[ 0; A(2:2:n-3)-A(4:2:n); 0 ]; 
    %
    % Inverse Fourier transform of B.
    % 
    [Yr,Yi] = datana( -5,Br,Bi );  
    %
    D(1:2:n-1) = (n-1)*Yr(1:(n-1)/2);  D(2:2:n) = (n-1)*Yi(1:(n-1)/2);
    %
    % Cosine transform coefficients.
    % 
    A0   = 2*sum( A(2:2:n-1) );
    S(1) = 2*T*( D(1) + A0 );  t = pi*(1:n-2)/(n-1);
    S(2:n-1) = T*( ( D(2:n-1) + D(n-1:-1:2) ) - ...
                   ( D(2:n-1) - D(n-1:-1:2) )./( 2*sin( t ) ) );
    S(n) = 2*T*( D(1) - A0 );
    %
end
S = S(:);
%
% end sncs
