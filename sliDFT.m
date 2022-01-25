function y = sliDFT( x, pad0 )
%SLIDFT  Computes the inverse discrete Fourier transform (DFT) of a signal
%        using Fast Fourier transform.
% 
%        Y = SLIDFT(X)  computes the inverse DFT of the signal X.
%        If length(X) is not a power of 2, the trailing part of X
%        is truncated so that the truncated X has a length equal to
%        the largest power of 2.
%
%        Y = SLIDFT(X,PAD0)  has an additional input argument:
%
%        PAD0 specifies how sequences whose length is not a power of 2
%        should be dealt with:
%        PAD0 = 0: truncate the trailing part and use a number of data
%                  points equal to the largest power of 2;
%        PAD0 = 1: pad the trailing part with 0 till the next power of 2.
%        Default:  PAD0 = 0.
%
%        See also SLDFT
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, Dec. 2002.
%        Revised: Mar. 2009.
%

ni = nargin;
%
if ni < 1,
    error('Usage: Y = sliDFT(X,PAD0)')
elseif ni == 1,
    pad0 = 0;
end
%
if isreal( x ) == 1,
    n = length( x );
    if pad0 == 0,
        xr = x(1:2:n-1);  xi = x(2:2:n);
    else
        xr = x(1:2:n);  xi = x(2:2:n);  if mod( n,2 ) ~= 0,  xi(n+1) = 0;  end
    end    
    [ yr, yi ] = datana( -5, xr, xi, pad0 );
    no = 2*length( yr );  y(1:2:no) = yr;  y(2:2:no) = yi;  y = y(:);  
else
    xr = real( x );  xi = imag( x );
    [ yr, yi ] = datana( -4, xr, xi, pad0 );
    y = yr + 1i*yi;
end
%
% end sliDFT
