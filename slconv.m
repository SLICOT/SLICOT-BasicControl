function c = slconv( a, b, meth, pad0 )
%SLCONV  Computes the convolution of two real signals of equal length
%        using either FFT or Hartley transform.
% 
%        C = SLCONV(A,B)  computes the convolution of the signals A
%        and B using Hartley transform.
%        If length(A) is not a power of 2, the trailing parts of A
%        and B are truncated so that the truncated signals have a length
%        equal to the largest power of 2.
%
%        C = SLCONV(A,B,METH,PAD0)  has additional input arguments:
%
%        METH specifies the transform to be used.
%        METH = 0: Hartley transform;
%        METH = 1: FFT transform.
%        Default:  METH = 0.
%
%        PAD0 specifies how sequences whose length is not a power of 2
%        should be dealt with:
%        PAD0 = 0: truncate the trailing parts and use a number of data
%                  points equal to the largest power of 2;
%        PAD0 = 1: pad the trailing parts with 0 till the next power of 2.
%        Default:  PAD0 = 0.
%
%        See also SLDECONV, SLDFT, SLIDFT, SLHART
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, Dec. 2002.
%        Revised: Mar. 2009.
%

ni = nargin;
%
if ni < 2,
    error('Usage: C = slconv(A,B,METH,PAD0)')
elseif ni == 2,
    meth = 0;
    pad0 = 0;
elseif ni == 3,
    pad0 = 0;
end
%
if isreal( a ) == 1 && isreal( b ) == 1 ,
    if meth == 0, 
        c = datana( 2, a, b, pad0 );
    else 
        c = datana( 1, a, b, pad0 );
    end
else
    error('A and B should be real vectors')
end
%
% end slconv
