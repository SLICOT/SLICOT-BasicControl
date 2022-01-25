function y = slHart( x, scr, pad0 )
%SLHART  Computes the discrete Hartley transform of a real signal.
% 
%        Y = SLHART(X)  computes the discrete Hartley transform of
%        the signal X. The signal is not scrambled.
%        If length(X) is not a power of 2, the trailing part of X
%        is truncated so that the truncated X has a length equal to
%        the largest power of 2.
%
%        Y = SLHART(X,SCR,PAD0)  has additional input arguments:
%
%        SCR specifies how the signal should be processed:
%        SCR = 0:  the signal is not scrambled;
%        SCR = 1:  the input  signal is bit-reversed;
%        SCR = 2:  the output signal is bit-reversed.
%        Default:  SCR = 0.
%
%        PAD0 specifies how sequences whose length is not a power of 2
%        should be dealt with:
%        PAD0 = 0: truncate the trailing part and use a number of data
%                  points equal to the largest power of 2;
%        PAD0 = 1: pad the trailing part with 0 till the next power of 2.
%        Default:  PAD0 = 0.
%
%        See also SLDFT, SLIDFT
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
    error('Usage: Y = slHart(X,SCR,PAD0)')
elseif ni == 1,
    scr  = 0;
    pad0 = 0;
elseif ni == 2,
    pad0 = 0;
end
%
if isreal( x ) == 1,
    if scr == 0,
        job = 0;
    elseif scr == 1,
        job = -6;
    else
        job = 6;
    end
    y = datana(job,x,pad0);
else
    error('A should be a real vector')
end
%
% end slHart
