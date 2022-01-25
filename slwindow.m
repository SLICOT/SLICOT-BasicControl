function y = slwindow( x, wndw )
%SLWINDOW  Computes the anti-aliasing window applied to a real signal.
% 
%        Y = SLWINDOW(X)  computes the anti-aliasing Hamming window.
%
%        Y = SLWINDOW(X,WNDW)  has an additional input argument:
%
%        WNDW specifies how the signal should be processed:
%        WNDW = 1:  Hamming window;
%        WNDW = 2:  Hann window;
%        WNDW = 3:  Quadratic window.
%        Default:   WNDW = 1.
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
    error('Usage: Y = slwindow(X,WNDW)')
elseif ni == 1,
    wndw = 1;
end
%
if isreal( x ) == 1,
    y = datana(7,x,wndw);
else
    error('X should be a real vector')
end
%
% end slwindow
