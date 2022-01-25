function s = sincos( a, t, job, pad0 )
%SINCOS  Computes the sine transform or cosine transform of a real
%        signal using Fast Fourier transform.
% 
%        S = SINCOS(A)  computes the sine transform of the signal A.
%        If length(A) is not 2^m+1, with integer m, the trailing part
%        of A is truncated so that the truncated A has a length equal
%        to the largest integer 2^m+1.
%
%        S = SINCOS(A,T,JOB,PAD0)  has additional input arguments:
%
%        T specifies the sampling time of the signal.
%        Default:  T = 1.
%
%        JOB specifies which transform should be computed.
%        JOB = 0:  compute the sine transform;
%        JOB = 1:  compute the cosine transform.
%        Default:  JOB = 0.
%
%        PAD0 specifies how sequences whose length is not 2^m+1 should
%        be dealt with:
%        PAD0 = 0: truncate the trailing part and use a number of data
%                  points equal to 2^m+1;
%        PAD0 = 1: pad the trailing part with 0 till the next value 2^m+1.
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
    error('Usage: S = sincos(A,T,JOB,PAD0)')
elseif ni == 1,
    t    = 1;
    job  = 0;
    pad0 = 0;
elseif ni == 2,
    job  = 0;
    pad0 = 0;
elseif ni == 3,
    pad0 = 0;
end
%
if isreal( a ) == 1,
    if job == 0, 
        s = datana( -3, a, t, pad0 );
    else 
        s = datana( 3, a, t, pad0 );
    end
else
    error('A should be a real vector')
end
%
% end sincos
