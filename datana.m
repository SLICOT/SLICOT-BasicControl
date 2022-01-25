% DATANA.F - MEX-function for data analysis calculations using
%            SLICOT routines DE01OD, DE01PD, DF01MD, DG01MD, DG01ND,
%            DG01OD, and DK01MD.
%
%   [C(,D)] = datana(job,A(,B)(,T)(,window)(,pad0))
%
%   [C]     = datana(job,A(,pad0))        |job| = 0, 6;
%   [C]     = datana(job,A,B(,pad0))      |job| = 1, 2;
%   [C]     = datana(job,A(,T)(,pad0))    |job| = 3;
%   [YR,YI] = datana(job,XR,XI(,pad0))    |job| = 4, 5;
%   [C]     = datana(job,A(,window))       job  = 7.
%
%  DATANA performs various transforms of real or complex vectors, used in
%  data analysis calculations.
%
%   Description of input parameters:
%   job    - option parameter indicating the task to be performed.
%            =-6 :  scrambled discrete Hartley transform of a real
%                   signal A (the input signal is bit-reversed);
%            =-5 :  inverse discrete Fourier transform of a real signal;
%            =-4 :  inverse discrete Fourier transform of a complex
%                   signal XR+i*XI;
%            =-3 :  sine transform of a real signal A;
%            =-2 :  deconvolution of two real signals A and B using
%                   Hartley transform;
%            =-1 :  deconvolution of two real signals A and B using FFT;
%            = 0 :  discrete Hartley transform of a real signal A
%                   (the signal is not scrambled);
%            = 1 :  convolution of two real signals A and B using FFT,
%                   defined in MATLAB by real(ifft(fft(A).*fft(B)));
%            = 2 :  convolution of two real signals A and B using
%                   Hartley transform;
%            = 3 :  cosine transform of a real signal A;
%            = 4 :  discrete Fourier transform of a complex signal
%                   XR+i*XI;
%            = 5 :  discrete Fourier transform of a real signal X;
%            = 6 :  scrambled discrete Hartley transform of a real
%                   signal A (the output transform is bit-reversed);
%            = 7 :  anti-aliasing window applied to a real signal A.
%   A      - the n-vector A.
%   B      - (optional) if |job| = 1 or 2, the n-vector B.
%   T      - (optional) if |job| = 3, the sampling time of the
%            signal; otherwise, it is not used.
%            Default:  T = 1.
%   XR     - if job =  4, the n-vector containing the real part of the
%                         complex signal X.
%            if job = -4, the n-vector containing the real part of the
%                         discrete Fourier transform.
%            if job =  5, the n-vector containing the odd part of the
%                         real signal X.
%            if job = -5, the n+1-vector containing the real part of
%                         the discrete Fourier transform.
%   XI     - if job =  4, the n-vector containing the imaginary part of
%                         the complex signal X.
%            if job = -4, the n-vector containing the imaginary part of
%                         the discrete Fourier transform.
%            if job =  5, the n-vector containing the even part of the
%                         real signal X.
%            if job = -5, the n+1-vector containing the imaginary
%                         part of the discrete Fourier transform.
%   window - (optional) if job = 7, integer specifying the type of
%            window to use:
%            = 1: Hamming window;
%            = 2: Hann window;
%            = 3: Quadratic window.
%            Default:  window = 1.
%   pad0   - (optional) if job <> 7, integer specifying how sequences 
%            whose length is not a power of 2 should be dealt with:
%            = 0: truncate the trailing part and use a number of data
%                 points corresponding to the largest power of 2, i.e.,
%                 2^m, if |job| <> 3, or 2^m+1, if |job| = 3;
%            = 1: pad the trailing part with 0 till the next power
%                 of 2 (plus 1, if |job| = 3).
%            Default:  pad0 = 0.
%
%   Description of output parameters:
%   C      - the n-vector of results computed according to job.
%   YR     - if job =  4, the n-vector containing the real part of
%                         the computed discrete Fourier transform.
%            if job = -4, the n-vector containing the real part of
%                         inverse discrete Fourier transform.
%            if job =  5, the n+1-vector containing the real part of
%                         the discrete Fourier transform.
%            if job = -5, the odd part of the inverse discrete
%                         Fourier transform.
%   YI     - if job =  4, the n-vector containing the imaginary part
%                         of the computed discrete Fourier transform.
%            if job = -4, the n-vector containing the imaginary part
%                         of the inverse discrete Fourier transform.
%            if job =  5, the n+1-vector containing the imaginary part
%                         of the discrete Fourier transform.
%            if job = -5, the even part of the inverse discrete
%                         Fourier transform.
%
% Further Comments:
%   1) Except for job = 7, this function essentially works on signals
%      whose length is a power of 2, 2^m (or 2^m+1, if |job| = 3).  
%      For |job| = 3, m >= 2.  
%   2) For job = 5, this function computes the first n+1 elements
%      of the discrete Fourier transform.  The remaining n-1 elements
%      can be obtained using the MATLAB command conj( Y(n:-1:2) ),
%      where Y = YR + i*YI.
%   3) For job = -5, this function uses as input the first
%      n+1 elements of the discrete Fourier transform.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2002.
%
% Revisions:
%   -
%
