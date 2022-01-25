C DATANA.F - Gateway function for data analysis calculations using
C            SLICOT routines DE01OD, DE01PD, DF01MD, DG01MD, DG01ND,
C            DG01OD, and DK01MD.
C
C RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [C(,D)] = datana(job,A(,B)(,T)(,window)(,pad0))
C
C   [C]     = datana(job,A(,pad0))        |job| = 0, 6;
C   [C]     = datana(job,A,B(,pad0))      |job| = 1, 2;
C   [C]     = datana(job,A(,T)(,pad0))    |job| = 3;
C   [YR,YI] = datana(job,XR,XI(,pad0))    |job| = 4, 5;
C   [C]     = datana(job,A(,window))       job  = 7.
C
C Purpose:
C  To perform various transforms of real or complex vectors, used in
C  data analysis calculations.
C
C Input parameters:
C   job    - option parameter indicating the task to be performed.
C            =-6 :  scrambled discrete Hartley transform of a real
C                   signal A (the input signal is bit-reversed);
C            =-5 :  inverse discrete Fourier transform of a real signal;
C            =-4 :  inverse discrete Fourier transform of a complex
C                   signal XR+i*XI;
C            =-3 :  sine transform of a real signal A;
C            =-2 :  deconvolution of two real signals A and B using
C                   Hartley transform;
C            =-1 :  deconvolution of two real signals A and B using FFT;
C            = 0 :  discrete Hartley transform of a real signal A
C                   (the signal is not scrambled);
C            = 1 :  convolution of two real signals A and B using FFT,
C                   defined in MATLAB by real(ifft(fft(A).*fft(B)));
C            = 2 :  convolution of two real signals A and B using
C                   Hartley transform;
C            = 3 :  cosine transform of a real signal A;
C            = 4 :  discrete Fourier transform of a complex signal
C                   XR+i*XI;
C            = 5 :  discrete Fourier transform of a real signal X;
C            = 6 :  scrambled discrete Hartley transform of a real
C                   signal A (the output transform is bit-reversed);
C            = 7 :  anti-aliasing window applied to a real signal A.
C   A      - the n-vector A.
C   B      - (optional) if |job| = 1 or 2, the n-vector B.
C   T      - (optional) if |job| = 3, the sampling time of the
C            signal; otherwise, it is not used.
C            Default:  T = 1.
C   XR     - if job =  4, the n-vector containing the real part of the
C                         complex signal X.
C            if job = -4, the n-vector containing the real part of the
C                         discrete Fourier transform.
C            if job =  5, the n-vector containing the odd part of the
C                         real signal X.
C            if job = -5, the n+1-vector containing the real part of
C                         the discrete Fourier transform.
C   XI     - if job =  4, the n-vector containing the imaginary part of
C                         the complex signal X.
C            if job = -4, the n-vector containing the imaginary part of
C                         the discrete Fourier transform.
C            if job =  5, the n-vector containing the even part of the
C                         real signal X.
C            if job = -5, the n+1-vector containing the imaginary
C                         part of the discrete Fourier transform.
C   window - (optional) if job = 7, integer specifying the type of
C            window to use:
C            = 1: Hamming window;
C            = 2: Hann window;
C            = 3: Quadratic window.
C            Default:  window = 1.
C   pad0   - (optional) if job <> 7, integer specifying how sequences
C            whose length is not a power of 2 should be dealt with:
C            = 0: truncate the trailing part and use a number of data
C                 points corresponding to the largest power of 2, i.e.,
C                 2^m, if |job| <> 3, or 2^m+1, if |job| = 3;
C            = 1: pad the trailing part with 0 till the next power
C                 of 2 (plus 1, if |job| = 3).
C            Default:  pad0 = 0.
C
C Output parameters:
C   C      - the n-vector of results computed according to job.
C   YR     - if job =  4, the n-vector containing the real part of
C                         the computed discrete Fourier transform.
C            if job = -4, the n-vector containing the real part of
C                         inverse discrete Fourier transform.
C            if job =  5, the n+1-vector containing the real part of
C                         the discrete Fourier transform.
C            if job = -5, the odd part of the inverse discrete
C                         Fourier transform.
C   YI     - if job =  4, the n-vector containing the imaginary part
C                         of the computed discrete Fourier transform.
C            if job = -4, the n-vector containing the imaginary part
C                         of the inverse discrete Fourier transform.
C            if job =  5, the n+1-vector containing the imaginary part
C                         of the discrete Fourier transform.
C            if job = -5, the even part of the inverse discrete
C                         Fourier transform.
C
C Further Comments:
C   1) Except for job = 7, this function essentially works on signals
C      whose length is a power of 2, 2^m (or 2^m+1, if |job| = 3).
C      For |job| = 3, m >= 2.
C   2) For job = 5, this function computes the first n+1 elements
C      of the discrete Fourier transform.  The remaining n-1 elements
C      can be obtained using the MATLAB command conj( Y(n:-1:2) ),
C      where Y = YR + i*YI.
C   3) For job = -5, this function uses as input the first
C      n+1 elements of the discrete Fourier transform.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2002.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009,
C   Dec. 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C
C .. Mex-file interface parameters ..
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         CONV, INDI, SCR, SICO, WGHT, WIND
      INTEGER           INFO, N
      DOUBLE PRECISION  DT
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE :: A(:), B(:), DWORK(:), W(:),
     $                                 XI(:), XR(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           I, IJOB, IP, IPAD, IWIN, JOB, LA, M, NZ
      DOUBLE PRECISION  DUM(1), TEMP
C
C .. External subroutines ..
      EXTERNAL          DCOPY,  DE01OD, DE01PD, DF01MD, DG01MD, DG01ND,
     $                  DG01OD, DK01MD, DSCAL
C
C .. Intrinsic functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'DATANA requires at least 2 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'DATANA requires at least 1 output argument' )
      END IF
C
C   job, A(n), (B(n)) (,T) (,window) (,pad0).
C
      M = mxGetM( PRHS(2) )
      N = mxGetN( PRHS(2) )
      I = MIN( M, N )
      N = MAX( M, N )
      M = I
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR.
     $     mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      JOB  = TEMP
      IJOB = ABS( JOB )
      IF ( JOB.LT.-6 ) THEN
         CALL mexErrMsgTxt
     $           ( 'JOB must be larger than or equal to -6' )
      ELSE IF ( JOB.GT.7 ) THEN
         CALL mexErrMsgTxt
     $           ( 'JOB must be less than or equal to 7' )
      END IF
C
C Recheck for proper number of arguments.
C
      IF ( ( IJOB.EQ.1 .OR. IJOB.EQ.2 .OR. IJOB.EQ.4 .OR. IJOB.EQ.5 )
     $     .AND. NRHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'DATANA requires at least 3 input arguments' )
      END IF
      IF ( ( IJOB.EQ.4 .OR. IJOB.EQ.5 ) .AND. NLHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'DATANA requires 2 output arguments' )
      END IF
C
      IF ( ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(2) ).EQ.1 ) .OR. M.NE.1 ) THEN
         IF ( IJOB.EQ.4 .OR. IJOB.EQ.5 ) THEN
            CALL mexErrMsgTxt( 'XR must be a real vector' )
         ELSE
            CALL mexErrMsgTxt( 'A must be a real vector' )
         END IF
      END IF
C
      IP = 3
      IF ( IJOB.EQ.1 .OR. IJOB.EQ.2 ) THEN
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(3) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'B must be a real vector' )
         END IF
         IF ( MAX( mxGetM( PRHS(3) ), mxGetN( PRHS(3) ) ).NE.N .AND.
     $        MIN( mxGetM( PRHS(3) ), mxGetN( PRHS(3) ) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'B must have the same size as A' )
         END IF
         IP = 4
C
      ELSE IF ( IJOB.EQ.4 .OR. IJOB.EQ.5 ) THEN
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(3) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'XI must be a real vector' )
         END IF
         IF ( MAX( mxGetM( PRHS(3) ), mxGetN( PRHS(3) ) ).NE.N .AND.
     $        MIN( mxGetM( PRHS(3) ), mxGetN( PRHS(3) ) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'XI must have the same size as XR' )
         END IF
         IP = 4
C
      ELSE IF ( IJOB.EQ.3 ) THEN
         IF ( NRHS.GE.3 ) THEN
C
C   T
C
            IF ( mxGetM( PRHS(3) ).NE.1 .OR.
     $           mxGetN( PRHS(3) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'T must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(3) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'T must be a real scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), DT, 1 )
            IP = IP + 1
         ELSE
            DT = ONE
         END IF
C
      ELSE IF ( JOB.EQ.7 ) THEN
         IF ( NRHS.EQ.3 ) THEN
C
C   window
C
            IF ( mxGetM( PRHS(3) ).NE.1 .OR.
     $           mxGetN( PRHS(3) ).NE.1 ) THEN
               CALL mexErrMsgTxt( 'WINDOW must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(3) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'WINDOW must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TEMP, 1 )
            IWIN = TEMP
            IF ( IWIN.LT.0 .OR. IWIN.GT.3 ) THEN
               CALL mexErrMsgTxt
     $            ( 'WINDOW has 1, 2, or 3 the only admissible values' )
            END IF
         ELSE
            IWIN = 1
         END IF
      END IF
C
      IF ( JOB.NE.7 .AND. NRHS.GE.IP ) THEN
C
C   pad0
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'PAD0 must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'PAD0 must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IPAD = TEMP
         IF ( IPAD.LT.0 .OR. IPAD.GT.1 ) THEN
            CALL mexErrMsgTxt
     $         ( 'PAD0 has 0 or 1 the only admissible values' )
         END IF
      ELSE
         IPAD = 0
      END IF
C
      IF ( IJOB.EQ.1 .OR. IJOB.EQ.2 ) THEN
         IF ( JOB.LT.0 ) THEN
            CONV = 'D'
         ELSE
            CONV = 'C'
         END IF
         IF ( IJOB.EQ.2 )
     $      WGHT = 'N'
      ELSE IF ( IJOB.EQ.3 ) THEN
         IF ( JOB.LT.0 ) THEN
            SICO = 'S'
         ELSE
            SICO = 'C'
         END IF
      ELSE IF ( IJOB.EQ.4 .OR. IJOB.EQ.5 ) THEN
         IF ( JOB.LT.0 ) THEN
            INDI = 'I'
         ELSE
            INDI = 'D'
         END IF
         IF ( JOB.EQ.-5 )
     $      N = N - 1
C
      ELSE IF ( IJOB.EQ.6 ) THEN
         IF ( JOB.LT.0 ) THEN
            SCR = 'I'
         ELSE
            SCR = 'O'
         END IF
         WGHT = 'N'
      ELSE IF ( IJOB.EQ.0 ) THEN
         SCR  = 'N'
         WGHT = 'N'
C
      ELSE IF ( JOB.EQ.7 ) THEN
         IF ( IWIN.EQ.1 ) THEN
            WIND = 'M'
         ELSE IF ( IWIN.EQ.2 ) THEN
            WIND = 'N'
         ELSE
            WIND = 'Q'
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C First check if for job <> 7, n is power of 2, n = 2^m
C (or 2^m+1, for |job| = 3).
C
      IF ( IJOB.LT.7 ) THEN
         I = 1
         M = 0
         IF ( N.GT.1 ) THEN
C           WHILE 2**M < N DO
   10       CONTINUE
               I = I*2
               M = M + 1
               IF ( I.LT.N )
     $            GO TO 10
C           END WHILE 10
         END IF
         IF ( I.GT.N ) THEN
            IF ( IPAD.EQ.1 ) THEN
               IF ( IJOB.EQ.3 ) THEN
                  IF ( I/2 + 1.NE.N ) THEN
                     LA = I + 1
                  ELSE
                     LA = N
                  END IF
               ELSE
                  LA = I
               END IF
            ELSE
               IF ( IJOB.EQ.3 ) THEN
                  N = I/2 + 1
               ELSE
                  N = I/2
               END IF
               LA = N
            END IF
         ELSE
            IF ( IJOB.EQ.3 ) THEN
               IF ( IPAD.EQ.1 ) THEN
                  LA = I + 1
               ELSE
                  N  = I/2 + 1
                  LA = N
               END IF
            ELSE
               LA = N
            END IF
         END IF
         NZ = LA - N
         IF ( JOB.EQ.-5 .AND. NZ.GT.0 )
     $      NZ = NZ - 1
      ELSE
         NZ = 0
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      IF ( IJOB.EQ.3 .AND. LA.LT.5 ) THEN
         CALL mexErrMsgTxt( 'N must be at least 5' )
      ELSE IF ( ( IJOB.EQ.1 .OR. IJOB.EQ.4 .OR. IJOB.EQ.5 ) .AND.
     $          LA.LT.2 ) THEN
         CALL mexErrMsgTxt( 'N must be at least 2' )
      ELSE IF ( JOB.EQ.7 .AND. N.LT.1 ) THEN
         CALL mexErrMsgTxt( 'N must be at least 1' )
      END IF
C
      IF ( IJOB.EQ.1 ) THEN
         ALLOCATE ( A( LA ), B( LA ) )
      ELSE IF ( IJOB.EQ.2 ) THEN
         ALLOCATE ( A( LA ), B( LA ), W( LA-M ) )
      ELSE IF ( IJOB.EQ.3 ) THEN
         ALLOCATE ( A( LA ), DWORK( LA+1 ) )
      ELSE IF ( IJOB.EQ.4 ) THEN
         ALLOCATE ( XI( LA ), XR( LA ) )
      ELSE IF ( IJOB.EQ.5 ) THEN
         ALLOCATE ( XI( LA+1 ), XR( LA+1 ) )
      ELSE IF ( IJOB.EQ.6 .OR. IJOB.EQ.0 ) THEN
         ALLOCATE ( A( LA ), W( LA-M ) )
      ELSE
         ALLOCATE ( A( N ) )
      END IF
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      M = N
      IF ( JOB.EQ.-5 )
     $   M = N + 1
      IF ( IJOB.LE.3 .OR. IJOB.GE.6 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, N )
      ELSE
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), XR, M )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), XI, M )
      END IF
      IF ( IJOB.EQ.1 .OR. IJOB.EQ.2 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), B, N )
      END IF
C
      IF ( NZ.GT.0 ) THEN
C
C Pad the trailing part with zero and reset N.
C
         DUM(1) = ZERO
         IF ( IJOB.LE.3 .OR. IJOB.GE.6 ) THEN
            CALL DCOPY( NZ, DUM(1), 0, A(N+1), 1 )
         ELSE
            CALL DCOPY( NZ, DUM(1), 0, XR(M+1), 1 )
            CALL DCOPY( NZ, DUM(1), 0, XI(M+1), 1 )
         END IF
         IF ( IJOB.LE.2 ) THEN
            CALL DCOPY( NZ, DUM(1), 0, B(N+1), 1 )
         END IF
         N = LA
      END IF
      M = N
      IF ( JOB.EQ.5 )
     $   M = N + 1
C
C Do the actual computations.
C
      IF ( IJOB.EQ.1 ) THEN
         CALL DE01OD( CONV, N, A, B, INFO )
      ELSE IF ( IJOB.EQ.2 ) THEN
         CALL DE01PD( CONV, WGHT, N, A, B, W, INFO )
      ELSE IF ( IJOB.EQ.3 ) THEN
         CALL DF01MD( SICO, N, DT, A, DWORK, INFO )
      ELSE IF ( IJOB.EQ.4 ) THEN
         CALL DG01MD( INDI, N, XR, XI, INFO )
         IF ( JOB.LT.0 ) THEN
            TEMP = ONE/DBLE( N )
            CALL DSCAL( N, TEMP, XR, 1 )
            CALL DSCAL( N, TEMP, XI, 1 )
         END IF
      ELSE IF ( IJOB.EQ.5 ) THEN
         CALL DG01ND( INDI, N, XR, XI, INFO )
         IF ( JOB.LT.0 ) THEN
            TEMP = ONE/DBLE( 2*N )
            CALL DSCAL( N, TEMP, XR, 1 )
            CALL DSCAL( N, TEMP, XI, 1 )
         END IF
      ELSE IF ( IJOB.EQ.6 .OR. IJOB.EQ.0 ) THEN
         CALL DG01OD( SCR, WGHT, N, A, W, INFO )
      ELSE
         CALL DK01MD( WIND, N, A, INFO )
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( IJOB.LE.3 .OR. IJOB.GE.6 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( N, 1, 0 )
         CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(1) ), N )
      ELSE
         PLHS(1) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( XR, mxGetPr( PLHS(1) ), M )
         PLHS(2) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( XI, mxGetPr( PLHS(2) ), M )
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      IF ( IJOB.EQ.1 ) THEN
         DEALLOCATE ( A, B )
      ELSE IF ( IJOB.EQ.2 ) THEN
         DEALLOCATE ( A, B, W )
      ELSE IF ( IJOB.EQ.3 ) THEN
         DEALLOCATE ( A, DWORK )
      ELSE IF ( IJOB.EQ.4 .OR. IJOB.EQ.5 ) THEN
         DEALLOCATE ( XI, XR )
      ELSE IF ( IJOB.EQ.6 .OR. IJOB.EQ.0 ) THEN
         DEALLOCATE ( A, W )
      ELSE
         DEALLOCATE ( A )
      END IF
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( IJOB.EQ.1 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DE01OD'')'
     $           ) INFO
         ELSE IF ( IJOB.EQ.2 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DE01PD'')'
     $           ) INFO
         ELSE IF ( IJOB.EQ.3 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DF01MD'')'
     $           ) INFO
         ELSE IF ( IJOB.EQ.4 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DG01MD'')'
     $           ) INFO
         ELSE IF ( IJOB.EQ.5 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DG01ND'')'
     $           ) INFO
         ELSE IF ( IJOB.EQ.6 .OR. IJOB.EQ.0 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DG01OD'')'
     $           ) INFO
         ELSE
            WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DK01MD'')'
     $           ) INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of DATANA ***
      END
