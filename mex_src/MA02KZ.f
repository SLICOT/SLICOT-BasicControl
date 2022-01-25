      SUBROUTINE MA02KZ( UPLO, M, N, I1, J1, I2, J2, A, LDA, B, LDB )
C
C     SLICOT RELEASE 5.6.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To copy the (M,N) submatrix starting in the row I1 and column J1
C     of the real matrix A to the complex matrix B starting in the row
C     I2 and column J2.
C
C     This is used for calls involving %VAL constructs.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Specifies the part of the submatrix of A to be copied
C             into B, as follows:
C             = 'U': Upper triangular part;
C             = 'L': Lower triangular part;
C             Otherwise:  All of the submatrix of A.
C
C     Input/Output Parameters
C
C     M      (input) INTEGER
C            The number of rows of the submatrix of A.  M >= 0.
C
C     N      (input) INTEGER
C            The number of columns of the submatrix of A.  N >= 0.
C
C     I1,    (input) INTEGER
C     J1     The row and column indices, respectively, of the submatrix
C            of A.
C
C     I2,    (input) INTEGER
C     J2     The row and column indices, respectively, of the submatrix
C            of B.
C
C     A      (input) DOUBLE PRECISION array, dimension (LDA,*)
C            The given matrix A.  If JOB = 'U', only part of the upper
C            triangle or trapezoid is accessed; if JOB = 'L', only part
C            of the lower triangle or trapezoid is accessed.
C
C     LDA    INTEGER
C            The leading dimension of the array A.
C            LDA >= max(1,M+I1-1).
C
C     B      (output) COMPLEX*16 array, dimension (LDB,*)
C            real(B) = A in the locations specified by UPLO, I1, J1, I2,
C            and J2.
C
C     LDB    INTEGER
C            The leading dimension of the array B.
C            LDB >= max(1,M+I2-1).
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2012.
C
C     REVISIONS
C
C     V. Sima, Nov. 2012.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            I1, I2, J1, J2, LDA, LDB, M, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
      COMPLEX*16         B( LDB, * )
C     ..
C     .. Local scalars ..
      INTEGER            IA, IB, JA, JB
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, MIN
C     ..
C     .. Executable Statements ..
C
C     For efficiency, the input parameters are not checked.
C
      IA = I1 - 1
      IB = I2 - 1
      JA = J1 - 1
      JB = J2 - 1
C
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( IB+I, JB+J ) = DCMPLX( A( IA+I, JA+J ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( IB+I, JB+J ) = DCMPLX( A( IA+I, JA+J ) )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( IB+I, JB+J ) = DCMPLX( A( IA+I, JA+J ) )
   50       CONTINUE
   60    CONTINUE
      END IF
C
      RETURN
C *** Last line of MA02KZ ***
      END
