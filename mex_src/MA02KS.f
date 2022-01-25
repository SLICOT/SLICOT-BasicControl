      SUBROUTINE MA02KS( N, A, LDA )
C
C     SLICOT RELEASE 5.6.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To symmetrize an N-by-N matrix.
C
C     This is used for calls involving %VAL constructs.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N      (input) INTEGER
C            The order of the matrix A.  N >= 0.
C
C     A      (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C            On entry, the given full matrix A.
C            On exit, the symmetrized matrix.
C
C     LDA    INTEGER
C            The leading dimension of the array A.
C            LDA >= max(1,N).
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2015.
C
C     REVISIONS
C
C     -
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  HALF, ONE
      PARAMETER         ( HALF = 0.5D0, ONE = 1.0D0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            LDA, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C     .. Local Scalars ..
      INTEGER            I
C     ..
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DSCAL
C     .. Executable Statements ..
C
C     For efficiency, the input parameters are not checked.
C
      DO 10 I = 1, N-1
         CALL DAXPY( N-I, ONE, A(I,I+1), LDA, A(I+1,I), 1 )
         CALL DSCAL( N-I, HALF, A(I+1,I), 1 )
         CALL DCOPY( N-I, A(I+1,I), 1, A(I,I+1), LDA )
   10 CONTINUE
C
      RETURN
C *** Last line of MA02KS ***
      END
