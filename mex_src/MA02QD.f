      INTEGER FUNCTION MA02QD( N, ALPHAR, ALPHAI, BETA )
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To count the number of complex conjugate eigenvalues of a real
C     matrix pencil, given their real and imaginary parts, as well as
C     the denominator.
C
C     This is useful for calls involving %VAL constructs.
C
C     FUNCTION VALUE
C
C     MA02QD  INTEGER
C             The number of complex conjugate eigenvalues.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of eigenvalues.  N >= 0.
C
C     ALPHAR  (input) DOUBLE PRECISION array, dimension (N)
C             The real parts of each scalar alpha defining an eigenvalue
C             of the pencil.
C
C     ALPHAI  (input) DOUBLE PRECISION array, dimension (N)
C             The imaginary parts of each scalar alpha defining an
C             eigenvalue of the pencil.
C             If ALPHAI(j) is zero, then the j-th eigenvalue is real.
C
C     BETA    (input) DOUBLE PRECISION array, dimension (N)
C             The scalars beta that define the eigenvalues of the
C             pencil.
C             Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
C             beta = BETA(j) represent the j-th eigenvalue of the pencil
C             in the form lambda = alpha/beta. Since lambda may overflow
C             the ratios should not, in general, be computed.
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2020.
C
C     REVISIONS
C
C     -
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C
C     .. Scalar Arguments ..
      INTEGER            N
      DOUBLE PRECISION   ALPHAI( * ), ALPHAR( * ), BETA( * )
C
C     .. Local Arguments ..
      INTEGER            I, J
C
C     .. Executable Statements ..
C
C     For efficiency, the input parameters are not checked.
C
      I = 0
      J = 1
C     WHILE( J.LE.N ) DO
   10 CONTINUE
      IF ( J.LE.N ) THEN
         IF( ALPHAR( J ).NE.ZERO .AND. BETA( J ).NE.ZERO .AND.
     $       ALPHAI( J ).NE.ZERO ) THEN
            I = I + 2
            J = J + 2
         ELSE
            J = J + 1
         END IF
         GO TO 10
      END IF
C     END WHILE 10
C
      MA02QD = I
C
      RETURN
C *** Last line of MA02QD ***
      END
