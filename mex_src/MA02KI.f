      SUBROUTINE MA02KI( N, I, INCA, J, INCB, A, B )
C
C     SLICOT RELEASE 5.6.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To copy N entries starting at the index I of the vector A, with
C     increment INCA, to the vector B starting at the index J, with
C     increment INCB.
C
C     This is used for calls involving %VAL constructs.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N      (input) INTEGER
C            The number of entries to copy.  N >= 0.
C
C     I      (input) INTEGER
C            The starting index in the vector A.
C
C     INCA   (input) INTEGER
C            The increment between the entries of the vector A.
C
C     J      (input) INTEGER
C            The starting index in the vector B.
C
C     INCB   (input) INTEGER
C            The increment between the entries of the vector B.
C
C     A      (input) DOUBLE PRECISION array, dimension (*)
C            The given vector A.
C
C     B      (output) DOUBLE PRECISION array, dimension (*)
C            The resulted vector B:
C            B(J) = A(I), B(J+INCB) = A(I+INCA), ..., 
C            B(J+(N-1)*INCB) = A(I+(N-1)*INCA).
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2013.
C
C     REVISIONS
C
C     -
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER            I, INCA, INCB, J, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( * ), B( * )
C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY
C
C     .. Executable Statements ..
C
C     For efficiency, the input parameters are not checked.
C
      CALL DCOPY( N, A(I), INCA, B(J), INCB )
C
      RETURN
C *** Last line of MA02KI ***
      END
