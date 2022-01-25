      SUBROUTINE MA02KD( UPLO, M, N, I1, J1, I2, J2, A, LDA, B, LDB )
C
C     SLICOT RELEASE 5.6.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To copy the (M,N) submatrix starting in the row I1 and column J1
C     of the matrix A to the matrix B starting in the row I2 and
C     column J2.
C
C     This is a wrapper to the LAPACK Library routine DLACPY, used for
C     calls involving %VAL constructs.
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
C     B      (output) DOUBLE PRECISION array, dimension (LDB,*)
C            B = A in the locations specified by UPLO, I1, J1, I2, J2.
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
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
C     .. External Subroutines ..
      EXTERNAL           DLACPY
C
C     .. Executable Statements ..
C
C     For efficiency, the input parameters are not checked.
C
      CALL DLACPY( UPLO, M, N, A(I1,J1), LDA, B(I2,J2), LDB )
C
      RETURN
C *** Last line of MA02KD ***
      END
