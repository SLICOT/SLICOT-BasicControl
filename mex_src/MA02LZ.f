      SUBROUTINE MA02LZ( UPLO, M, N, I, J, ALPHA, BETA, A, LDA )
C
C     SLICOT RELEASE 5.6.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To set the (M,N) submatrix starting in the row I and column J
C     of the matrix A to BETA on the diagonal and ALPHA elsewhere.
C
C     This is a wrapper to the LAPACK Library routine ZLASET, used for
C     calls involving %VAL constructs.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Specifies the part of the submatrix of A to be set, as
C             follows:
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
C     I,     (input) INTEGER
C     J      The row and column indices, respectively, of the submatrix
C            of A.
C
C     ALPHA  (input) DOUBLE PRECISION
C            The constant to which the offdiagonal elements are to be
C            set.
C
C     BETA   (input) DOUBLE PRECISION
C            The constant to which the diagonal elements are to be set.
C
C     A      (input) DOUBLE PRECISION array, dimension (LDA,*)
C            The given matrix A.  If JOB = 'U', only part of the upper
C            triangle or trapezoid is accessed; if JOB = 'L', only part
C            of the lower triangle or trapezoid is accessed.
C
C     LDA    INTEGER
C            The leading dimension of the array A.
C            LDA >= max(1,M+I-1).
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
      INTEGER            I, J, LDA, M, N
      COMPLEX*16         ALPHA, BETA
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
C
C     .. External Subroutines ..
      EXTERNAL           ZLASET
C
C     .. Executable Statements ..
C
C     For efficiency, the input parameters are not checked.
C
      CALL ZLASET( UPLO, M, N, ALPHA, BETA, A(I,J), LDA )
C
      RETURN
C *** Last line of MA02LZ ***
      END
