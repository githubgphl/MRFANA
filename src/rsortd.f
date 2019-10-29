C     Singleton, Richard C. "Algorithm 347: an efficient algorithm for
C     sorting with minimal storage [M1]." Communications of the ACM 12.3
C     (1969): 185-186.
C=======================================================================--------
      SUBROUTINE   R SORT D ( A, N, TAG )
C
C        PURPOSE
C
C          SORT A REAL VECTOR INTO DECREASING ORDER FROM A(1) TO A(N).
C          VECTOR TAG IS PERMUTED THE SAME AS VECTOR A.
C
C        DESCRIPTION OF ARGUMENTS
C
C           A    - THE N-VECTOR OF REAL NUMBERS TO BE SORTED.
C           N    - THE NUMBER OF ELEMENTS IN THE VECTOR TO BE SORTED.
C           TAG  - THE N-VECTOR CONTAINING THE TAG FIELDS.
C
C        REMARKS
C
C           THE PROCEDURE REQUIRES TWO ADDITIONAL ARRAYS IU(K) AND
C           IL(K) WHICH PERMIT SORTING UP TO 2**(K+1)-1 ELEMENTS.
C           THESE ARRAYS ARE SUPPLIED BY THE SUBROUTINE WITH K = 32.
C           THIS ALLOWS SORTING A MAXIMUM OF 8589934591 ELEMENTS.
C
C        EXAMPLE
C
C           IF N RANDOM NUMBERS, UPON GENERATION, WERE ASSIGNED
C           CONSECUTIVE TAGS OF 1 THROUGH N, AND THE ARRAY OF RANDOM
C           NUMBERS SORTED, THE TAG ARRAY COULD THEN BE USED TO
C           DETERMINE THAT THE RANDOM NUMBER NOW IN POSITION A(I) WAS
C           ORIGINALLY IN POSITION A(J) BECAUSE TAG(I) = J.
C
C        METHOD
C
C           -AN EFFICIENT ALGORITHM FOR SORTING WITH MINIMAL STORAGE-
C           BY RICHARD C. SINGLETON
C           PREPARED FOR:
C           INSTITUTE RESEARCH AND DEVELOPMENT
C           STANFORD RESEARCH INSTITUTE PROJECT 387531-132
C           SEPTEMBER 1968.
C
C
C  ..................................................................
C
C
      IMPLICIT NONE
C
C     Argument variables
C
      INTEGER N
C
      REAL A(N)
C
C     Local variables
C
      INTEGER I, IJ, IL(32), IU(32), J, K, L, M
C
      REAL T, TT
C
      INTEGER      TAG(N), TG
C
      TT=0
      M=1
      I=1
      J=N
    5 IF(I .GE. J) GO TO 70
   10 K=I
      IJ=(J+I)/2
      T=A(IJ)
      IF(A(I) .GE. T) GO TO 20
      A(IJ)=A(I)
      A(I)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(I)
      TAG(I)=TG
   20 L=J
      IF(A(J) .LE. T) GO TO 40
      A(IJ)=A(J)
      A(J)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(J)
      TAG(J)=TG
      IF(A(I) .GE. T) GO TO 40
      A(IJ)=A(I)
      A(I)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(I)
      TAG(I)=TG
      GO TO 40
   30 A(L)=A(K)
      A(K)=TT
      TG=TAG(L)
      TAG(L)=TAG(K)
      TAG(K)=TG
   40 L=L-1
      IF(A(L) .LT. T) GO TO 40
      TT=A(L)
   50 K=K+1
      IF(A(K) .GT. T) GO TO 50
      IF(K .LE. L) GO TO 30
      IF(L-I .LE. J-K) GO TO 60
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GO TO 80
   60 IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GO TO 80
   70 M=M-1
      IF(M.EQ.0) RETURN
      I=IL(M)
      J=IU(M)
   80 IF(J-I.GE.1) GOTO 10
      IF(I.EQ.1) GOTO 5
      I=I-1
   90 I=I+1
      IF(I .EQ. J) GO TO 70
      T=A(I+1)
      IF(A(I) .GE. T) GO TO 90
      TG=TAG(I+1)
      K=I
  100 A(K+1)=A(K)
      TAG(K+1)=TAG(K)
      K=K-1
      IF(T .GT. A(K)) GO TO 100
      A(K+1)=T
      TAG(K+1)=TG
      GO TO 90
      END
