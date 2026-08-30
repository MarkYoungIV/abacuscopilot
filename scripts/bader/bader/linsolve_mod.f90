!**************************************************************
!*       SOLVE A LINEAR SYSTEM BY DIRECT FACTORIZATION        *
!* ---------------------------------------------------------- *
!* SAMPLE RUN:                                                *
!* (Solve linear system A X = B, where:                       *
!*                                                            *
!*       1  0  0  0  0  1              1                      *
!*       1  1  0  0  0 -1              0                      *
!*  A = -1  1  1  0  0  1          B = 1                      *
!*       1 -1  1  1  0 -1              0                      *
!*      -1  1 -1  1  1  1              1                      *
!*       1 -1  1 -1  1 -1              0 )                    *
!*                                                            *
!* LINEAR SYSTEM AX = B:                                      *
!*                                                            *
!*  1.0000  0.0000  0.0000  0.0000  0.0000  1.0000    1.0000  *
!*  1.0000  1.0000  0.0000  0.0000  0.0000 -1.0000    0.0000  *
!* -1.0000  1.0000  1.0000  0.0000  0.0000  1.0000    1.0000  *
!*  1.0000 -1.0000  1.0000  1.0000  0.0000 -1.0000    0.0000  *
!* -1.0000  1.0000 -1.0000  1.0000  1.0000  1.0000    1.0000  *
!*  1.0000 -1.0000  1.0000 -1.0000  1.0000 -1.0000    0.0000  *
!*                                                            *
!* SOLUTION:                                                  *
!*  0.343750000000000                                         *
!*  0.312500000000000                                         *
!*  0.375000000000000                                         *
!*  0.250000000000000                                         *
!*  0.500000000000000                                         *
!*  0.656250000000000                                         *
!*                                                            *
!* ---------------------------------------------------------- *
!* Ref.: From Numath Library By Tuan Dang Trong in Fortran 77 *
!*       [BIBLI 18].                                          *
!*                                                            *
!*                           F90 Release By J-P Moreau, Paris *
!*                                  (www.jpmoreau.fr)         *
!**************************************************************  
!     PROGRAM TEST_DLITTL
!     
!     parameter(N=6)
!     
!     real*8 A(N,N),B(N),X(N),W(N,N),Z(N)
!     
!     DATA A /1.d0, 1.d0,-1.d0, 1.d0,-1.d0, 1.d0, &
!             0.d0, 1.d0, 1.d0,-1.d0, 1.d0,-1.d0, &
!                 0.d0, 0.d0, 1.d0, 1.d0,-1.d0, 1.d0, &
!             0.d0, 0.d0, 0.d0, 1.d0, 1.d0,-1.d0, &
!                 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 1.d0, &
!                     1.d0,-1.d0, 1.d0,-1.d0, 1.d0,-1.d0/              
!     
!     DATA B /1.d0,0.d0,1.d0,0.d0,1.d0,0.d0/
!     
!     print *,' '
!     print *,' LINEAR SYSTEM AX=B:'
!     print *,' '
!     do I=1,N
!       write(*,10) (A(I,J),J=1,N), '   ', B(I)
!     end do
!     
!     call DLITTL(A,B,N,N,X,W,Z)
!     
!     print *,' '
!     print *,' SOLUTION:'
!     do I=1,N
!       write(*,*) X(I)
!     end do
!     
!     print *,' '
!     
!     10 format(6F8.4,A3,F8.4)
!     
!     END

      MODULE linsolve_mod
        IMPLICIT NONE

        PUBLIC :: DLITTL

        CONTAINS 


        SUBROUTINE DLITTL(A,B,LDA,N,X,W,Z)
!----------------------------------------------------------------------
!     SOLVE A LINEAR SYSTEM BY DIRECT FACTORIZATION
!     A*X = L*U*X = B
!     L(LOWER TRIANGULAR MATRIX WITH L(I,I)=1)
!     U(UPPER TRIANGULAR MATRIX),THE SUB-DIAGONAL TERMS OF L AND THE
!     TERMS U(I,J) ARE STORED IN THE SAME MATRIX, W. THE RESOLUTION
!     IS OBTAINED BY SOLVING TWO TRIANGULAR SYSTEMS:
!     L*Z = B  AND U*X = Z.
!     THE PARTIAL PIVOTING IS MADE BY CHOOSING THE BIGGEST ELEMENT OF
!     EACH COLUMN OF THE TRANSFORMATION MATRIX.
!     INPUTS:
!     A = TABLE OF SIZE (LDA,N)                                    R*8
!     B = SECOND MEMBER VECTOR (N)                                 R*8
!     LDA = FIRST DIMENSION OF A AND W IN MAIN                     I*4
!     N = ORDER OF LINEAR SYSTEM                                   I*4
!     OUTPUTS:
!     X = SOLUTION VECTOR (N)                                      R*8
!     WORKING ZONE:
!     W = TABLE OF SIZE (LDA,N)                                    R*8
!     Z = AUXILIARY VECTOR (N)                                     R*8
!     NOTE:
!     MESSAGE '** DLITTL ** NO UNIQUE SOLUTION' IS GIVEN WHEN A 
!     NULL PIVOT IS FOUND.
!-----------------------------------------------------------------------
        INTEGER :: P, I, J, K, N, LDA
        REAL*8 A(LDA,*),W(LDA,*),B(*),X(*),Z(*),SUM,AMAX,T

        P=1
        AMAX=ABS(A(1,1))
!     SEEK MAX. ELEMENT OF FIRST COLUMN
        DO 2 J=2,N
        IF(ABS(A(J,1)).LT.AMAX) GO TO 2
        AMAX=A(J,1)
        P=J
      2 CONTINUE
        IF(AMAX.EQ.0.D0) GO TO 32
        IF(P.NE.1) THEN
!     EXCHANGE LINES OF MATRIX A AND ORDER OF UNKNOWNS
!     LINKED TO B
        DO 4 J=1,N
        T=A(1,J)
        A(1,J)=A(P,J)
        A(P,J)=T
      4 CONTINUE
        T=B(1)
        B(1)=B(P)
        B(P)=T
        ENDIF
!     FIRST LINE OF U AND FIRST COLUMN OF L
        W(1,1)=A(1,1)
        DO 6 J=2,N
        W(1,J)=A(1,J)
        W(J,1)=A(J,1)/W(1,1)
      6 CONTINUE
!     SECOND LINE OF U AND SECOND COLUMN OF L
!     (N-1)TH LINE OF U AND (N-1)TH COLUMN OF L
        DO 20 I=2,N-1
        DO 10 J=I,N
!     SEEK PIVOT
        AMAX=0.D0
        SUM=0.D0
        DO 8 K=1,I-1
        SUM=SUM+W(J,K)*W(K,I)
      8 CONTINUE
        T=ABS(A(J,I)-SUM)
        IF(T.LT.AMAX) GO TO 10
        AMAX=T
        P=J
     10 CONTINUE
        IF(AMAX.EQ.0.D0) GO TO 32
        IF(P.NE.I) THEN
!     EXCHANGE LINES OF MATRICES A , W AND ORDER OF UNKNOWNS
!     LINKED TO B
        DO 12 J=1,N
        T=A(I,J)
        A(I,J)=A(P,J)
        A(P,J)=T
        T=W(I,J)
        W(I,J)=W(P,J)
        W(P,J)=T
     12 CONTINUE
        T=B(I)
        B(I)=B(P)
        B(P)=T
        ENDIF
!     CALCULATE THE U(I,I)
        SUM=0.D0
        DO 14 K=1,I-1
        SUM=SUM+W(I,K)*W(K,I)
     14 CONTINUE
        W(I,I)=A(I,I)-SUM
        IF(W(I,I).EQ.0.D0) THEN
        GO TO 32
        ENDIF
        DO 20 J=I+1,N
        SUM=0.D0
        DO 16 K=1,I-1
        SUM=SUM+W(I,K)*W(K,J)
     16 CONTINUE
!     CALCULATE THE U(I,J)
        W(I,J)=A(I,J)-SUM
        SUM=0.D0
        DO 18 K=1,I-1
        SUM=SUM+W(J,K)*W(K,I)
     18 CONTINUE
!     CALCULATE THE L(I,J)
        W(J,I)=(A(J,I)-SUM)/W(I,I)
     20 CONTINUE
        SUM=0.D0
        DO 22 K=1,N-1
        SUM=SUM+W(N,K)*W(K,N)
     22 CONTINUE
!     CALCULATE  U(N,N)
        W(N,N)=A(N,N)-SUM
        IF(W(N,N).EQ.0.D0) GO TO 32
!     SOLVE SYSTEM  Z = L(-1)*B
        Z(1)=B(1)
        DO 26 I=2,N
        SUM=0.D0
        DO 24 K=1,I-1
        SUM=SUM+W(I,K)*Z(K)
     24 CONTINUE
        Z(I)=B(I)-SUM
     26 CONTINUE
!     SOLVE SYSTEM  X = U(-1)*Z
        X(N)=Z(N)/W(N,N)
        DO 30 I=N-1,1,-1
        SUM=0.D0
        DO 28 K=I+1,N
        SUM=SUM+W(I,K)*X(K)
     28 CONTINUE
        X(I)=(Z(I)-SUM)/W(I,I)
     30 CONTINUE
!     END OF RESOLUTION OF  A*X = L*(U*X) = B
        RETURN
     32 WRITE(*,'(1X,A)') '** DLITTL ** NO UNIQUE SOLUTION'
        RETURN
        END SUBROUTINE
     END MODULE

! End of file tdlittl.f90
