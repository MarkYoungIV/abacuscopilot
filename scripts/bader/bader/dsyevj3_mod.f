      MODULE dsyevj3_mod

        CONTAINS 
* ----------------------------------------------------------------------------
* Numerical diagonalization of 3x3 matrcies
* Copyright (C) 2006  Joachim Kopp
* ----------------------------------------------------------------------------
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
* ----------------------------------------------------------------------------


* ----------------------------------------------------------------------------
        SUBROUTINE DSYEVJ3(A, Q, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
* matrix A using the Jacobi algorithm.
* The upper triangular part of A is destroyed during the calculation,
* the diagonal elements are read but not destroyed, and the lower
* triangular elements are not referenced at all.
* ----------------------------------------------------------------------------
* Parameters:
*   A: The symmetric input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
*     .. Arguments ..
        DOUBLE PRECISION A(3,3)
        DOUBLE PRECISION Q(3,3)
        DOUBLE PRECISION W(3)

*     .. Parameters ..
        INTEGER          N
        PARAMETER        ( N = 3 )
    
*     .. Local Variables ..
        DOUBLE PRECISION SD, SO
        DOUBLE PRECISION S, C, T
        DOUBLE PRECISION G, H, Z, THETA
        DOUBLE PRECISION THRESH
        INTEGER          I, X, Y, R

*     Initialize Q to the identitity matrix
*     --- This loop can be omitted if only the eigenvalues are desired ---
        DO 10 X = 1, N
          Q(X,X) = 1.0D0
          DO 11, Y = 1, X-1
            Q(X, Y) = 0.0D0
            Q(Y, X) = 0.0D0
   11     CONTINUE
   10   CONTINUE

*       Initialize W to diag(A)
        DO 20 X = 1, N
          W(X) = A(X, X)
   20   CONTINUE

*       Calculate SQR(tr(A))  
        SD = 0.0D0
        DO 30 X = 1, N
          SD = SD + ABS(W(X))
   30   CONTINUE
        SD = SD**2
 
*       Main iteration loop
        DO 40 I = 1, 50
*         Test for convergence
          SO = 0.0D0
          DO 50 X = 1, N
            DO 51 Y = X+1, N
              SO = SO + ABS(A(X, Y))
   51       CONTINUE
   50     CONTINUE
          IF (SO .EQ. 0.0D0) THEN
            RETURN
          END IF

          IF (I .LT. 4) THEN
            THRESH = 0.2D0 * SO / N**2
          ELSE
            THRESH = 0.0D0
          END IF

*         Do sweep
          DO 60 X = 1, N
            DO 61 Y = X+1, N
              G = 100.0D0 * ( ABS(A(X, Y)) )
              IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X))
     $                    .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
                A(X, Y) = 0.0D0
              ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
*               Calculate Jacobi transformation
                H = W(Y) - W(X)
                IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                  T = A(X, Y) / H
                ELSE
                  THETA = 0.5D0 * H / A(X, Y)
                  IF (THETA .LT. 0.0D0) THEN
                    T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                  ELSE
                    T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                  END IF
                END IF

                C = 1.0D0 / SQRT( 1.0D0 + T**2 )
                S = T * C
                Z = T * A(X, Y)
                
*               Apply Jacobi transformation
                A(X, Y) = 0.0D0
                W(X)    = W(X) - Z
                W(Y)    = W(Y) + Z
                DO 70 R = 1, X-1
                  T       = A(R, X)
                  A(R, X) = C * T - S * A(R, Y)
                  A(R, Y) = S * T + C * A(R, Y)
   70           CONTINUE
                DO 80, R = X+1, Y-1
                  T       = A(X, R)
                  A(X, R) = C * T - S * A(R, Y)
                  A(R, Y) = S * T + C * A(R, Y)
   80           CONTINUE
                DO 90, R = Y+1, N
                  T       = A(X, R)
                  A(X, R) = C * T - S * A(Y, R)
                  A(Y, R) = S * T + C * A(Y, R)
   90           CONTINUE

*               Update eigenvectors
*               --- This loop can be omitted if only the eigenvalues are desired ---
                DO 100, R = 1, N
                  T       = Q(R, X)
                  Q(R, X) = C * T - S * Q(R, Y)
                  Q(R, Y) = S * T + C * Q(R, Y)
  100           CONTINUE
              END IF
   61       CONTINUE
   60     CONTINUE
   40   CONTINUE

        PRINT *, "DSYEVJ3: No convergence."
              
        END SUBROUTINE
* End of subroutine DSYEVJ3
      END MODULE
