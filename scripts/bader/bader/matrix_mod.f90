!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module containing matrix functions
!-----------------------------------------------------------------------------------!
MODULE matrix_mod
  USE kind_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: inverse, adjoint, matrix_volume, triple_product
  PUBLIC :: cross_product, determinant, eigenvectors, v2area
  PUBLIC :: NormalizeLatticeVectors
  CONTAINS

!-----------------------------------------------------------------------------------!
! inverse: return the inverse of A(3,3)
!-----------------------------------------------------------------------------------!

  FUNCTION inverse(A)

    REAL(q2), INTENT(IN), DIMENSION(3,3) :: A
    REAL(q2), DIMENSION(3,3) :: inverse
    REAL(q2) :: det

!    write(*,*) A(1,1),A(1,2),A(1,3)
!    write(*,*) A(2,1),A(2,2),A(2,3)
!    write(*,*) A(3,1),A(3,2),A(3,3)
    det = determinant(A)
    IF (det == 0) STOP 'Divide by zero in matrix inverse'
    inverse = adjoint(A)/det

  RETURN
  END FUNCTION inverse

!-----------------------------------------------------------------------------------!
! adjoint: return the adjoint of A(3,3)
!-----------------------------------------------------------------------------------!

  FUNCTION adjoint(A)

    REAL(q2), INTENT(IN), DIMENSION(3,3) :: A
    REAL(q2), DIMENSION(3,3) :: adjoint

    adjoint(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
    adjoint(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
    adjoint(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)

    adjoint(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
    adjoint(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
    adjoint(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)

    adjoint(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
    adjoint(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
    adjoint(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

  RETURN
  END FUNCTION adjoint

!-----------------------------------------------------------------------------------!
! matrix_volume: function returning the triple product of the lattice vectors
!-----------------------------------------------------------------------------------!

  FUNCTION matrix_volume(h)

    REAL(q2),INTENT(IN),DIMENSION(3,3) :: h
    REAL(q2) :: matrix_volume

    matrix_volume = h(1,1)*(h(2,2)*h(3,3) - h(2,3)*h(3,2))  &
    &             - h(1,2)*(h(2,1)*h(3,3) - h(3,1)*h(2,3))  &
    &             + h(1,3)*(h(2,1)*h(3,2) - h(3,1)*h(2,2))
    matrix_volume = ABS(matrix_volume)

  RETURN
  END FUNCTION matrix_volume
  
!-----------------------------------------------------------------------------------!
! triple_product 
!-----------------------------------------------------------------------------------!

  FUNCTION triple_product(a,b,c)

    REAL(q2), INTENT(IN), DIMENSION(:) :: a,b,c
    REAL(q2) :: triple_product

    triple_product = c(1)*(a(2)*b(3) - a(3)*b(2))  & 
                   + c(2)*(a(3)*b(1) - a(1)*b(3))  &
                   + c(3)*(a(1)*b(2) - a(2)*b(1))

    RETURN
  END FUNCTION
 
!-----------------------------------------------------------------------------------!
! cross_product: calculate the cross product of two vectors
!-----------------------------------------------------------------------------------!

  FUNCTION cross_product(A,B)

    REAL(q2), INTENT(IN), DIMENSION(3) :: A,B
    REAL(q2), DIMENSION(3) :: cross_product

    cross_product(1) = A(2)*B(3) - A(3)*B(2)
    cross_product(2) = A(3)*B(1) - A(1)*B(3)
    cross_product(3) = A(1)*B(2) - A(2)*B(1)

    RETURN
  END FUNCTION

!-----------------------------------------------------------------------------------!
! determinant: of a 3x3 matrix
!-----------------------------------------------------------------------------------!

  FUNCTION determinant(A)

    REAL(q2), INTENT(IN), DIMENSION(3,3) :: A
    REAL(q2) :: determinant

    determinant = A(1,1)*A(2,2)*A(3,3) &
                - A(1,1)*A(2,3)*A(3,2) &
                - A(1,2)*A(2,1)*A(3,3) &
                + A(1,2)*A(2,3)*A(3,1) &
                + A(1,3)*A(2,1)*A(3,2) &
                - A(1,3)*A(2,2)*A(3,1)

    RETURN
    END FUNCTION

!----------------------------------------------------------------------
! eigen_vectors : find the eigenvectors
! ***************** KING RAYCHARD :***********************************
! *************** THIS MAY NOT BE NEEDED *****************************
!----------------------------------------------------------------------

  SUBROUTINE eigenvectors(yita2,iDM,dM,s1,s2,v1,v2,v3)

    REAL(q2), INTENT(IN) :: yita2
    REAL(q2), INTENT(IN), DIMENSION(3) :: s1, s2, v1
    REAL(q2), INTENT(IN), DIMENSION(3,3) :: iDM, dM
    REAL(q2), INTENT(OUT), DIMENSION(3) :: v2, v3
    REAL(q2), DIMENSION(3) :: u1, u2, w1
    REAL(q2), DIMENSION(3,3) :: tempMat
    REAL(q2) :: norm

    tempMat = iDM*yita2
    tempMat = dM - tempMat
    u1 = MATMUL(tempMat,s1)
    u2 = MATMUL(tempMat,s2)
    norm = 1.0_q2/SQRT(SUM(u1(:)**2))
    w1 = u1*norm
    v2 = cross_product(w1,v1)
    v3 = cross_product(v1,v2)

    RETURN
  END SUBROUTINE

!-----------------------------------------------------------------------------------!
! v2_area : find the area of a parallelogram defined by two vectors
!-----------------------------------------------------------------------------------!

  FUNCTION v2area(v1,v2)

    REAL(q2), INTENT(IN), DIMENSION(3) :: v1,v2
    REAL(q2) :: v2area
    REAL(q2), DIMENSION(3) :: tempvec
    tempvec = cross_product(v1, v2)
    v2area = SQRT(SUM(tempvec(:)**2))

    RETURN
  END FUNCTION

!-----------------------------------------------------------------------------------!
! NormalizeLatticeVectors
! This function takes the lattice vectors and normalize each of them.
!-----------------------------------------------------------------------------------!  
  FUNCTION NormalizeLatticeVectors(v1)
    ! v1 should be ions%lattice
    REAL(q2),INTENT(IN),DIMENSION(3,3):: v1
    REAL(q2),DIMENSION(3,3) :: NormalizeLatticeVectors
    REAL(q2) :: normcoeff
    INTEGER :: i,j
    i = 1

    DO WHILE (i <= 3)
      normcoeff = 0
      j = 1
      DO WHILE (j <= 3)
        normcoeff = normcoeff + v1(i,j)**2 
        j = j + 1
      END DO
      j = 1
      DO WHILE (j <= 3)
        NormalizeLatticeVectors(i,j) = v1(i,j) / SQRT(normcoeff)
        j = j + 1
      END DO
      i = i + 1  
    END DO
    ! BECAUSE FORTRAN ROW/COLUMN ORDER MIX UP
    NormalizeLatticeVectors = TRANSPOSE(NormalizeLatticeVectors)
    RETURN 
  END FUNCTION

END MODULE matrix_mod

