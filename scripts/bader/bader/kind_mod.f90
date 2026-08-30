!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!   Module specifying precision and constants
!-----------------------------------------------------------------------------------!
MODULE kind_mod
  IMPLICIT NONE

  PUBLIC

! Public parameters
  INTEGER,PARAMETER :: q1=SELECTED_REAL_KIND(6,30)         ! single precision
  INTEGER,PARAMETER :: q2=SELECTED_REAL_KIND(15,305)       ! double precision
  REAL(q2),PARAMETER :: pi=3.141592653589793238462643_q2    

!  INTEGER,PARAMETER :: QSORT_THRESHOLD = 16

END MODULE kind_mod

