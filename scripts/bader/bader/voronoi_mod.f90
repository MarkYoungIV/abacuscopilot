!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge with a voronoi analysis
!-----------------------------------------------------------------------------------!
MODULE voronoi_mod
  USE kind_mod
  USE matrix_mod
  USE ions_mod
  USE charge_mod
  USE options_mod
  USE io_mod
  IMPLICIT NONE

! Public, allocatable variables
  TYPE voronoi_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: vorchg
  END TYPE

  PRIVATE
  PUBLIC :: voronoi_obj
  PUBLIC :: voronoi

  CONTAINS

!-----------------------------------------------------------------------------------!
!  voronoi:  Calculate the charge density populations based upon Voronoi polyhedra.
!    In this scheme each element of charge density is associated with the atom that
!    it is closest to.
!-----------------------------------------------------------------------------------!

  SUBROUTINE voronoi(vor,ions,chg)

    TYPE(voronoi_obj) :: vor
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2),DIMENSION(3) :: r_lat, dr_lat, dr_car
    REAL(q2) :: dist, min_dist, vol
    INTEGER :: i, n1, n2, n3, closest, tenths_done, cr, count_max, t1, t2
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: ionvol

    CALL SYSTEM_CLOCK(t1, cr, count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING VORONOI CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    ALLOCATE(vor%vorchg(ions%nions))

    vor%vorchg = 0._q2
    tenths_done = 0

    ALLOCATE(ionvol(ions%nions))
    ionvol = 0.0

    DO n1 = 1, chg%npts(1)
      r_lat(1) = REAL(n1,q2)
      IF ((n1*10/chg%npts(1)) > tenths_done) THEN
        tenths_done = (n1*10/chg%npts(1))
        WRITE(*,'(A,$)') '**'
      END IF
      DO n2 = 1, chg%npts(2)
        r_lat(2) = REAL(n2,q2)
        DO n3 = 1, chg%npts(3)
          r_lat(3) = REAL(n3,q2)
          closest = 1
!          r_car = lat2car(chg, r_lat)
!          dr_car = r_car - ions%r_car(1,:)
          dr_lat = r_lat - ions%r_lat(1,:)
          dr_car = MATMUL(chg%lat2car, dr_lat)
          CALL dpbc_car(ions, dr_car)
          min_dist = DOT_PRODUCT(dr_car,dr_car)
          DO i=2,ions%nions
!            dr_car = r_car - ions%r_car(i,:)
            dr_lat = r_lat - ions%r_lat(i,:)
            dr_car = MATMUL(chg%lat2car, dr_lat)
            CALL dpbc_car(ions, dr_car)
            dist = DOT_PRODUCT(dr_car, dr_car)
            IF (dist < min_dist) THEN
              min_dist = dist
              closest = i
            END IF
          END DO
          ionvol(closest) = ionvol(closest) + 1._q2
          vor%vorchg(closest) = vor%vorchg(closest) + rho_val(chg,n1,n2,n3)
        END DO
      END DO
    END DO

    vor%vorchg(:) = vor%vorchg(:)/REAL(chg%nrho,q2)

    vol = matrix_volume(ions%lattice)
    vol = vol/chg%nrho
    DO i = 1, ions%nions
      ionvol(i) = ionvol(i)*vol
    END DO

    CALL SYSTEM_CLOCK(t2, cr, count_max)
    WRITE(*,'(A12,F7.2,A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

    WRITE(*,'(/,2X,A,/)') 'VORONOI ANALYSIS RESULT'
    WRITE(*,556) '#','X','Y','Z','CHARGE','ATOMIC VOL'
    556 FORMAT(4X,1A1,9X,1A1,2(11X,1A1),8X,1A6,6X,1A10)

    WRITE(*,'(A)') '  ----------------------------------------------------------------------'
    DO i=1,ions%nions
       WRITE(*,777) i,ions%r_car(i,:),vor%vorchg(i),ionvol(i)
       777 FORMAT(1I5,4F12.4,3X,1F12.4)
    END DO
    WRITE(*,'(A)') '  -----------------------------------------------------------------------'

    WRITE(*,'(2x,A,2X,1F12.5)')  '         NUMBER OF ELECTRONS: ', &
    &                                        SUM(vor%vorchg(1:ions%nions))

    DEALLOCATE(ionvol)
    
  RETURN
  END SUBROUTINE voronoi

!-----------------------------------------------------------------------------------!

END MODULE voronoi_mod
