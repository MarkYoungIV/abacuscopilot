!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for reading and writing VASP CHGCAR files
!-----------------------------------------------------------------------------------!

MODULE chgcar_mod
  USE kind_mod
  USE matrix_mod
  USE ions_mod 
  USE charge_mod 
  USE options_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_charge_chgcar, write_charge_chgcar

  CONTAINS

!-----------------------------------------------------------------------------------!
! read_charge_chgcar: Reads the charge density from a file in vasp format
!-----------------------------------------------------------------------------------!

  SUBROUTINE read_charge_chgcar(ions,chg,chargefile,opts)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    CHARACTER(LEN=128) :: chargefile
    CHARACTER(LEN=128) :: line
    TYPE(options_obj) :: opts

    REAL(q2) :: scalefactor
    REAL(q2), DIMENSION(3) :: dlat,dcar
    INTEGER :: i,n1,n2,n3,d1,d2,d3
    INTEGER, DIMENSION(110) :: nionlist=0

    OPEN(100,FILE=chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
    WRITE(*,'(/,1A11,1A20)') 'OPEN ... ',chargefile
    READ(100,'(/,1F20.16)') scalefactor
    READ(100,'(3F13.6)') (ions%lattice(i,1:3) , i=1,3)

    IF(opts%in_opt == opts%in_chgcar5) THEN
      READ(100,'(1X,330A)') ions%name_ion
      WRITE(*,'(2x,A)') 'VASP5-STYLE INPUT FILE'
    ELSE
      WRITE(*,'(2x,A)') 'VASP4-STYLE INPUT FILE'
    END IF

!GH: vasp switched to I6 from I4 at 5.2.11; to be compatible with both
!    use free format.  But, will fail for >999 atoms with vasp4.x

!ORIGINAL
!    READ(100,'(110I4)') nionlist
!    READ(100,*)

!NEXT: this version fails for gfortran on linux
!    READ(100,*,ERR=999) nionlist
!    WRITE(*,*) "nionlist: ",nionlist
!999 CONTINUE

!WORKING:
    READ(100,'(A)') line
    READ(line,*,END=999) nionlist
999 CONTINUE
    DO i = 1, 110
     IF(nionlist(i).eq.0) EXIT
    ENDDO
    ions%niontypes = i - 1
    ALLOCATE(ions%num_ion(ions%niontypes))
    DO i = 1, ions%niontypes
      ions%num_ion(i) = nionlist(i)
    END DO
    ions%scalefactor = scalefactor
    ions%nions = SUM(ions%num_ion)
    ions%lattice = scalefactor*ions%lattice
    ions%dir2car = TRANSPOSE(ions%lattice)
    ions%car2dir = inverse(ions%dir2car)
    ALLOCATE(ions%r_dir(ions%nions,3))
    ALLOCATE(ions%r_car(ions%nions,3))
    READ(100,*) !GH: read Direct line
    DO i = 1, ions%nions
      READ(100,'(3F10.6)') ions%r_dir(i,:)
      ions%r_car(i,:) = MATMUL(ions%dir2car, ions%r_dir(i,:))
    END DO
    READ(100,*) 
    READ(100,*) chg%npts(:)
    chg%nrho = PRODUCT(chg%npts(:))
    chg%i_npts = 1.0_q2/REAL(chg%npts,q2)
    ALLOCATE(chg%rho(chg%npts(1),chg%npts(2),chg%npts(3)))
    READ(100,*) (((chg%rho(n1,n2,n3), &
    &  n1 = 1,chg%npts(1)), n2 = 1,chg%npts(2)), n3 = 1,chg%npts(3))
    WRITE(*,'(2X,1A13,1I5,1A2,1I4,1A2,1I4)') &
    &  'DENSITY-GRID:',chg%npts(1),'x',chg%npts(2),'x',chg%npts(3)
    WRITE(*,'(2x,A,1A20)') 'CLOSE ... ', chargefile
    CLOSE(100)
    DO i=1,3
      chg%lat2car(:,i) = ions%dir2car(:,i)/REAL(chg%npts(i),q2)
    END DO
    chg%car2lat = inverse(chg%lat2car)

    ! origin of the lattice is at chg(1,1,1)
    chg%org_lat = (/1._q2,1._q2,1._q2/)
    chg%org_car = (/0._q2,0._q2,0._q2/)

    ! ion positions in grid points
    ALLOCATE(ions%r_lat(ions%nions,3))
    DO i = 1,ions%nions
      ions%r_lat(i,:) = MATMUL(chg%car2lat, ions%r_car(i,:))
      ions%r_lat(i,:) = ions%r_lat(i,:) + chg%org_lat
      CALL pbc_r_lat(ions%r_lat(i,:), chg%npts)
    END DO

    ! distance between neighboring points
    DO d1 = -1, 1
      dlat(1) = REAL(d1,q2)
      DO d2 = -1, 1
        dlat(2) = REAL(d2,q2)
        DO d3 = -1, 1
          dlat(3) = REAL(d3,q2)
          dcar = MATMUL(chg%lat2car, dlat)
          chg%lat_dist(d1, d2, d3) = SQRT(SUM(dcar*dcar))
          IF ((d1 == 0).AND.(d2 == 0).AND.(d3 == 0)) THEN
            chg%lat_i_dist(d1, d2, d3) = 0._q2
          ELSE
            chg%lat_i_dist(d1, d2, d3) = 1._q2/chg%lat_dist(d1, d2, d3)
          END IF
        END DO
      END DO
    END DO
  
  RETURN
  END SUBROUTINE read_charge_chgcar

!-----------------------------------------------------------------------------------!
! write_charge_chgcar: Write the charge density from a file in vasp format
!-----------------------------------------------------------------------------------!
    
  SUBROUTINE write_charge_chgcar(ions,chg,chargefile,opts)
    
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    CHARACTER(LEN=128) :: chargefile
    TYPE(options_obj) :: opts

    INTEGER :: i, n1, n2, n3, it

    OPEN(100, FILE=chargefile(1:LEN_TRIM(ADJUSTL(chargefile))), STATUS='replace')
    WRITE(100,*)'Bader charge density file'
    WRITE(100,*)'1.00'
    WRITE(100,'(3F13.6)') (ions%lattice(i,1:3) , i=1,3)
    IF(opts%in_opt == opts%in_chgcar5) THEN
      ions%name_ion = ADJUSTL(ions%name_ion)
      it = LEN_TRIM(ions%name_ion)
      WRITE(100,'(3X,A)') ions%name_ion(1:it)
    END IF
    WRITE(100,'(110I4)') ions%num_ion
    WRITE(100,*)'Direct'
    WRITE(100,'(3(2X,1F8.6))') (ions%r_dir(i,:), i = 1, ions%nions)
    WRITE(100,*)
    WRITE(100,*) chg%npts(1), chg%npts(2), chg%npts(3)
    WRITE(100,'(5(1X,E17.11))') (((chg%rho(n1,n2,n3), &
    &  n1=1,chg%npts(1)), n2=1,chg%npts(2)), n3=1,chg%npts(3))
    CLOSE(100)

  RETURN
  END SUBROUTINE write_charge_chgcar

!-----------------------------------------------------------------------------------!

END MODULE chgcar_mod
