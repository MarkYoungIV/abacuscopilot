!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for reading and writing charge density data
!-----------------------------------------------------------------------------------!

MODULE cube_mod
  USE kind_mod
  USE matrix_mod
  USE ions_mod
  USE charge_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_charge_cube, write_charge_cube

  CONTAINS

!-----------------------------------------------------------------------------------!
! read_charge_cube: Reads the charge density from a file in vasp or Gaussian cube 
!   format, by first reading the header, and then charges
!-----------------------------------------------------------------------------------!

  SUBROUTINE read_charge_cube(ions,chg,chargefile)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    CHARACTER(LEN=128) :: chargefile

    REAL(q2) :: vol
    REAL(q2), DIMENSION(3) :: dlat, dcar
    INTEGER :: i, n1, n2, n3, d1, d2, d3

    OPEN(100,FILE=chargefile(1:LEN_TRIM(ADJUSTL(chargefile))), &
    &    STATUS='old',ACTION='read')
    WRITE(*,'(/,1A11,1A20)') 'OPEN ... ',chargefile
    WRITE(*,'(1A27)') 'GAUSSIAN-STYLE INPUT FILE'
    ! Skip the first two lines
    READ(100,'(/)') 
    READ(100,*) ions%nions, chg%org_car

    ALLOCATE(ions%r_car(ions%nions,3))
    ALLOCATE(ions%r_dir(ions%nions,3))
    ALLOCATE(ions%ion_chg(ions%nions))
    ALLOCATE(ions%atomic_num(ions%nions))

    DO i=1,3
      READ(100,*) chg%npts(i),ions%lattice(i,1:3)
    END DO
    chg%lat2car = TRANSPOSE(ions%lattice)

    ! This should really indicate the units (Bohr/Ang)
    IF(chg%npts(1)<0) chg%npts(1) = (-1)*chg%npts(1)
    chg%i_npts = 1._q2/REAL(chg%npts,q2)
    ! The -1 is to account for having points at the edge of the cube

    DO i=1,3
      ions%lattice(i,:) = ions%lattice(i,:)*REAL(chg%npts(i),q2)
    END DO
    ! GH: still not sure about this -1.  Should it factor into the
    !     conversion matricies between lat and dir/car, as it is now?
    !     What about the places in which we multiply by npts, or i_npts?

    ions%dir2car = TRANSPOSE(ions%lattice)
    ions%car2dir = inverse(ions%dir2car)
    chg%car2lat = inverse(chg%lat2car)

    vol = matrix_volume(ions%lattice)
    DO i = 1, ions%nions
      READ(100,*) ions%atomic_num(i),ions%ion_chg(i),ions%r_car(i,:)
      ions%r_dir(i,:) = MATMUL(ions%car2dir,ions%r_car(i,:)-chg%org_car(:))
    END DO

    ! Note: this is only for isolated atoms.  For periodic systems, this shift
    ! might not be appropriate

    ! origin of the lattice is at chg(0.5,0.5,0.5)
!    chg%org_lat=(/0.5_q2,0.5_q2,0.5_q2/)
    chg%org_lat = (/1._q2,1._q2,1._q2/)
!    CALL matrix_vector(ions%car2dir,chg%org_car,chg%org_dir)

    ALLOCATE(ions%r_lat(ions%nions,3))
    DO i=1,ions%nions
!      ions%r_lat(i,:)=dir2lat(chg,ions%r_dir(i,:))
      ions%r_lat(i,:) = MATMUL(chg%car2lat, ions%r_car(i,:) - chg%org_car(:))
      ions%r_lat(i,:) = ions%r_lat(i,:) + chg%org_lat
      CALL pbc_r_lat(ions%r_lat(i,:),chg%npts)
    END DO

    chg%nrho = PRODUCT(chg%npts(:))
    ALLOCATE(chg%rho(chg%npts(1),chg%npts(2),chg%npts(3)))
    READ(100,*) (((chg%rho(n1,n2,n3),  &
    &            n3 = 1,chg%npts(3)), n2 = 1,chg%npts(2)),n1 = 1,chg%npts(1))
! GH: for some reason this is not working; replace with loop
!    chg%rho=chg%rho*vol
    DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
          chg%rho(n1,n2,n3) = vol*chg%rho(n1,n2,n3)
        END DO
      END DO
    END DO

    WRITE(*,'(2X,1A13,1I5,1A2,1I4,1A2,1I4)') 'DENSITY-GRID:',  &
    &         chg%npts(1),'x',chg%npts(2),'x',chg%npts(3)
    WRITE(*,'(2x,A,1A20)') 'CLOSE ... ', chargefile
    CLOSE(100)

    ! distance between neighboring points
    DO d1 = -1, 1
      dlat(1) = REAL(d1,q2)
      DO d2 = -1, 1
        dlat(2) = REAL(d2,q2)
        DO d3 = -1, 1
          dlat(3) = REAL(d3,q2)
          dcar = MATMUL(chg%lat2car, dlat)
          chg%lat_dist(d1,d2,d3) = SQRT(SUM(dcar*dcar))
         IF ((d1 == 0).AND.(d2 == 0).AND.(d3 == 0)) THEN
            chg%lat_i_dist(d1,d2,d3) = 0._q2
          ELSE
            chg%lat_i_dist(d1,d2,d3) = 1._q2/chg%lat_dist(d1,d2,d3)
          END IF
        END DO
      END DO
    END DO
 
    !add ions%num_ion as the total atom number in case some want to write 
    !a chgcar file from a cube file
!    ALLOCATE(ions%num_ion(1))
!    ions%num_ion(1)=ions%nions

  RETURN
  END SUBROUTINE read_charge_cube

!-----------------------------------------------------------------------------------!
! write_charge_cube: Write out a Gaussian cube type file
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_charge_cube(ions,chg,chargefile)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    CHARACTER(LEN=128) :: chargefile
    CHARACTER(LEN=128) :: cubefile
    
    INTEGER :: i,n1,n2,n3
    REAL(q2),DIMENSION(3,3) :: voxel
    REAL(q2) :: vol

!   EM: change the file extension to something more descriptive   
    cubefile = chargefile(1:(index(chargefile,'.'))) // 'cube'

    OPEN(100,FILE = cubefile(1:LEN_TRIM(ADJUSTL(cubefile))), STATUS='replace')

    WRITE(100,*) 'Gaussian cube file'
    WRITE(100,*) 'Bader charge'

!   EM: Formats changed to comply with G03 cube standard
!    WRITE(100,'(1I5,3(3X,1F9.6))') ions%nions,chg%org_car
    WRITE(100,'(I5,3(F12.6))') ions%nions, chg%org_car

    DO i = 1, 3
      voxel(i,:) = ions%lattice(i,:)/REAL(chg%npts(i),q2)
      WRITE(100,'(1I5,3(3X,1F9.6))') chg%npts(i), voxel(i,:)
    END DO

    DO i = 1, ions%nions
      WRITE(100,'(I5,F12.6,3(F12.6))') & 
      &     ions%atomic_num(i), ions%ion_chg(i), ions%r_car(i,:)
    END DO

    vol = matrix_volume(ions%lattice)
    DO n1 = 1, chg%npts(1)
      DO n2 = 1, chg%npts(2)
!   EM: Voxel data is printed with newlines after inner z-loop
!       (see ex. http://local.wasp.uwa.edu.au/~pbourke/dataformats/cube/
!        DO n3=1,chg%npts(3)
!          chg%rho(n1,n2,n3)=chg%rho(n1,n2,n3)/vol
!        END DO
        WRITE(100,'(6E13.5)') ((chg%rho(n1,n2,n3))/vol, n3 = 1,chg%npts(3))
      END DO  
    END DO

    CLOSE(100)
 
  RETURN
  END SUBROUTINE write_charge_cube

END MODULE cube_mod

