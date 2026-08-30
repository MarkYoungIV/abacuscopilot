!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for reading wavefunctions
!-----------------------------------------------------------------------------------!

MODULE wave_mod
  USE kind_mod , ONLY : q2
  USE matrix_mod
  USE ions_mod , ONLY : ions_obj
  IMPLICIT NONE

! Public, allocatable variables

  TYPE :: wave_obj
    COMPLEX(q2), ALLOCATABLE,DIMENSION(:,:) :: ac,dc
    COMPLEX(q2), ALLOCATABLE,DIMENSION(:) :: coef
    REAL(q2),DIMENSION(3,3) :: A
    REAL(q2),DIMENSION(3) :: kpt

    COMPLEX(q2) :: eval
    REAL(q2) :: emax,overlap,fweight,gweight
    INTEGER :: nkpt,npw,nband,ic
    INTEGER :: i,j,k,iband,ikpt,nxn
    CHARACTER(LEN=5) :: code
  END TYPE


  PRIVATE
  PUBLIC :: wave_obj,read_wavecar

  CONTAINS

!---------------------------------------------------------------------
! Read in WAVECAR
!---------------------------------------------------------------------

  SUBROUTINE read_wavecar(wave)
    TYPE(wave_obj) :: wave

    open(unit=12,file="WAVECAR",status="old",form="unformatted")

    ! Read in from WAVECAR:
    ! the number of k-points
    ! the number of bands
    ! the energy maximum
    ! the cell dimensions
    ! the type of VASP code

    read(12) nkpt,nband,emax,((A(i,j),i=1,3),j=1,3)

    !????????????????????????????????????????????????????????
    !      write(19,*) 
    !      write(19,*) 'nkpt  =',nkpt
    !      write(19,*) 'nband =',nband
    !      write(19,*) 'emax  =',emax
    !      write(19,*) 'A='
    !      write(19,'(3X,3(1X,f8.3))') (A(i,1),i=1,3)
    !      write(19,'(3X,3(1X,f8.3))') (A(i,2),i=1,3)
    !      write(19,'(3X,3(1X,f8.3))') (A(i,3),i=1,3)
    !????????????????????????????????????????????????????????

    read(12) code
    DO ikpt=1, nkpt

      ! Read in from WAVECAR:
      ! the number of pw wave fns
      ! the k-point
      
      read(12) npw, kpt(1:3)

      !????????????????????????????????????????????????????????
      ! write(19,*)
      ! write(19,'("k-point #",I3,":  (",3f7.4,")    npw=",I6)') 
      !		& ikpt, (kpt(i),i=1,3),npw 
      ! write(19,*) 'kpt =',kpt
      ! write(19,*) 'npw =',npw
      ! write(19,*) "  band       energy        weight"
      !
      !????????????????????????????????????????????????????????

      ALLOCATE(coef(npw))
      ALLOCATE(ac(npw,200))
      ic = 1 
      DO iband = 1, nband
            
        ! Read in from WAVECAR:
        ! the energy (eval)
        ! the occupation (fweight)

        read(12) eval, fweight, (coef(i),i=1,npw)

        !????????????????????????????????????????????????????????
        ! write(19,'(5X,I3,5X,f8.4,5x,f8.4)') iband, dreal(eval), fweight
        !????????????????????????????????????????????????????????

        ! Create matrix (ac):
        ! row: coefficients for plane waves
        ! column: coeff of 1st plane wave for bands interested in 

!        IF(iband >= nbandmin.and.iband.le.nbandmax) then
        DO i = 1,npw
          ac(i,ic) = coef(i)
        ENDDO 
          ic = ic+1
!       end if

      ENDDO
      DEALLOCATE(coef)
    ENDDO
    close(unit=12)

! Read second WAVECAR into memory
!---------------------------------------------------------------------

    open(unit=12,file="WAVECARNEW",status="old",form="unformatted")

! Read in from WAVECAR:
! the number of k-points
! the number of bands
! the energy maximum
! the cell dimensions
! the type of VASP code

    read(12) nkpt,nband,emax,((A(i,j),i=1,3),j=1,3)

    !????????????????????????????????????????????????????????
    !      write(19,*)
    !      write(19,*) 'nkpt  =',nkpt
    !      write(19,*) 'nband =',nband
    !      write(19,*) 'emax  =',emax
    !      write(19,*) 'A='
    !      write(19,'(3X,3(1X,f8.3))') (A(i,1),i=1,3)
    !      write(19,'(3X,3(1X,f8.3))') (A(i,2),i=1,3)
    !      write(19,'(3X,3(1X,f8.3))') (A(i,3),i=1,3)
    !????????????????????????????????????????????????????????

    read(12) code

    do ikpt=1, nkpt

! Read in from WAVECAR:
! the number of pw wave fns?
! the k-point?

      read(12) npw, kpt(1:3)

      !????????????????????????????????????????????????????????
      ! write(19,*)
      !write(19,'("k-point #",I3,":  (",3f7.4,")    npw=",I6)')
      !               & ikpt, (kpt(i),i=1,3),npw
      ! write(19,*) 'kpt =',kpt
      ! write(19,*) 'npw =',npw
      ! write(19,*) "  band       energy        weight"
      !
      !????????????????????????????????????????????????????????

        allocate(coef(npw))
        allocate(dc(npw,200))
        ic = 1

        do iband = 1, nband

! Read in from WAVECAR:
! the energy (eval)
! the occupation

          read(12) eval, fweight, (coef(i),i=1,npw)

      !????????????????????????????????????????????????????????
      ! write(19,'(5X,I3,5X,f8.4,5x,f8.4)') iband, dreal(eval), fweight
      !????????????????????????????????????????????????????????

! Create matrix (dc):
! row: coefficients for plane waves
! column: coeff of 1st plane wave for bands interested in

          if(iband.ge.nbandmin.and.iband.le.nbandmax) then
            do i = 1,npw
              dc(i,ic) = coef(i)
            enddo
            ic = ic+1
          end if

        enddo
        deallocate(coef)

      enddo

      close(unit=12)

! THESE OVERLAPS SHOULD BE DIVIDED BY 2

!---------------------------------------------------------------------

! Overlap of "old wavecar" bands with other "old wavecar" bands
! matrix nxn by nxn:
! row1: band1 overlap with bands bandbegin through bandend

        write(13,*)"old wavecar"
        nxn=ndiff+1
        do i = 1,nxn
          do j = 1,nxn
            overlap = 0.
            do k = 1,npw
              overlap=overlap+conjg(ac(k,i))*ac(k,j)
            end do
            write(13,'(E25.18)', ADVANCE='NO') overlap
          end do
            write(13,*)" "
        end do
!--------------------------------------------------------

! Overlap of "new wavecar" bands with other "old wavecar" bands
! matrix nxn by nxn:
! row1: band1 overlap with bands bandbegin through bandend

    write(13,*)"new wavecar"
        nxn=ndiff+1
         do i = 1,nxn
          do j = 1,nxn
            overlap = 0.
            do k = 1,npw
              overlap=overlap+conjg(dc(k,i))*dc(k,j)
            end do
            write(13,'(E25.18)', ADVANCE='NO') overlap
          end do
            write(13,*)" "
        end do
!--------------------------------------------------------

! Overlap of "old wavecar" bands with other "old wavecar" bands
! matrix nxn by nxn:
! row1: band1 overlap with bands bandbegin through bandend

	    write(13,*)"mix wavecars"
	nxn=ndiff+1
	do i = 1,nxn
	  do j = 1,nxn
	    overlap = 0.
	    do k = 1,npw
	      overlap=overlap+conjg(ac(k,i))*dc(k,j)
	    end do
            write(13,'(E25.18)', ADVANCE='NO') overlap
	  end do
	    write(13,*)" "
	end do
!--------------------------------------------------------

! Calculate overlap we actually use (real### files)
! matrix nxn by nxn:

	    write(13,*)"first term"
	nxn=ndiff+1
	do i = 1,nxn
	  do j = 1,nxn
	    overlap = 0.
	    do k = 1,npw
	      overlap=overlap+conjg(ac(k,i))*dc(k,j) 
	      overlap=overlap-conjg(dc(k,i))*ac(k,j) 
	    end do
            write(13,'(E25.18)', ADVANCE='NO') overlap
	  end do
	    write(13,*)" "
	end do
!--------------------------------------------------------

! Yet more overlap
! matrix nxn by nxn:

        write(13,*)"second term"
    nxn=ndiff+1
    do i = 1,nxn
      do j = 1,nxn
        overlap = 0.
        do k = 1,npw
          overlap=overlap+conjg(dc(k,i))*dc(k,j) 
          overlap=overlap-conjg(ac(k,i))*ac(k,j) 
        end do
            write(13,'(E25.18)', ADVANCE='NO') overlap
      end do
        write(13,*)" "
    end do
!--------------------------------------------------------
    end

