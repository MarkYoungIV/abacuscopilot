!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for writing results
!-----------------------------------------------------------------------------------!

MODULE results_mod
  USE kind_mod , ONLY : q1,q2
  USE matrix_mod
  USE options_mod 
  USE ions_mod
  USE charge_mod
  USE chgcar_mod
  USE cube_mod
  USE bader_mod
  USE voronoi_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: output

  CONTAINS 
!------------------------------------------------------------------------------------!
! output: Write out a summary of the bader analysis.
!         AtomVolumes.dat: Stores the 'significant' Bader volumes associated with
!                          each atom.
!         ACF.dat        : Stores the main output to the screen.
!         BCF.dat        : Stores 'significant' Bader volumes, their coordinates and
!                          charge, atom associated and distance to it. 
!------------------------------------------------------------------------------------!

  SUBROUTINE output(bdr,vor,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(voronoi_obj) :: vor
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg


    REAL(q2) :: sum_ionchg
    INTEGER :: i,bdimsig,mib,mab,cc,j,nmax
    INTEGER,DIMENSION(bdr%nvols) :: rck

    mab=MAXVAL(bdr%nnion)
    mib=MINVAL(bdr%nnion)
    OPEN(100,FILE='AtomVolumes.dat',STATUS='replace',ACTION='write')
    WRITE(100,'(A)') '   Atom                     Volume(s)   '
    WRITE(100,'(A,A)') '-----------------------------------------------------------',&
  &                    '-------------'

    bdr%tol=1.0e-4_q2

    DO i=mib,mab
      cc=0
      rck=0
      nmax=0
      DO j=1,bdr%nvols
        IF (bdr%volchg(j) > bdr%tol) THEN
          nmax=nmax+1
          IF(bdr%nnion(j) == i) THEN
            cc=cc+1
            rck(cc)=nmax
          END IF
        END IF
      END DO
      IF (cc == 0) CYCLE
      WRITE(100,'(2X,1I4,2X,A,2X,10000I5)') i,' ... ',rck(1:cc)
    END DO
    CLOSE(100)

    OPEN(100,FILE='ACF.dat',STATUS='replace',ACTION='write')
    WRITE(*,555) '#','X','Y','Z','VORONOI','BADER','%','MIN DIST'
    WRITE(100,555) '#','X','Y','Z','VORONOI','BADER','%','MIN DIST'
    555 FORMAT(/,4X,1A,9X,1A1,2(11X,1A1),8X,1A7,5X,1A5,9X,1A1,6X,1A10)
    WRITE(*,666) '----------------------------------------------------------------', &
  &              '------------------------------'
    WRITE(100,667) '---------------------------------------------------------------',&
  &              '------------------------------'
    666 FORMAT(1A66,1A26)
    667 FORMAT(1A65,1A27)


    sum_ionchg=SUM(bdr%ionchg)
!    write(*,*) sum_ionchg
    DO i=1,ions%nions

!     write(*,*) i
!     write(*,*) ions%r_car(i,:)
!     write(*,*) bdr%ionchg(i)
!     write(*,*) 100.0_q2*bdr%ionchg(i)/sum_ionchg
!     write(*,*) bdr%minsurfdist(i)


      WRITE(*,'(1I5,7F12.4)') i,ions%r_car(i,:),vor%vorchg(i),bdr%ionchg(i),         &
  &                           100.*bdr%ionchg(i)/sum_ionchg,bdr%minsurfdist(i)
      WRITE(100,'(1I5,7F12.4)') i,ions%r_car(i,:),vor%vorchg(i),bdr%ionchg(i),       &
  &                           100.*bdr%ionchg(i)/sum_ionchg,bdr%minsurfdist(i)
    END DO
    CLOSE(100)

    bdimsig=0
    OPEN(200,FILE='BCF.dat',STATUS='replace',ACTION='write')

    WRITE(200,556) '#','X','Y','Z','CHARGE','NEAREST ATOM','DISTANCE'
    556 FORMAT(/,4X,1A1,13X,1A1,2(16X,1A1),13X,1A7,2X,1A14,2X,1A10)

    WRITE(200,668) '---------------------------------------------------------------',&
  &                '----------------------------------------'
    668 FORMAT(1A65,1A37)

    bdr%tol=1.0e-4_q2

    DO i=1,bdr%nvols
        IF(bdr%volchg(i) > bdr%tol) THEN
           bdimsig=bdimsig+1
           WRITE(200,777) bdimsig,bdr%volchg(i),bdr%nnion(i),bdr%iondist(i)
           777 FORMAT(1I5,4(5X,F12.4),5X,1I5,5X,1F12.4)
        END IF
    END DO
    CLOSE(200)

!    OPEN(300,FILE='dipole.dat',STATUS='replace',ACTION='write')
!    WRITE(300,557) '#','X','Y','Z','MAGNITUDE'
!    557 FORMAT(/,4X,1A1,10X,1A1,2(15X,1A1),10X,1A10)
!    WRITE(300,*) '--------------------------------------------------------------------'
!    DO i=1,ndim
!      WRITE(300,888) i,dipole(i,:)*4.803_q2,                                         &
!  &                  sqrt(DOT_PRODUCT(dipole(i,:),dipole(i,:)))*4.803_q2
!!      888 FORMAT(1I5,4ES16.5)
!      888 FORMAT(1I5,4F16.6)
!    END DO
!    CLOSE(300)

    WRITE(*,'(/,2x,A,6X,1I8)')     'NUMBER OF BADER MAXIMA FOUND: ',bdr%nvols
    WRITE(*,'(2x,A,6X,1I8)')       '    SIGNIFICANT MAXIMA FOUND: ',bdimsig
    WRITE(*,'(2x,A,2X,1F12.5)')  '         NUMBER OF ELECTRONS: ',                 &
  &                                          SUM(bdr%volchg(1:bdr%nvols))

  RETURN
  END SUBROUTINE output

!------------------------------------------------------------------------------------!

END MODULE results_mod
