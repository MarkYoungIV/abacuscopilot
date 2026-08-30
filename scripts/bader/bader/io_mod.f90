!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for reading and writing charge density data
!-----------------------------------------------------------------------------------!

MODULE io_mod
  USE kind_mod
  USE matrix_mod
  USE options_mod 
  USE ions_mod
  USE charge_mod
  USE chgcar_mod
  USE cube_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_charge, read_charge_ref, write_charge

  CONTAINS

!-----------------------------------------------------------------------------------!
! read_charge: Reads the charge density from a file in vasp or Gaussian cube format,
!    by first reading the header, and then charges
!-----------------------------------------------------------------------------------!

  SUBROUTINE read_charge(ions,chg,opts)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts

    CHARACTER(LEN=128) :: chargefile
    CHARACTER(LEN=7) :: text1,text2
    INTEGER :: cr,count_max,t1,t2,it1,it2

    CALL SYSTEM_CLOCK(t1,cr,count_max)

    chargefile=opts%chargefile
    IF (opts%in_opt == opts%in_auto .OR. opts%in_opt==opts%in_chgcar) THEN
      ! Try to guess the file type
      OPEN(100,FILE=chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
      READ(100,'(6/,1A7)') text1
      READ(100,'(1A7)') text2
!      WRITE(*,*) text1; WRITE(*,*) text2
      CLOSE(100)
      text1=ADJUSTL(text1); text2=ADJUSTL(text2)
      it1=LEN_TRIM(text1); it2=LEN_TRIM(text2)
      IF (text1(1:it1) == 'Direct') THEN
        opts%in_opt=opts%in_chgcar4
      ELSE IF (text2(1:it2) == 'Direct') THEN
        opts%in_opt=opts%in_chgcar5
      ELSE 
        opts%in_opt=opts%in_cube
      ENDIF
    ENDIF

    IF (opts%in_opt == opts%in_chgcar4) THEN
      ! make the output in the same format as input
      opts%out_opt=opts%out_chgcar4
      CALL read_charge_chgcar(ions,chg,chargefile,opts)
    ELSEIF  (opts%in_opt == opts%in_chgcar5) THEN
      opts%out_opt=opts%out_chgcar5
      CALL read_charge_chgcar(ions,chg,chargefile,opts)
    ELSEIF (opts%in_opt == opts%in_cube) THEN
      opts%out_opt=opts%out_cube
      CALL read_charge_cube(ions,chg,chargefile)
    ENDIF

    CALL SYSTEM_CLOCK(t2,cr,count_max)
    WRITE(*,'(1A12,1F7.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'
   
  RETURN
  END SUBROUTINE read_charge

!-----------------------------------------------------------------------------------!
! read_charge_ref: Reads the reference charge density from a file in vasp or
!    Gaussian cube format, by first reading the header, and then the charge density
!-----------------------------------------------------------------------------------!

  SUBROUTINE read_charge_ref(ions,chg,opts)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts

    CHARACTER(LEN=128) :: chargefile
    INTEGER :: cr,count_max,t1,t2

    CALL SYSTEM_CLOCK(t1,cr,count_max)

    chargefile=opts%refchgfile
    IF (opts%in_opt == opts%in_chgcar4 .OR. opts%in_opt == opts%in_chgcar5) THEN
      CALL read_charge_chgcar(ions,chg,chargefile,opts)
    ELSEIF (opts%in_opt == opts%in_cube) THEN
      CALL read_charge_cube(ions,chg,chargefile)
    ENDIF

    CALL SYSTEM_CLOCK(t2,cr,count_max)
    WRITE(*,'(1A12,1F7.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE read_charge_ref

!-----------------------------------------------------------------------------------!
! write_charge: Writes the charge density to a file in vasp or Gaussian cube format
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_charge(ions,chg,opts,chargefile)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts

    CHARACTER(LEN=128) :: chargefile

    IF ( opts%out_opt == opts%out_auto ) THEN
      ! write output file as the same type as input
      opts%out_opt=opts%in_opt
    ENDIF
    IF (opts%out_opt == opts%out_chgcar4 .OR. opts%out_opt == opts%out_chgcar5) THEN
      CALL write_charge_chgcar(ions,chg,chargefile,opts)
    ELSEIF (opts%out_opt == opts%out_cube) THEN 
      CALL write_charge_cube(ions,chg,chargefile)
    ENDIF
    
  RETURN
  END SUBROUTINE write_charge

!-----------------------------------------------------------------------------------!

END MODULE io_mod
