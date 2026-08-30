  MODULE automatedcp_mod
    USE kind_mod
    USE matrix_mod
    USE bader_mod
    USE charge_mod 
    USE options_mod
    USE ions_mod
    USE io_mod
    USE ions_mod
    USE weight_mod
    USE dsyevj3_mod
    USE critpoint_mod
    IMPLICIT NONE
    CONTAINS 

    SUBROUTINE AutoCP(bdr,chg,opts,ions,stat)
      TYPE (bader_obj) :: bdr
      TYPE (charge_obj) :: chg
      TYPE (ions_obj)  :: ions
      TYPE (options_obj) :: opts
      INTEGER :: stat
      LOGICAL :: isExhausted

      !stat 1 means pass 0 means fail
      isExhausted = .FALSE.
      DO WHILE ( .NOT. isExhausted .AND. stat /= 1) 
        CALL critpoint_find(bdr,chg,opts,ions,stat)
        !CALL critpoint_find(bdr,chg,opts,ions)
        !PRINT *, "Stat came back as ", stat
        isExhausted = CheckExhaustion(opts)
        !PRINT *, "isExhausted came back as ", isExhausted
        IF (stat /= 1 .AND. .NOT. isExhausted) THEN
          CALL AdjustParameters(opts)
        END IF
      END DO
    END SUBROUTINE AutoCP

    FUNCTION CheckExhaustion(opts)
      ! This function checks if the parameter options are exhausted
      TYPE (options_obj) :: opts
      LOGICAL :: CheckExhaustion
      CheckExhaustion = .FALSE.
      IF (opts%autocp_flag) THEN
        IF (opts%cp_search_radius<=2 .AND. opts%par_newtonr <= 0.000000001 .AND. &
            opts%par_gradfloor <= 0.000000001) THEN
          IF (opts%enableDensityDescend  .AND. &
              opts%gradMode ) THEN
            CheckExhaustion = .TRUE.
            PRINT *, "The strictest setting has been used."
          ELSE IF (.NOT. opts%enableDensityDescend ) THEN
            opts%enableDensityDescend = .TRUE.
            ! In case parameters has been unnecessarily strict, reset everything
            opts%par_newtonr = 0.0001
            opts%par_gradfloor = 0.0001
            opts%cp_search_radius = 3
          ELSE IF (.NOT. opts%gradMode ) THEN
            opts%gradMode = .TRUE.
            opts%par_newtonr = 0.0001
            opts%par_gradfloor = 0.0001
            opts%cp_search_radius = 3
          ! Smoothening seems to be not worth the effort. Someone is doing it
          ! better than us.
          !ELSE IF (.NOT. opts%enableCHGCARSmoothening ) THEN
          !  PRINT *, "Enabling CHGCAR Smoothening"
          !  !PRINT *, "WARNING :: This function is only designed to work with &
          !  !  AFLOW"
          !  PRINT *, "Setting other parameters back to default values"
          !  opts%enableCHGCARSmoothening = .TRUE.
          !  opts%par_newtonr = 0.0001
          !  opts%par_gradfloor = 0.0001
          !  opts%cp_search_radius = 3
          END IF

        END IF
      ELSE 
        PRINT *, "Heuristic features not turned on. Exiting"
        CheckExhaustion = .TRUE.
      END IF
      RETURN 
    END FUNCTION checkExhaustion

    SUBROUTINE AdjustParameters(opts)
      TYPE (options_obj) :: opts
      IF (opts%par_newtonr > 0.000000001) THEN
        opts%par_newtonr = opts%par_newtonr * 0.01
      END IF
      IF ( opts%par_gradfloor > 0.000000001) THEN
        opts%par_gradfloor = opts%par_gradfloor * 0.01
      END IF
      IF ( opts%par_gradfloor <= 0.000000001 .AND. &
           opts%par_newtonr <= 0.000000001 .AND. opts%cp_search_radius > 2) THEN
        opts%cp_search_radius = opts%cp_search_radius - 1
      END IF
      PRINT *, "Search failed to satisfy constrains! Adjusting parameters:"
      PRINT *, "parnewtonr is now ", opts%par_newtonr
      PRINT *, "pargradfloor is now ", opts%par_gradfloor
      PRINT *, "CP search radius is now ", opts%cp_search_radius
    END SUBROUTINE AdjustParameters

    ! Output the last set of parameters used. If a calculationo is restarted, it
    ! doesn't go through the heuristic searching process again.
    SUBROUTINE OutputParameters()
     
    END SUBROUTINE OutPutParameters

  END MODULE
