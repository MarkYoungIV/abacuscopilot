!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for specifying input options
!-----------------------------------------------------------------------------------!
  MODULE options_mod
    USE kind_mod , ONLY : q2
    IMPLICIT NONE

    TYPE :: options_obj
      CHARACTER(LEN=128) :: chargefile, refchgfile
      REAL(q2) :: badertol, stepsize, vacval 
      INTEGER :: out_opt, out_auto = 0, out_cube = 1, out_chgcar4 = 2, out_chgcar5 = 3
      INTEGER :: in_opt, in_auto = 0, in_cube = 1, in_chgcar = 2, in_chgcar4 = 3,in_chgcar5 = 4
      INTEGER :: bader_opt, bader_offgrid = 0, bader_ongrid = 1, bader_neargrid = 2, bader_weight = 3
      INTEGER :: quit_opt, quit_max = 0, quit_known = 1
      INTEGER :: refine_edge_itrs
! refine_edge_itrs=-1 check points around the reassigned points during refinement
! refine_edge_itrs=-2 check every edge point during refinement
      INTEGER :: selanum, selbnum, sumanum, sumbnum
      INTEGER,ALLOCATABLE,DIMENSION(:) :: selavol, selbvol,sumavol,sumbvol
      LOGICAL :: vac_flag
      LOGICAL :: bader_flag, voronoi_flag, dipole_flag, ldos_flag
      LOGICAL :: print_all_bader, print_all_atom
      LOGICAL :: print_sel_bader, print_sel_atom
      LOGICAL :: print_sum_bader, print_sum_atom
      LOGICAL :: print_bader_index, print_atom_index
      LOGICAL :: verbose_flag, ref_flag, find_critpoints_flag
      LOGICAL :: print_surfaces_atoms
      REAL(q2) :: par_tem, par_gradfloor
      REAL(q2) :: cp_min_dinstance
      REAL(q2) :: par_newtonr, cp_search_radius
      LOGICAL :: ismolecule
      LOGICAL :: iscrystal
      LOGICAL :: debugMode
      LOGICAL :: dohes
      LOGICAL :: enableDensityDescend,enableCHGCARSmoothening
      LOGICAL :: gradMode
      !gradMode enables the gradient descend subroutine to replace the Newton's Method subroutine
      LOGICAL :: GD_magMode 
      LOGICAL :: autocp_flag ! Turns on heuristic features
      LOGICAL :: static_search ! Only initiate searches from coordinates read
      LOGICAL :: static_check ! reads in a list of CPs and runs through checks.
      LOGICAL :: aflow_sym ! Reads aflow.fgroup.relax.out to get symmetry info
      LOGICAL :: ignore_cp_conflict ! This forces ReduceCP to always finish.
   END TYPE options_obj

    PRIVATE
    PUBLIC :: get_options, options_obj

    CONTAINS

!-----------------------------------------------------------------------------------!
! get_options: Read any input flags and the charge density file name
!-----------------------------------------------------------------------------------!

    SUBROUTINE get_options(opts)

      TYPE(options_obj) :: opts
      LOGICAL :: existflag
      LOGICAL :: readchgflag
      INTEGER :: n, i, ip, m, it
      INTEGER :: sel, itmp, istart, iend
      CHARACTER(LEN=128) :: p
      CHARACTER*128 :: inc
      INTEGER :: COMMAND_ARGUMENT_COUNT

! Default values
      opts%par_tem = 0
      opts%cp_search_radius = 3
      opts%cp_min_dinstance = 1.5
      ! cp_min_dinstance is the criteria for determine if two points are identical
      ! due to spacial proximity. This is in lattice units.
      opts%par_newtonr = 0.1
      opts%par_gradfloor = 0.1
     ! opts%par_GDgradfloor = 0.2
      opts%ismolecule = .FALSE.
      opts%iscrystal = .FALSE.
      opts%dohes = .FALSE.
      opts%out_opt = opts%out_chgcar4
      opts%in_opt = opts%in_auto
      opts%debugMode = .FALSE.
    !   print options
      opts%vac_flag = .FALSE.
!      opts%vac_flag = .TRUE.
      opts%vacval = 1E-3
      opts%print_all_atom = .FALSE.
      opts%print_all_bader = .FALSE.
      opts%print_sel_atom = .FALSE.
      opts%print_sel_bader = .FALSE.
      opts%print_sum_atom = .FALSE.
      opts%print_sum_bader = .FALSE.
      opts%print_bader_index = .FALSE.
      opts%print_atom_index = .FALSE.
      ! end of print options
      opts%bader_opt = opts%bader_neargrid
      opts%quit_opt = opts%quit_known
      opts%refine_edge_itrs = -1
      opts%bader_flag = .TRUE.
      opts%voronoi_flag = .FALSE.
      opts%dipole_flag = .FALSE.
      opts%ldos_flag = .FALSE.
      opts%verbose_flag = .FALSE.
      opts%badertol = 1E-3
      opts%stepsize = 0.0_q2
      opts%ref_flag = .FALSE.
      opts%find_critpoints_flag = .FALSE.
      opts%autocp_flag = .FALSE.
      opts%print_surfaces_atoms = .FALSE.
      opts%enableDensityDescend = .FALSE.
      opts%enableCHGCARSmoothening = .FALSE.
      opts%gradMode = .FALSE. 
      opts%GD_magMode = .TRUE.
      opts%static_search = .FALSE.
      opts%aflow_sym = .FALSE.
      opts%ignore_cp_conflict = .FALSE.
!      n=IARGC()
      n=COMMAND_ARGUMENT_COUNT()
      IF (n == 0) THEN
        call write_options()
        STOP
      END IF

      ! Loop over all arguments
      m=0
      readchgflag = .FALSE.
      readopts: DO WHILE(m<n)
210        m=m+1
!        CALL GETARG(m,p)
        CALL GET_COMMAND_ARGUMENT(m,p)
        p=ADJUSTL(p)
        ip=LEN_TRIM(p)
        i=INDEX(p,'-')

        IF (i /= 1) THEN

          ! Not a flag, so read the charge density file name
          IF (readchgflag) THEN
            WRITE(*,'(A,A,A)') ' Option "',p(1:ip),'" is not valid'
            CALL write_options()
            STOP
          END IF
          opts%chargefile=p
          INQUIRE(FILE=opts%chargefile,EXIST=existflag)
          IF (.NOT. existflag) THEN
            WRITE(*,'(2X,A,A)') opts%chargefile(1:ip),' does not exist'
            STOP
          END IF
          readchgflag = .TRUE.

        ! Help
        ELSEIF (p(1:ip) == '-h') THEN
          CALL write_help()
          STOP
        ! Verbose
        ELSEIF (p(1:ip) == '-v') THEN
          opts%verbose_flag = .TRUE.
        ! critpoint options
        ELSEIF (p(1:ip) == '-cp') THEN
          opts%find_critpoints_flag = .TRUE.
          opts%vac_flag = .TRUE. ! set this to be default
        ELSEIF (p(1:ip) == '-autocp') THEN
          opts%autocp_flag = .TRUE.
        ELSEIF (p(1:ip) == '-static_search') THEN
          opts%static_search = .TRUE.
        ELSEIF (p(1:ip) == '-static_check') THEN
          opts%static_check = .TRUE.
        ELSEIF (p(1:ip) == '-molecule') THEN
          PRINT *, 'This system is assigned as a molecule'
          opts%ismolecule = .TRUE.
        ELSEIF (p(1:ip) == '-crystal') THEN
          PRINT *, 'This system is assigned as a crystal'
          opts%iscrystal = .TRUE.  
        ELSEIF (p(1:ip) == '-ignore_cp_conflict') THEN
          PRINT *, 'Reduce CP will run its full course even for problematic cases. '
          opts%ignore_cp_conflict = .TRUE.  
        ELSEIF (p(1:ip) == '-dohes') THEN
          ! Use diagonal terms in the hessian
          opts%dohes = .TRUE.  
        ELSEIF (p(1:ip) == '-debug') THEN
          PRINT *, 'Entering debug mode.'
          PRINT *, 'Reading debugConfig for instructions'
          opts%debugMode = .TRUE.  
        ELSEIF (p(1:ip) == '-gradMode') THEN
          PRINT *, 'Prioritizing gradient descent as default heuristic.'
          opts%gradMode = .TRUE.
        ELSEIF (p(1:ip) == '-smoothenCHGCAR') THEN
          opts%enableCHGCARSmoothening = .TRUE.  
        ELSEIF (p(1:ip) == '-aflow_sym') THEN
          opts%aflow_sym = .TRUE.  
        ELSEIF (p(1:ip) == '-DensityDescend') THEN
          opts%enableDensityDescend = .TRUE.  
        ELSEIF (p(1:ip) == '-partem') THEN
          m = m + 1
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc = ADJUSTL(inc)
          it=LEN_TRIM(inc)
          PRINT *, 'instructed to read partem'
          READ(inc,*) opts%par_tem
        ELSEIF (p(1:ip) == '-cp_search_radius') THEN
          m = m + 1
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc = ADJUSTL(inc)
          it=LEN_TRIM(inc)
          PRINT *, 'instructed to read partem'
          READ(inc,*) opts%cp_search_radius
        ELSEIF (p(1:ip) == '-cp_min_distance') THEN
          m = m + 1
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc = ADJUSTL(inc)
          it=LEN_TRIM(inc)
          READ(inc,*) opts%cp_min_dinstance
        ELSEIF (p(1:ip) == '-parnewtonr') THEN
          m = m + 1
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc = ADJUSTL(inc)
          it=LEN_TRIM(inc)
          READ(inc,*) opts%par_newtonr
        ELSEIF (p(1:ip) == '-pargradfloor') THEN
          m = m + 1
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc = ADJUSTL(inc)
          it=LEN_TRIM(inc)
          READ(inc,*) opts%par_gradfloor
        ! Vacuum options
        ELSEIF (p(1:ip) == '-vac') THEN
          m=m+1
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'AUTO' .OR. inc(1:it) == 'auto') THEN
            opts%vac_flag = .TRUE.
          ELSEIF (inc(1:it) == 'OFF' .OR. inc(1:it) == 'off') THEN
            opts%vac_flag = .FALSE.
          ELSE
             READ(inc,*) opts%vacval
             opts%vac_flag = .TRUE.
          END IF
        ! Bader options
        ELSEIF (p(1:ip) == '-b') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'OFFGRID' .OR. inc(1:it) == 'offgrid') THEN
            opts%bader_opt = opts%bader_offgrid
          ELSEIF (inc(1:it) == 'ONGRID' .OR. inc(1:it) == 'ongrid') THEN
            opts%bader_opt = opts%bader_ongrid
          ELSEIF (inc(1:it) == 'NEARGRID' .OR. inc(1:it) == 'neargrid') THEN
            opts%bader_opt = opts%bader_neargrid
          ELSEIF (inc(1:it) == 'WEIGHT' .OR. inc(1:it) == 'weight') THEN
            opts%bader_opt = opts%bader_weight
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF

        ! Quit options
        ELSEIF (p(1:ip) == '-m') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'MAX' .OR. inc(1:it) == 'max') THEN
            opts%quit_opt = opts%quit_max
          ELSEIF (inc(1:it) == 'KNOWN' .OR. inc(1:it) == 'known') THEN
            opts%quit_opt = opts%quit_known
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF

        ! Print options
        ELSEIF (p(1:ip) == '-p') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'BADER_INDEX' .OR. inc(1:it) == 'bader_index') THEN
            opts%print_bader_index = .TRUE.
          ELSEIF (inc(1:it) == 'ATOM_INDEX' .OR. inc(1:it) == 'atom_index') THEN
            opts%print_atom_index = .TRUE.
          ELSEIF (inc(1:it) == 'ALL_BADER' .OR. inc(1:it) == 'all_bader') THEN
            opts%print_all_bader = .TRUE.
          ELSEIF (inc(1:it) == 'ALL_ATOM' .OR. inc(1:it) == 'all_atom') THEN
            opts%print_all_atom = .TRUE.
          ELSEIF (inc(1:it) == 'SURFACES_ATOMS' .OR. inc(1:it) == 'surfaces_atoms') THEN
            opts%print_surfaces_atoms = .TRUE.
          ELSEIF (inc(1:it) == 'SEL_BADER' .OR. inc(1:it) == 'sel_bader') THEN
            opts%print_sel_bader = .TRUE.
            opts%selbnum=0
            m=m+1
            CALL GET_COMMAND_ARGUMENT(m,inc)
            inc=ADJUSTL(inc)
            it=LEN_TRIM(inc)
            itmp = INDEX(inc,"-")
            IF (itmp .GT. 0) THEN
              READ (inc(1:itmp-1),'(I10)',ERR=110) istart
              READ (inc(itmp+1:it),'(I10)',ERR=110) iend
              ALLOCATE(opts%selbvol(iend-istart+1))
              DO sel = istart, iend
                opts%selbnum=opts%selbnum+1
                opts%selbvol(opts%selbnum)=sel
              END DO
              IF(m==n) EXIT readopts
              GO TO 210
   110        m=m-1
              EXIT
            ELSE
              m=m-1
              ALLOCATE(opts%selbvol(n))
              DO
                m=m+1
                CALL GET_COMMAND_ARGUMENT(m,inc)
                inc=ADJUSTL(inc)
                it=LEN_TRIM(inc)
                READ (inc(1:it),'(I10)',ERR=120) sel
                opts%selbnum=opts%selbnum+1
                opts%selbvol(opts%selbnum)=sel
                IF(m==n) EXIT readopts
                CYCLE
     120        m=m-1
                EXIT
              END DO
            END IF
          ELSEIF (inc(1:it) == 'SEL_ATOM' .OR. inc(1:it) == 'sel_atom') THEN
            opts%print_sel_atom = .TRUE.
            opts%selanum=0
            m=m+1
            CALL GET_COMMAND_ARGUMENT(m,inc)
            inc=ADJUSTL(inc)
            it=LEN_TRIM(inc)
            itmp = INDEX(inc,"-")
            IF (itmp .GT. 0) THEN
              READ (inc(1:itmp-1),'(I10)',ERR=130) istart
              READ (inc(itmp+1:it),'(I10)',ERR=130) iend
              ALLOCATE(opts%selavol(iend-istart+1))
              DO sel = istart, iend
                opts%selanum=opts%selanum+1
                opts%selavol(opts%selanum)=sel
              END DO
              IF(m==n) EXIT readopts
              GO TO 210
   130        m=m-1
              EXIT
            ELSE
              m=m-1
              ALLOCATE(opts%selavol(n))
              DO
                m=m+1
                CALL GET_COMMAND_ARGUMENT(m,inc)
                inc=ADJUSTL(inc)
                it=LEN_TRIM(inc)
                READ (inc(1:it),'(I10)',ERR=140) sel
                opts%selanum=opts%selanum+1
                opts%selavol(opts%selanum)=sel
                IF(m==n) EXIT readopts
                CYCLE
     140        m=m-1
                EXIT
              END DO
            END IF
          ELSEIF (inc(1:it) == 'SUM_ATOM' .OR. inc(1:it) == 'sum_atom') THEN
            opts%print_sum_atom = .TRUE.
            opts%sumanum=0
            m=m+1
            CALL GET_COMMAND_ARGUMENT(m,inc)
            inc=ADJUSTL(inc)
            it=LEN_TRIM(inc)
            itmp = INDEX(inc,"-")
            IF (itmp .GT. 0) THEN
              READ (inc(1:itmp-1),'(I10)',ERR=150) istart
              READ (inc(itmp+1:it),'(I10)',ERR=150) iend
              ALLOCATE(opts%sumavol(iend-istart+1))
              DO sel = istart, iend
                opts%sumanum=opts%sumanum+1
                opts%sumavol(opts%sumanum)=sel
              END DO
              IF(m==n) EXIT readopts
              GO TO 210
   150        m=m-1
              EXIT
            ELSE
              m=m-1
              ALLOCATE(opts%sumavol(n))
              DO
                m=m+1
                CALL GET_COMMAND_ARGUMENT(m,inc)
                inc=ADJUSTL(inc)
                it=LEN_TRIM(inc)
                READ (inc(1:it),'(I10)',ERR=160) sel
                opts%sumanum=opts%sumanum+1
                opts%sumavol(opts%sumanum)=sel
                IF(m==n) EXIT readopts
                CYCLE
     160        m=m-1
                EXIT
              END DO
            END IF
          ELSEIF (inc(1:it) == 'SUM_BADER' .OR. inc(1:it) == 'sum_bader') THEN
            opts%print_sum_bader = .TRUE.
            opts%sumbnum=0
            m=m+1
            CALL GET_COMMAND_ARGUMENT(m,inc)
            inc=ADJUSTL(inc)
            it=LEN_TRIM(inc)
            itmp = INDEX(inc,"-")
            IF (itmp .GT. 0) THEN
              READ (inc(1:itmp-1),'(I10)',ERR=170) istart
              READ (inc(itmp+1:it),'(I10)',ERR=170) iend
              ALLOCATE(opts%sumbvol(iend-istart+1))
              DO sel = istart, iend
                opts%sumbnum=opts%sumbnum+1
                opts%sumbvol(opts%sumbnum)=sel
              END DO
              IF(m==n) EXIT readopts
              GO TO 210
   170        m=m-1
              EXIT
            ELSE
              m=m-1
              ALLOCATE(opts%sumbvol(n))
              DO
                m=m+1
                CALL GET_COMMAND_ARGUMENT(m,inc)
                inc=ADJUSTL(inc)
                it=LEN_TRIM(inc)
                READ (inc(1:it),'(I10)',ERR=180) sel
                opts%sumbnum=opts%sumbnum+1
                opts%sumbvol(opts%sumbnum)=sel
                IF(m==n) EXIT readopts
                CYCLE
     180        m=m-1
                EXIT
              END DO
            END IF
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF
        ! Output file type
!        ELSEIF (p(1:ip) == '-o') THEN
!          m=m+1
!          CALL GETARG(m,inc)
!          CALL GET_COMMAND_ARGUMENT(m,inc)
!          inc=ADJUSTL(inc)
!          it=LEN_TRIM(inc)
!          IF (inc(1:it) == 'CUBE' .OR. inc(1:it) == 'cube') THEN
!            opts%out_opt = opts%out_cube
!          ELSEIF (inc(1:it) == 'CHGCAR' .OR. inc(1:it) == 'chgcar') THEN
!            opts%out_opt = opts%out_chgcar
!          ELSE
!            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
!            STOP
!          END IF  

        ! Calculation options
        ELSEIF (p(1:ip) == '-c') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'ALL' .OR. inc(1:it) == 'all') THEN
            opts%bader_flag = .TRUE.
            opts%voronoi_flag = .TRUE.
            opts%dipole_flag = .TRUE.
            opts%ldos_flag = .TRUE.
          ELSEIF (inc(1:it) == 'BADER' .OR. inc(1:it) == 'bader') THEN
            opts%bader_flag = .TRUE.
          ELSEIF (inc(1:it) == 'VORONOI' .OR. inc(1:it) == 'voronoi') THEN
            opts%voronoi_flag = .TRUE.
          ELSEIF (inc(1:it) == 'DIPOLE' .OR. inc(1:it) == 'dipole') THEN
            opts%dipole_flag = .TRUE.
          ELSEIF (inc(1:it) == 'LDOS' .OR. inc(1:it) == 'ldos') THEN
            opts%ldos_flag = .TRUE.
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF
        ELSEIF (p(1:ip) == '-n') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'ALL' .OR. inc(1:it) == 'all') THEN
            opts%bader_flag = .FALSE.
            opts%voronoi_flag = .FALSE.
            opts%dipole_flag = .FALSE.
            opts%ldos_flag = .FALSE.
          ELSEIF (inc(1:it) == 'BADER' .OR. inc(1:it) == 'bader') THEN
            opts%bader_flag = .FALSE.
          ELSEIF (inc(1:it) == 'VORONOI' .OR. inc(1:it) == 'voronoi') THEN
            opts%voronoi_flag = .FALSE.
          ELSEIF (inc(1:it) == 'DIPOLE' .OR. inc(1:it) == 'dipole') THEN
            opts%dipole_flag = .FALSE.
          ELSEIF (inc(1:it) == 'LDOS' .OR. inc(1:it) == 'ldos') THEN
            opts%ldos_flag = .FALSE.
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          ENDIF

        ! Input file type
        ELSEIF (p(1:ip) == '-i') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'CUBE' .OR. inc(1:it) == 'cube') THEN
            opts%in_opt=opts%in_cube
          ELSEIF (inc(1:it) == 'CHGCAR' .OR. inc(1:it) == 'chgcar') THEN
            opts%in_opt=opts%in_chgcar
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF

        ! Bader tolerance
        ELSEIF (p(1:ip) == '-t') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          READ(inc,*) opts%badertol

        ! Refine edge iterations  -- change this to a flag once working
        ELSEIF (p(1:ip) == '-r') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc) 
          IF (inc(1:it) == 'AUTO' .OR. inc(1:it) == 'auto') THEN
            opts%refine_edge_itrs=-1
          ELSE
            READ(inc,*) opts%refine_edge_itrs
          END IF

        ! Step size
        ELSEIF (p(1:ip) == '-s') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          READ(inc,*) opts%stepsize

        ! Do analysis with a reference charge
        ELSEIF (p(1:ip) == '-ref') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'NONE' .OR. inc(1:it) == 'none') THEN
            opts%ref_flag = .FALSE.
          ELSE
            opts%ref_flag = .TRUE.
            opts%refchgfile = inc(1:it)
          END IF

        ! Unknown flag
        ELSE
          WRITE(*,'(A,A,A)') ' Unknown option flag "',p(1:ip),'"'
          STOP
        END IF

      END DO readopts

    ! If no file name, we die
    IF (.NOT. readchgflag) THEN
      WRITE(*,*) ' ERROR: Did not read a charge file name in the arguments'
      CALL write_options()
      STOP
    ENDIF

    ! Default to no edge refinement for the ongrid algorithm
    IF (opts%bader_opt==opts%bader_ongrid) THEN
      opts%refine_edge_itrs=0
    END IF
   
    IF (opts%print_sel_atom .AND. opts%selanum==0) THEN
      WRITE(*,'(/,A)') 'NO ATOMIC VOLUMES SELECTED'
      STOP
    END IF

    IF (opts%print_sel_bader .AND. opts%selbnum==0) THEN
      WRITE(*,'(/,A)') 'NO BADER VOLUMES SELECTED'
      STOP
    END IF

    IF (opts%print_sum_bader .AND. opts%sumbnum==0) THEN
      WRITE(*,'(/,A)') 'NO BADER VOLUMES SELECTED'
      STOP
    END IF

    IF (opts%print_sum_atom .AND. opts%sumanum==0) THEN
      WRITE(*,'(/,A)') 'NO BADER VOLUMES SELECTED'
      STOP
    END IF
    

    RETURN
    END SUBROUTINE get_options

!-----------------------------------------------------------------------------------!
! write_opts: write flag options
!-----------------------------------------------------------------------------------!

    SUBROUTINE write_options()

      WRITE(*,*) ''
      WRITE(*,*) 'Usage:'
      WRITE(*,*) '   bader [ -c bader | voronoi ]'
      WRITE(*,*) '         [ -n bader | voronoi ]'
      WRITE(*,*) '         [ -b neargrid | ongrid | weight]'
      WRITE(*,*) '         [ -r refine_edge_method ]'
      WRITE(*,*) '         [ -ref reference_charge ]'
      WRITE(*,*) '         [ -vac off | auto | vacuum_density ]'
      WRITE(*,*) '         [ -m known | max ]'
      WRITE(*,*) '         [ -p all_atom | all_bader ]'
      WRITE(*,*) '         [ -p sel_atom | sel_bader ] [ volume list or range ]'
      WRITE(*,*) '         [ -p sum_atom | sum_bader ] [ volume list or range ]'
      WRITE(*,*) '         [ -p atom_index | bader_index ]'
      WRITE(*,*) '         [ -p surfaces_atoms ]'
      WRITE(*,*) '         [ -i cube | chgcar ]'
      WRITE(*,*) '         [ -h ] [ -v ]'
      WRITE(*,*) '         chargefile'
      WRITE(*,*) ''

    END SUBROUTINE write_options

!-----------------------------------------------------------------------------------!
! write_help: write help
!-----------------------------------------------------------------------------------!

    SUBROUTINE write_help()

      WRITE(*,*) ''
      WRITE(*,*) 'Description of flags'
      WRITE(*,*) ''
!      WRITE(*,*) '   -c | -n  < bader | voronoi | dipole | ldos >'
      WRITE(*,*) '   -c | -n  < bader | voronoi >'
      WRITE(*,*) '        Turn on [-c] or off [-n] the following calculations'
      WRITE(*,*) '           bader: Bader atoms in molecules (default)'
      WRITE(*,*) '           voronoi: population analysis based on distance'
!      WRITE(*,*) '           dipole: multiple moments in Bader volumes'
!      WRITE(*,*) '           ldos: local density of states in Bader volumes'
      WRITE(*,*) ''
      WRITE(*,*) '   -b < neargrid | ongrid | weight >'
      WRITE(*,*) '        Use the default near-grid bader partitioning, the'
      WRITE(*,*) '        original on-grid based algorithm, or the weight method'
      WRITE(*,*) '        of Yu and Trinkle'
      WRITE(*,*) ''
!      WRITE(*,*) '   -s < stepsiz >'
!      WRITE(*,*) '        Steepest ascent trajectory step size.  This parameter is'
!      WRITE(*,*) '        (only) used for the default offgrid Bader analysis.  If'
!      WRITE(*,*) '        not specified, the stepsize is set to the minimum distance'
!      WRITE(*,*) '        between charge density grid points.'
!      WRITE(*,*) ''
      WRITE(*,*) '   -r < refine_edge_method >'
      WRITE(*,*) '        By default (-r -1) , only the points around reassigned'
      WRITE(*,*) '        points are checked during refinements. The old method, '
      WRITE(*,*) '        which checks every edge point during each refinement, can'
      WRITE(*,*) '        be enabled using the -r -2 switch:'
      WRITE(*,*) '           bader -r -2 CHGCAR'
      WRITE(*,*) ''
      WRITE(*,*) '   -ref < reference_charge >'
      WRITE(*,*) '        Use a reference charge file to do the Bader partitioning.'
      WRITE(*,*) '        This is the recommended way to analyze vasp output files:'
      WRITE(*,*) '           bader CHGCAR -ref CHGCAR_total'
      WRITE(*,*) '        where CHGCAR_total is the sum of AECCAR0 and AECCAR2.'
      WRITE(*,*) ''
      WRITE(*,*) '   -vac < off | auto | vacuum_density >'
      WRITE(*,*) '        Assign low density points to vacuum.'
      WRITE(*,*) '          auto: vacuum density cutoff is 1E-3 e/Ang^3 by default'
      WRITE(*,*) '          off: do not assign low density points to a vacuum volume'
      WRITE(*,*) '          vacuum_density: maximum density assigned to a vacuum volume'
      WRITE(*,*) ''
      WRITE(*,*) '   -m < known | max >'
      WRITE(*,*) '        Determines how trajectories terminate'
      WRITE(*,*) '           known: stop when a point is surrounded by known points' 
      WRITE(*,*) '           max: stop only when a charge density maximum is reached '
      WRITE(*,*) ''
      WRITE(*,*) '   -p < all_atom | all_bader >'
      WRITE(*,*) '        Print calculated Bader volumes'
      WRITE(*,*) '           all_atom: all atomic volumes'
      WRITE(*,*) '           all_bader: all Bader volumes'
      WRITE(*,*) ''
      WRITE(*,*) '   -p < sel_atom | sel_bader > <volume list or range>'
      WRITE(*,*) '        Print calculated Bader volumes'
      WRITE(*,*) '           sel_atom: atomic volume(s) around the selected atom(s)'
      WRITE(*,*) '           sel_bader: selected Bader volumes'
      WRITE(*,*) ''
      WRITE(*,*) '   -p < sum_atom | sum_bader > <volume list or range>'
      WRITE(*,*) '        Print calculated Bader volumes'
      WRITE(*,*) '           sum_atom: sum of atomic volume(s) around the selected atom(s)'
      WRITE(*,*) '           sum_bader: sum of selected Bader volumes'
      WRITE(*,*) ''
      WRITE(*,*) '   -p < atom_index | bader_index >'
      WRITE(*,*) '        Print index of atomic or Bader volumes'
      WRITE(*,*) '           atom_index: print atomic volume indicies'
      WRITE(*,*) '           bader_index: print Bader volume indicies'
      WRITE(*,*) ''
      WRITE(*,*) '   -p < surfaces_atoms >'
      WRITE(*,*) '        Print a file with a list of surfaces coordinates and atoms'
      WRITE(*,*) '        forming this surface'
      WRITE(*,*) ''
      WRITE(*,*) '   -i < cube | chgcar >'
      WRITE(*,*) '        Input charge density file type.  If not specified, the'
      WRITE(*,*) '        program will try to determine the charge density file type'
      WRITE(*,*) '        automatically.'
      WRITE(*,*) ''
      WRITE(*,*) '   -h'
      WRITE(*,*) '        Help.'
      WRITE(*,*) ''
      WRITE(*,*) '   -v'
      WRITE(*,*) '        Verbose output.'
      WRITE(*,*) ''
      WRITE(*,*) '   -cp'
      WRITE(*,*) '        Find critical points of the charge density'
      WRITE(*,*) '        Calculate eigenvalues and eigenvectors at those points'
      WRITE(*,*) '        Store results in file named CPF##.dat'
      WRITE(*,*) ''    
      !WRITE(*,*) '   -parnewtonr'
      !WRITE(*,*) '        Criteria for stopping Newton Raphson method trajectory for &
      !                    small stepsize. Default values should be OK.'
      !WRITE(*,*) '   -pargradfloor '
      !WRITE(*,*) '        Criterial for stopping Newtons method for small &
      !                    gradient. Default values should be OK.'
      !WRITE(*,*) '   -partem '
      !WRITE(*,*) '        Proximity criteria for initializing Newtons method &
      !                    trajectory. Default values should be OK.'
      !WRITE(*,*) '   -pardistance '
      !WRITE(*,*) '        Proximity criteria in Angstrom for declaring two CP &
      !                    duplicates of each other. Default values shoud &
      !                    be OK.'
      !WRITE(*,*) '   -parsr'
      !WRITE(*,*) '        Search radius. If a Newton Rhapson trajectory &
      !                    initiated at a point, no more trajectory will &
      !                    initiate within this radius. Default is 0. Higher &
      !                    values will improve performance. Value higher than &
      !                    3 will likely cause inaccurate results.'

    END SUBROUTINE write_help


!-----------------------------------------------------------------------------------!

  END module options_mod
