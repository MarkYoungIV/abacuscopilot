!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge with the Bader atom in molecules approach
!-----------------------------------------------------------------------------------!

MODULE bader_mod

  USE kind_mod
  USE matrix_mod
  USE options_mod
  USE ions_mod
  USE charge_mod
  USE io_mod
  USE chgcar_mod

  IMPLICIT NONE

! Variable descriptions:

! volnum: Bader volume number for each grid point
! volpos: position of maximum in each Bader volume
! colchg: integrated charge of each Bader volume
! ionchg: integrated charge over all Bader volumes associated with each ion
! iondist: distance from each Bader maximum to the nearest ion
! nnion: index of the nearst ion used to calculate iondist
! path: array of points along the current charge density maximization
! minsurfdist: minimum distance from the Bader surface to the included ion

  TYPE bader_obj
    REAL(q2), ALLOCATABLE, DIMENSION(:,:) :: volpos_lat, volpos_car, volpos_dir
    REAL(q2), ALLOCATABLE, DIMENSION(:) :: volchg, iondist, ionchg, minsurfdist, ionvol
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: volnum, known
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: path
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nnion
    REAL(q2) :: stepsize, tol
    REAL(q2) :: vacchg, vacvol
    INTEGER nvols, pnum, bnum, pdim, bdim, refine_edge_itrs
  END TYPE

  PRIVATE
  PUBLIC :: bader_obj
  PUBLIC :: bader_calc, bader_mindist, bader_output, write_all_atom, write_all_bader
  PUBLIC :: write_sel_atom, write_sel_bader, write_atom_index, write_bader_index
  PUBLIC :: write_sum_atom, write_sum_bader, assign_chg2atom, cal_atomic_vol
  PUBLIC :: reallocate_volpos, bader_check_partitioning
  PUBLIC :: write_surfaces_atoms, is_vol_edge

! scell_dir: flag if we have a doubling of the cell in direction x,y or z
  LOGICAL, DIMENSION(3) :: scell_dir = [.FALSE., .FALSE., .FALSE.]

  CONTAINS

!-----------------------------------------------------------------------------------!
! bader_calc: Calculate the Bader volumes and integrate to give the total
!   charge in each volume.
!-----------------------------------------------------------------------------------!

  SUBROUTINE bader_calc(bdr,ions,chgval,opts)

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chgval
    TYPE(options_obj) :: opts

    INTEGER, DIMENSION(3) :: p, ptemp
    INTEGER :: n1, n2, n3, i, path_volnum, tenths_done
    INTEGER :: cr, count_max, t1, t2
    INTEGER :: ref_itrs = 1

    REAL(q2), DIMENSION(3) :: voxlen
    REAL(q2) :: vol
    TYPE(charge_obj) :: chgtemp
    TYPE(ions_obj) :: ionstemp

    INTEGER, ALLOCATABLE :: sc_atom_map(:)
    TYPE(ions_obj) :: ionsscell
    TYPE(charge_obj) :: chgscell

    IF (opts%ref_flag) THEN
      CALL read_charge_ref(ionstemp,chgtemp,opts)
      ! If we are re-running the calcuation in a supercell, the reference must be expandet, too
      IF ( ANY(scell_dir(:)) ) THEN
        CALL build_scell( chgscell, ionsscell, chgtemp, ionstemp, scell_dir, sc_atom_map ) 
        chgtemp = chgscell
      END IF
      ! Assert that the charge and reference files are the same dimension
      IF ((chgval%npts(1) /= chgtemp%npts(1)) .OR. &
          (chgval%npts(2) /= chgtemp%npts(2)) .OR. &
          (chgval%npts(3) /= chgtemp%npts(3))) THEN
         WRITE(*,'(/,2x,A,/)') &
           'The dimensions of the primary and reference charge densities must be the same, stopping.'
         STOP
      END IF
    ELSE
      chgtemp = chgval
    END IF

    CALL SYSTEM_CLOCK(t1,cr,count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING BADER CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    ! copy bader variable from opts
    bdr%tol = opts%badertol
    bdr%stepsize = opts%stepsize
    IF (opts%stepsize == 0) THEN  ! check for transpose error
      DO i=1,3
        voxlen(i) = SQRT(SUM(chgval%lat2car(:,i)*chgval%lat2car(:,i)))
      END DO
      bdr%stepsize = MINVAL(voxlen(:))
    END IF

    bdr%bdim = 64
    bdr%pdim = 64
    ALLOCATE(bdr%volpos_lat(bdr%bdim,3)) ! will be expanded as needed
    ALLOCATE(bdr%path(bdr%pdim,3))
    ALLOCATE(bdr%volnum(chgval%npts(1), chgval%npts(2), chgval%npts(3)))
    ALLOCATE(bdr%known(chgval%npts(1), chgval%npts(2), chgval%npts(3)))
!    bdr%volchg=0._q2
    bdr%volnum = 0
    bdr%known = 0
    bdr%bnum = 0
    bdr%nvols = 0  ! true number of Bader volumes

    ! find vacuum points, get the charge and volume
    bdr%vacchg = 0.0_q2
    bdr%vacvol = 0.0_q2

    vol = matrix_volume(ions%lattice)
    IF (opts%vac_flag) THEN
      DO n1=1,chgval%npts(1)
        DO n2=1,chgval%npts(2)
          DO n3=1,chgval%npts(3)
            IF (ABS(rho_val(chgtemp,n1,n2,n3)/vol) <= opts%vacval) THEN
!            IF (ABS(rho_val(chgval,n1,n2,n3)/vol) <= opts%vacval) THEN
!            IF (rho_val(chgval,n1,n2,n3)/vol <= opts%vacval) THEN
               bdr%volnum(n1,n2,n3) = -1
               bdr%vacchg = bdr%vacchg + chgval%rho(n1,n2,n3)
               bdr%vacvol = bdr%vacvol + 1
            END IF
          END DO
        END DO
      END DO
    END IF
    bdr%vacchg = bdr%vacchg/REAL(chgval%nrho,q2)
    bdr%vacvol = bdr%vacvol*vol/chgval%nrho

    tenths_done = 0
    DO n1=1,chgval%npts(1)
      IF ((n1*10/chgval%npts(1)) > tenths_done) THEN
        tenths_done = (n1*10/chgval%npts(1))
        WRITE(*,'(A,$)') '**'
      END IF
      DO n2=1,chgval%npts(2)
        DO n3=1,chgval%npts(3)
          p = (/n1,n2,n3/)
          IF (bdr%volnum(p(1),p(2),p(3)) == 0) THEN
            IF (opts%bader_opt == opts%bader_offgrid) THEN
              CALL max_offgrid(bdr,chgtemp,p)
            ELSEIF (opts%bader_opt == opts%bader_ongrid) THEN
              CALL max_ongrid(bdr,chgtemp,p)
            ELSE 
              CALL max_neargrid(bdr,chgtemp,opts,p)
            END IF
            path_volnum = bdr%volnum(p(1),p(2),p(3))

            ! final point has not been assigned, assign new maximum
            IF (path_volnum == 0) THEN
!              print*,'new max:',p
              IF (bdr%bnum >= bdr%bdim) THEN
                CALL reallocate_volpos(bdr,bdr%bdim*2)
              END IF
              bdr%bnum = bdr%bnum+1
              path_volnum = bdr%bnum
              bdr%volpos_lat(bdr%bnum,:) = REAL(p,q2)
            END IF

            ! assign all points along the trajectory
            DO i = 1,bdr%pnum
              ptemp = (/bdr%path(i,1),bdr%path(i,2),bdr%path(i,3)/)
              IF(bdr%volnum(ptemp(1),ptemp(2),ptemp(3)) /= -1) THEN
                bdr%volnum(ptemp(1),ptemp(2),ptemp(3)) = path_volnum
              END IF
              IF (opts%bader_opt /= opts%bader_ongrid) THEN
                ! GH: how does this work?  Can it ever get here with bdr%known==2?
                ! If so, the path should have quit
                IF (bdr%known(ptemp(1),ptemp(2),ptemp(3)) /= 2) THEN
                  bdr%known(ptemp(1),ptemp(2),ptemp(3)) = 0
                END IF
              END IF
              IF (opts%quit_opt==opts%quit_known .AND. opts%bader_opt/=opts%bader_ongrid) THEN
                CALL assign_surrounding_pts(bdr,chgtemp,ptemp)
!                CALL assign_surrounding_pts2(bdr,chgtemp,ptemp)
              END IF
            END DO

          END IF
        END DO
      END DO
    END DO
    WRITE(*,*) ''

    IF (opts%vac_flag) THEN
      DO n1 = 1, chgval%npts(1)
        DO n2 = 1, chgval%npts(2)
          DO n3 = 1, chgval%npts(3)
            IF(bdr%volnum(n1,n2,n3) == -1) bdr%volnum(n1,n2,n3) = bdr%bnum + 1
          END DO
        END DO
      END DO
    END IF

    ! make a temporary variable which can be changed to indicate convergence
    opts%refine_edge_itrs = 0
    bdr%refine_edge_itrs = opts%refine_edge_itrs

!    IF(opts%refine_edge_itrs == -1 .OR. opts%refine_edge_itrs == -2) THEN
    IF(opts%refine_edge_itrs <= 0) THEN
      WRITE(*,'(/,2x,A)') 'REFINING AUTOMATICALLY'
      DO
        WRITE(*,'(2x,A,I2)') 'ITERATION:',ref_itrs
        CALL refine_edge(bdr,chgtemp,opts,ions,ref_itrs)
        IF (bdr%refine_edge_itrs == 0) EXIT
        ref_itrs = ref_itrs + 1
      END DO
    ELSEIF (opts%refine_edge_itrs > 0) THEN
      WRITE(*,'(/,2x,A)') 'REFINING EDGE'
      DO i=1,opts%refine_edge_itrs
        WRITE(*,'(2x,A,I2)') 'ITERATION:',i
        CALL refine_edge(bdr,chgtemp,opts,ions,i)
      END DO
    ENDIF

    ! total number of bader volumes is now known
    bdr%nvols = bdr%bnum
    CALL reallocate_volpos(bdr, bdr%nvols)
    ALLOCATE(bdr%volpos_dir(bdr%nvols,3))
    ALLOCATE(bdr%volpos_car(bdr%nvols,3))
    DO i=1,bdr%nvols
      bdr%volpos_dir(i,:) = lat2dir(chgtemp, bdr%volpos_lat(i,:))
      bdr%volpos_car(i,:) = lat2car(chgtemp, bdr%volpos_lat(i,:))
    END DO

    ! Sum up the charge included in each volume
    ALLOCATE(bdr%volchg(bdr%nvols))
    bdr%volchg = 0._q2
    DO n1 = 1, chgval%npts(1)
      DO n2 = 1, chgval%npts(2)
        DO n3 = 1, chgval%npts(3)
          IF (bdr%volnum(n1,n2,n3) == bdr%nvols+1) CYCLE
          bdr%volchg(bdr%volnum(n1,n2,n3)) = &
          &  bdr%volchg(bdr%volnum(n1,n2,n3)) + chgval%rho(n1,n2,n3)
        END DO
      END DO
    END DO
    bdr%volchg = bdr%volchg/REAL(chgval%nrho,q2)

    ALLOCATE(bdr%nnion(bdr%nvols), bdr%iondist(bdr%nvols), bdr%ionchg(ions%nions))
    CALL assign_chg2atom(bdr,ions,chgval)
    CALL cal_atomic_vol(bdr,ions,chgval)

    DEALLOCATE(bdr%path)

    CALL SYSTEM_CLOCK(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F10.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE bader_calc

!-----------------------------------------------------------------------------------!
! max_offgrid:  From the point p, do a maximization in the charge density
!-----------------------------------------------------------------------------------!

  SUBROUTINE max_offgrid(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p
    INTEGER,DIMENSION(3) :: pm
    REAL(q2) :: cp,cm 
    REAL(q2),DIMENSION(3) :: r

    bdr%pnum = 1
    bdr%path(bdr%pnum,:) = p
    IF (is_max(chg,p)) THEN
      write(*,*) '   max found (init)'
      RETURN
    END IF
    r = REAL(p,q2)

!    write(*,*) 'max: ',p
    DO
      cp = rho_val(chg,p(1),p(2),p(3))
      CALL step_offgrid(bdr,chg,r)
      pm = to_lat(chg,r)
      cm = rho_val(chg,pm(1),pm(2),pm(3))

      IF (cm < cp) EXIT
      p = pm

      ! if the point is new, add it to the path
      IF (.NOT.ALL(p(:) == bdr%path(bdr%pnum,:))) THEN
        IF (bdr%pnum >= bdr%pdim) THEN
          CALL reallocate_path(bdr,bdr%pdim*2)
        ENDIF
        bdr%pnum = bdr%pnum+1
        CALL pbc(p,chg%npts)
        bdr%path(bdr%pnum,:) = p
      END IF

      ! quit if this is a maximum or a known volume number
      IF (is_max(chg,p) .OR. known_volnum(bdr,chg,r) /= 0) EXIT

    END DO
    
  RETURN
  END SUBROUTINE max_offgrid

!-----------------------------------------------------------------------------------!
!  step_offgrid: step a distance of StepSize along rho_grad
!-----------------------------------------------------------------------------------!

  SUBROUTINE step_offgrid(bdr,chg,r)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    REAL(q2), DIMENSION(3), INTENT(INOUT) :: r

    REAL(q2), DIMENSION(3) :: grad,dr_car,dr_lat
    REAL(q2) :: rho

    grad = rho_grad(chg,r,rho)
    dr_car = grad*bdr%stepsize/SQRT(SUM(grad*grad))
!    CALL matrix_vector(chg%car2lat,dr_car,dr_lat)
    dr_lat = MATMUL(chg%car2lat,dr_car)
    r = r+dr_lat

  RETURN
  END SUBROUTINE step_offgrid

!-----------------------------------------------------------------------------------!
! max_ongrid:  From the point p do a maximization on the charge density grid and
!   assign the maximum found to the volnum array.
!-----------------------------------------------------------------------------------!

  SUBROUTINE max_ongrid(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER, DIMENSION(3), INTENT(INOUT) :: p

    bdr%pnum = 1
    bdr%path(bdr%pnum,:) = p
    DO
      CALL step_ongrid(chg,p)
      ! if we didn't move, we're at a maximum
      IF (ALL(p == bdr%path(bdr%pnum,:))) EXIT
      ! otherwise, add point to path
      IF (bdr%pnum >= bdr%pdim) THEN
        CALL reallocate_path(bdr,bdr%pdim*2)
      ENDIF
      bdr%pnum = bdr%pnum+1
      CALL pbc(p,chg%npts)
      bdr%path(bdr%pnum,:) = p
!      IF(bdr%volnum(p(1),p(2),p(3)) /= 0) EXIT
!GH: modified to deal with refine_edge
      IF (bdr%volnum(p(1),p(2),p(3)) > 0) EXIT
    END DO

  RETURN
  END SUBROUTINE max_ongrid

!-----------------------------------------------------------------------------------!
!  step_ongrid:  Do a single iteration of a maximization on the charge density 
!    grid from the point (px,py,pz).  Return a logical indicating if the current
!    point is a charge density maximum.
!-----------------------------------------------------------------------------------!

  SUBROUTINE step_ongrid(chg,p)

    TYPE(charge_obj) :: chg
    INTEGER, DIMENSION(3), INTENT(INOUT) :: p

    REAL(q2) :: rho_max,rho_tmp,rho_ctr
    INTEGER, DIMENSION(3) :: pt,pm
    INTEGER :: d1,d2,d3

    pm = p
    rho_ctr = rho_val(chg,p(1),p(2),p(3))
    rho_max = rho_ctr
    DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p+(/d1,d2,d3/)
          rho_tmp = rho_val(chg,pt(1),pt(2),pt(3))
          rho_tmp = rho_ctr+(rho_tmp-rho_ctr)*chg%lat_i_dist(d1,d2,d3)
          IF (rho_tmp > rho_max) THEN
            rho_max = rho_tmp
            pm = pt
          END IF
        END DO
      END DO
    END DO
    CALL pbc(pm,chg%npts)
    p = pm

  RETURN
  END SUBROUTINE step_ongrid

!-----------------------------------------------------------------------------------!
! max_neargrid:  From the point p do a maximization on the charge density grid and
!   assign the maximum found to the volnum array.
!-----------------------------------------------------------------------------------!

  SUBROUTINE max_neargrid(bdr,chg,opts,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts

    INTEGER, DIMENSION(3), INTENT(INOUT) :: p

    bdr%pnum=1
    bdr%path(bdr%pnum,:)=p
!    write(*,*) p 
    DO
      CALL step_neargrid(bdr,chg,p)
!      IF (pout) write(*,'(3I4)') p
      ! if we didn't move, we're at a maximum
      IF (ALL(p == bdr%path(bdr%pnum,:))) EXIT
      ! otherwise, add point to path
      IF (bdr%pnum >= bdr%pdim) THEN
        CALL reallocate_path(bdr,bdr%pdim*2)
      ENDIF
      bdr%pnum = bdr%pnum+1
      bdr%path(bdr%pnum,:) = p
!
!GH: change this to be a known point (all neighbor points assigned)
      IF (opts%quit_opt==opts%quit_known .AND. bdr%known(p(1),p(2),p(3))==2) THEN
!        IF (pout) write(*,*) 'quitting at known point' 
        EXIT
      END IF
    END DO

  RETURN
  END SUBROUTINE max_neargrid

!-----------------------------------------------------------------------------------!
!  step_neargrid:  Do a single iteration of a maximization on the charge density 
!    grid from the point (px,py,pz).
!-----------------------------------------------------------------------------------!

  SUBROUTINE step_neargrid(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER, DIMENSION(3), INTENT(INOUT) :: p
    INTEGER, DIMENSION(3) :: pm
    REAL(q2), DIMENSION(3) :: gradrl, dr=(/0._q2,0._q2,0._q2/)
    REAL(q2) :: coeff
    SAVE dr

    IF (bdr%pnum == 1) THEN
      dr = (/0._q2,0._q2,0._q2/)
    END IF

    gradrl = rho_grad_dir(chg,p) 

    IF (MAXVAL(ABS(gradrl)) < 1E-30) THEN
      IF (is_max(chg,p)) THEN
        dr = (/0._q2,0._q2,0._q2/)
        RETURN
      ELSE
        pm = p
        CALL step_ongrid(chg,pm)
        dr = (/0._q2,0._q2,0._q2/)
      END IF 
    ELSE
      coeff = 1._q2/MAXVAL(ABS(gradrl))
      gradrl = coeff*gradrl
      pm = p + ANINT(gradrl)
      dr = dr + gradrl - ANINT(gradrl)
      pm = pm + ANINT(dr)
      dr = dr - ANINT(dr)
    END IF
    bdr%known(p(1),p(2),p(3)) = 1

    CALL pbc(pm,chg%npts)
    IF (bdr%known(pm(1),pm(2),pm(3)) == 1) THEN
      pm = p
      CALL step_ongrid(chg,pm)
      dr = (/0._q2,0._q2,0._q2/)
    END IF

    p = pm

  RETURN
  END SUBROUTINE step_neargrid

!-----------------------------------------------------------------------------------!
! refine_edge: refine the grid points on the edge of the Bader volumes.
!-----------------------------------------------------------------------------------!

  SUBROUTINE refine_edge(bdr,chg,opts,ions,ref_itrs)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    INTEGER :: ref_itrs

    INTEGER, DIMENSION(3) :: p,pt
    INTEGER :: n1,n2,n3,path_volnum,bvolnum,i
    INTEGER :: num_edge,num_reassign,num_check
    INTEGER :: d1,d2,d3

     IF(opts%refine_edge_itrs==-2 .OR. ref_itrs==1) THEN
       num_edge = 0
       DO n1 = 1,chg%npts(1)
        DO n2 = 1,chg%npts(2)
          DO n3 = 1,chg%npts(3)
            p = (/n1,n2,n3/)
            ! change for calculating the vacuum volume
            IF (bdr%volnum(n1,n2,n3) == bdr%bnum+1) CYCLE
            IF (is_vol_edge(bdr,chg,p) .AND. (.NOT.is_max(chg,p))) THEN
              num_edge = num_edge + 1
              bdr%volnum(p(1),p(2),p(3)) = -bdr%volnum(p(1),p(2),p(3))
              IF (opts%quit_opt == opts%quit_known) THEN
                bdr%known(p(1),p(2),p(3)) = 0
                CALL reassign_volnum_ongrid2(bdr,chg,p)
              END IF 
            END IF
          END DO
        END DO
      END DO
      WRITE(*,'(2x,A,6x,1I8)') 'EDGE POINTS:',num_edge
    END IF

    IF(opts%refine_edge_itrs==-1 .AND. ref_itrs>1) THEN
      num_check=0
      DO n1 = 1,chg%npts(1)
        DO n2 = 1,chg%npts(2)
          DO n3 = 1,chg%npts(3)
            p = (/n1,n2,n3/)
            ! change for calculating the vacuum volume
            IF (bdr%volnum(n1,n2,n3)==bdr%bnum+1) CYCLE

            IF(bdr%volnum(n1,n2,n3) < 0 .AND. bdr%known(n1,n2,n3) /=-1) THEN
              DO d1 = -1,1
               DO d2 = -1,1
                DO d3 = -1,1
                  pt = p + (/d1,d2,d3/)
                  CALL pbc(pt,chg%npts)
                  ! change for calculating the vacuum volume
                  IF (bdr%volnum(pt(1),pt(2),pt(3)) == bdr%bnum+1) CYCLE
                  IF(.NOT.is_max(chg,pt)) THEN 
                    IF(bdr%volnum(pt(1),pt(2),pt(3)) > 0) THEN
                      bdr%volnum(pt(1),pt(2),pt(3)) = -bdr%volnum(pt(1),pt(2),pt(3))
                      bdr%known(pt(1),pt(2),pt(3)) = -1
                      num_check=num_check+1
                    ELSE IF(bdr%volnum(pt(1),pt(2),pt(3))<0 .AND. bdr%known(pt(1),pt(2),pt(3)) == 0) THEN
                      bdr%known(pt(1),pt(2),pt(3)) = -2
                      num_check = num_check + 1
                    END IF
                  END IF
                END DO
               END DO
              END DO
              num_check = num_check - 1
              IF (bdr%known(pt(1),pt(2),pt(3)) /= -2) THEN
                bdr%volnum(p(1),p(2),p(3)) = ABS(bdr%volnum(p(1),p(2),p(3)))
              END IF
              ! end of mark
            END IF

          END DO
        END DO
      END DO
      WRITE(*,'(2x,A,3x,1I8)') 'CHECKED POINTS:', num_check

      ! make the surrounding points unkown
      DO n1 = 1,chg%npts(1)
        DO n2 = 1,chg%npts(2)
          DO n3 = 1,chg%npts(3)
            p = (/n1,n2,n3/)
            bvolnum = bdr%volnum(n1,n2,n3)

            IF (bvolnum < 0) THEN
              DO d1 = -1,1
               DO d2 = -1,1
                DO d3 = -1,1
                  pt = p + (/d1,d2,d3/)
                  CALL pbc(pt,chg%npts)
                  IF(bdr%known(pt(1),pt(2),pt(3)) == 2) bdr%known(pt(1),pt(2),pt(3)) = 0
                END DO
               END DO
              END DO
            END IF

          END DO
        END DO
      END DO

    END IF

    num_reassign = 0
    DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
          p = (/n1,n2,n3/)
          bvolnum = bdr%volnum(n1,n2,n3)
          IF (bvolnum<0) THEN
            IF (opts%bader_opt == opts%bader_offgrid) THEN
              CALL max_offgrid(bdr,chg,p)
            ELSEIF (opts%bader_opt == opts%bader_ongrid) THEN
              CALL max_ongrid(bdr,chg,p)
            ELSE
              CALL max_neargrid(bdr,chg,opts,p)
            END IF
            path_volnum = bdr%volnum(p(1),p(2),p(3))
            IF (path_volnum<0 .OR. path_volnum>bdr%bnum) THEN
              WRITE(*,*) 'ERROR: should be no new maxima in edge refinement'
            END IF
            bdr%volnum(n1,n2,n3) = path_volnum
            IF (ABS(bvolnum) /= path_volnum) THEN
              num_reassign = num_reassign + 1
              IF(opts%refine_edge_itrs==-1 .OR. opts%refine_edge_itrs==-3) bdr%volnum(n1,n2,n3) = -path_volnum
            END IF
            DO i = 1,bdr%pnum
              pt = (/bdr%path(i,1),bdr%path(i,2),bdr%path(i,3)/)
              IF (bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
                bdr%known(pt(1),pt(2),pt(3)) = 0
              END IF
            END DO 
          END IF
        END DO
      END DO
    END DO
 
    WRITE(*,'(2x,A,1I8)') 'REASSIGNED POINTS:',num_reassign

    ! flag to indicate that we are done refining
    IF ((opts%refine_edge_itrs==-1 .OR. opts%refine_edge_itrs==-3) .AND. num_reassign==0) THEN
      bdr%refine_edge_itrs = 0
    END IF
    IF (opts%refine_edge_itrs==-2 .AND. num_reassign==0) THEN
      bdr%refine_edge_itrs = 0
    END IF

  RETURN
  END SUBROUTINE refine_edge

!-----------------------------------------------------------------------------------!
! assign_chg2atom: Assign an element of charge to a Bader atom.
!-----------------------------------------------------------------------------------!
  SUBROUTINE assign_chg2atom(bdr,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2), DIMENSION(3) :: dv, v
    REAL(q2) :: dsq, dminsq
    INTEGER :: i, j, dindex

    bdr%ionchg = 0._q2
    DO i = 1,bdr%nvols
      dv = bdr%volpos_dir(i,:) - ions%r_dir(1,:)
!      CALL dpbc_dir(ions,dv)
      v = MATMUL(ions%dir2car, dv)
      CALL dpbc_car(ions, v)
      dminsq = DOT_PRODUCT(v, v)
      dindex = 1
      DO j = 2,ions%nions
        dv = bdr%volpos_dir(i,:) - ions%r_dir(j,:)
!        CALL dpbc_dir(ions,dv)
        v = MATMUL(ions%dir2car, dv)
        CALL dpbc_car(ions, v)
        dsq = DOT_PRODUCT(v, v)
        IF (dsq < dminsq) THEN
          dminsq = dsq
          dindex = j
        END IF
      END DO
      bdr%iondist(i) = SQRT(dminsq)
      bdr%nnion(i) = dindex
      bdr%ionchg(dindex) = bdr%ionchg(dindex) + bdr%volchg(i)
    END DO

  RETURN
  END SUBROUTINE assign_chg2atom

!-----------------------------------------------------------------------------------!
! bader_mindist: Find the minimum distance from the surface of each volume to 
!                each ion
!-----------------------------------------------------------------------------------!

  SUBROUTINE bader_mindist(bdr,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2), DIMENSION(3) :: v, dv_dir, dv_car
    INTEGER, DIMENSION(3) :: p
    REAL :: dist
    INTEGER :: i, atom, n1, n2, n3
    INTEGER :: cr, count_max, t1, t2, tenths_done
    
    CALL SYSTEM_CLOCK(t1, cr, count_max)

 
    WRITE(*,'(/,2x,A)') 'CALCULATING MINIMUM DISTANCES TO ATOMS'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

!   Store the minimum distance and the vector
    ALLOCATE(bdr%minsurfdist(ions%nions))
    bdr%minsurfdist = 0._q2

    tenths_done = 0
    DO n1 = 1,chg%npts(1)
      IF ((n1*10/chg%npts(1)) > tenths_done) THEN
        tenths_done = (n1*10/chg%npts(1))
        WRITE(*,'(A,$)') '**'
      END IF
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
          p = (/n1,n2,n3/)
          ! change for calculating the vacuum volume
          IF (bdr%volnum(n1,n2,n3) == bdr%nvols+1) CYCLE

!         If this is an edge cell, check if it is the closest to the atom so far
          IF(is_atm_edge(bdr,chg,p,atom)) THEN 
            v = REAL((/n1,n2,n3/),q2)
            dv_dir = (v-chg%org_lat)/REAL(chg%npts,q2) - ions%r_dir(atom,:)
!            CALL dpbc_dir(ions,dv_dir)
!            CALL dpbc_dir_org(dv_dir)
!            CALL matrix_vector(ions%dir2car, dv_dir, dv_car)
            dv_car = MATMUL(ions%dir2car, dv_dir)
            CALL dpbc_car(ions, dv_car)
            dist = DOT_PRODUCT(dv_car, dv_car)
            IF ((bdr%minsurfdist(atom) == 0.0_q2) .OR.  &
            &   (bdr%minsurfdist(atom) > dist)) THEN
              bdr%minsurfdist(atom) = dist
            END IF
          END IF
        END DO
      END DO
    END DO

    DO i = 1,ions%nions
      bdr%minsurfdist(i) = SQRT(bdr%minsurfdist(i))
    END DO

    CALL SYSTEM_CLOCK(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F7.2,1A8)') 'RUN TIME: ', (t2-t1)/REAL(cr,q2), ' SECONDS'

  RETURN
  END SUBROUTINE bader_mindist

!-----------------------------------------------------------------------------------!
! write_all_bader: Write out a CHGCAR type file for each of the Bader volumes found.
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_all_bader(bdr, opts, ions, chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp
    INTEGER :: atomnum, badercur, tenths_done, t1, t2, cr, count_max
    INTEGER :: n1, n2, n3
    CHARACTER(LEN=128) :: atomfilename
    
    CALL SYSTEM_CLOCK(t1, cr,count_max)

    WRITE(*,'(/,2x,A)') 'WRITING BADER VOLUMES'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    
    tmp = chg
    atomnum = 0
    tenths_done = 0

    DO badercur = 1,bdr%nvols
      DO WHILE ((badercur*10/bdr%nvols) > tenths_done)
        tenths_done = tenths_done+1
        WRITE(*,'(A,$)') '**'
      ENDDO
      IF(bdr%volchg(badercur) > bdr%tol) THEN
        atomnum = atomnum+1
        WRITE(atomfilename,'(A4,I4.4,A4)') "Bvol",atomnum,".dat"
        tmp%rho = 0._q2
!        WHERE (bdr%volnum == badercur) tmp%rho = chg%rho
        ! Replace for ifort/x86_32
        DO n1 = 1,chg%npts(1)
          DO n2 = 1,chg%npts(2)
            DO n3 = 1,chg%npts(3)
              IF (bdr%volnum(n1,n2,n3) == badercur) tmp%rho(n1,n2,n3) = chg%rho(n1,n2,n3)
            END DO
          END DO
        END DO
        CALL write_charge(ions,tmp,opts,atomfilename)
      END IF
    END DO

    DEALLOCATE(tmp%rho)

    CALL SYSTEM_CLOCK(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F7.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE write_all_bader

!-----------------------------------------------------------------------------------!
! write_all_atom: Write out a CHGCAR type file for each atom where all Bader volumes
!              assigned to that atom are added together. This is only done if the 
!              atoms has any 'significant' bader volumes associated with it.
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_all_atom(bdr,opts,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp

    INTEGER :: j, b, mab, mib, ik, sc, cc
    INTEGER :: tenths_done, t1, t2, cr, count_max
    INTEGER :: n1, n2, n3
    INTEGER, DIMENSION(bdr%nvols) :: rck
    CHARACTER(LEN=128) :: atomfilename

    CALL SYSTEM_CLOCK(t1, cr, count_max)

    tmp = chg

    WRITE(*,'(/,2x,A)') 'WRITING ATOMIC VOLUMES '
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    tenths_done = 0
    mab = MAXVAL(bdr%nnion)
    mib = MINVAL(bdr%nnion)
    sc = 0

    DO ik = mib,mab
      cc = 0
      rck = 0
      DO j = 1,bdr%nvols
        IF (bdr%volchg(j) < bdr%tol) CYCLE
        IF (bdr%nnion(j) == ik) THEN
          cc = cc+1
          rck(cc) = j
        END IF
      END DO
      sc = sc + cc
      DO WHILE ((ik*10/(mab-mib+1)) > tenths_done)
        tenths_done = tenths_done+1
        WRITE(*,'(A,$)') '**'
      END DO
      IF (cc == 0) CYCLE
      WRITE(atomfilename,'(A4,I4.4,A4)') "BvAt",ik,".dat"

      tmp%rho = 0._q2
      DO b=1,cc
!        WHERE (bdr%volnum == rck(b)) tmp%rho = chg%rho
        ! this change is needed for ifort/x86_32
        DO n1 = 1,chg%npts(1)
          DO n2 = 1,chg%npts(2)
            DO n3 = 1,chg%npts(3)
              IF (bdr%volnum(n1,n2,n3) == rck(b)) tmp%rho(n1,n2,n3) = chg%rho(n1,n2,n3)
            END DO
          END DO
        END DO
      END DO 
      CALL write_charge(ions, tmp, opts, atomfilename)

    END DO
    DEALLOCATE(tmp%rho)

    CALL SYSTEM_CLOCK(t2, cr, count_max)
    WRITE(*,'(/,1A12,1F7.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE write_all_atom

!-----------------------------------------------------------------------------------!
! write_sel_atom: Write a CHGCAR file for the user specified atomic volumes.
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_sel_atom(bdr,opts,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp

    INTEGER :: n1, n2, n3, i, j, b, ik, sc, cc
    INTEGER :: tenths_done, t1, t2, cr, count_max
    INTEGER, DIMENSION(bdr%nvols) :: rck
    CHARACTER(LEN=128) :: atomfilename

    CALL SYSTEM_CLOCK(t1,cr,count_max)

    tmp=chg

    WRITE(*,'(/,2x,A)') 'WRITING SELECTED ATOMIC VOLUMES '
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    tenths_done = 0
    sc = 0

    DO i = 1,opts%selanum
      ik = opts%selavol(i)
      cc = 0
      rck = 0
      DO j = 1,bdr%nvols
        IF (bdr%volchg(j) < bdr%tol) CYCLE
        IF (bdr%nnion(j) == ik) THEN
          cc = cc+1
          rck(cc) = j
        END IF
      END DO
      sc = sc + cc
      DO WHILE ((i*10/opts%selanum) > tenths_done)
        tenths_done = tenths_done+1
        WRITE(*,'(A,$)') '**'
      END DO
      IF (cc == 0) CYCLE
      WRITE(atomfilename,'(A4,I4.4,A4)') "BvAt",ik,".dat"

      tmp%rho = 0._q2
      DO b = 1,cc
!        WHERE (bdr%volnum == rck(b)) tmp%rho = chg%rho
        DO n1 = 1,chg%npts(1)
          DO n2 = 1,chg%npts(2)
            DO n3 = 1,chg%npts(3)
              IF (bdr%volnum(n1,n2,n3) == rck(b))  &
              &   tmp%rho(n1,n2,n3) = chg%rho(n1,n2,n3)
            END DO
          END DO
        END DO

      END DO 
      CALL write_charge(ions,tmp,opts,atomfilename)

    END DO
    DEALLOCATE(tmp%rho)

    CALL SYSTEM_CLOCK(t2, cr, count_max)
    WRITE(*,'(/,1A12,1F7.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE write_sel_atom

!-----------------------------------------------------------------------------------!
! write_sum_atom: Write a CHGCAR file for the charge summed over the user specified 
! atomic volumes.
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_sum_atom(bdr,opts,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp

    INTEGER :: n1, n2, n3, i, j, b, ik, sc, cc
    INTEGER :: tenths_done, t1, t2, cr, count_max
    INTEGER, DIMENSION(bdr%nvols) :: rck
    CHARACTER(LEN=128) :: atomfilename

    CALL SYSTEM_CLOCK(t1, cr, count_max)

    tmp = chg
!   Set temporary density to zero at start
    tmp%rho = 0._q2

    WRITE(*,'(/,2x,A)') 'WRITING SUM OVER SELECTED ATOMIC VOLUMES '
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    tenths_done = 0
    sc = 0

!   Here we only need name one file
    WRITE(atomfilename,'(A)') "BvAt_summed.dat"

!   Loop over number of selected atoms (opts%selanum)
    DO i = 1,opts%sumanum
      ik = opts%sumavol(i)
      cc = 0
      rck = 0
      DO j = 1,bdr%nvols
        IF (bdr%volchg(j) < bdr%tol) CYCLE
        IF (bdr%nnion(j) == ik) THEN
          cc = cc+1
          rck(cc) = j
        END IF
      END DO
      sc = sc + cc
      DO WHILE ((i*10/opts%sumanum) > tenths_done)
        tenths_done = tenths_done+1
        WRITE(*,'(A,$)') '**'
      END DO
      IF (cc == 0) CYCLE
      DO b = 1,cc
!        WHERE (bdr%volnum == rck(b)) tmp%rho = chg%rho
        DO n1 = 1,chg%npts(1)
          DO n2 = 1,chg%npts(2)
            DO n3 = 1,chg%npts(3)
              IF (bdr%volnum(n1,n2,n3) == rck(b))  &
              &   tmp%rho(n1,n2,n3) = chg%rho(n1,n2,n3)
            END DO
          END DO
        END DO
      END DO

    END DO

!   Write out charge to file
    CALL write_charge(ions, tmp, opts, atomfilename)
    DEALLOCATE(tmp%rho)

    CALL SYSTEM_CLOCK(t2, cr, count_max)
    WRITE(*,'(/,1A12,1F7.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE write_sum_atom

!-----------------------------------------------------------------------------------!
! write_sel_bader: Write a CHGCAR file for the user specified Bader volumes.
!              Volumes associated with a atom can be read from AtomVolumes.dat
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_sel_bader(bdr, opts, ions, chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp
    CHARACTER(LEN=128) :: atomfilename
    INTEGER, DIMENSION(bdr%nvols,2) :: volsig
    INTEGER :: n1, n2, n3, cr, count_max, t1, t2
    INTEGER :: tenths_done, i, bdimsig, bvolnum

    CALL SYSTEM_CLOCK(t1, cr, count_max)

    tmp = chg

! Correlate the number for each significant bader volume to its real number
    bdimsig = 0
    volsig = 0

    DO i = 1,bdr%nvols
      IF (bdr%volchg(i) > bdr%tol) THEN
        bdimsig = bdimsig+1
        volsig(bdimsig,1) = bdimsig
        volsig(bdimsig,2) = i
      END IF
    END DO

    WRITE(*,'(/,2x,A)') 'WRITING SELECTED BADER VOLUMES '
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    tenths_done = 0

    tmp%rho = 0._q2

    DO i = 1,opts%selbnum
      DO WHILE ((i*10/opts%selbnum) > tenths_done)
        tenths_done = tenths_done+1
        WRITE(*,'(A,$)') '**'
      END DO

      bvolnum = opts%selbvol(i)
      WRITE(atomfilename,'(A4,I4.4,A4)') "Bvol",bvolnum,".dat"
      tmp%rho = 0._q2
!      WHERE (bdr%volnum == volsig(bvolnum,2)) tmp%rho = chg%rho
      DO n1 = 1,chg%npts(1)
        DO n2 = 1,chg%npts(2)
          DO n3 = 1,chg%npts(3)
            IF (bdr%volnum(n1,n2,n3) == volsig(bvolnum,2))  &
            &  tmp%rho(n1,n2,n3) = chg%rho(n1,n2,n3)
          END DO
        END DO
      END DO

      CALL write_charge(ions, tmp, opts, atomfilename)
    END DO

    DEALLOCATE(tmp%rho)

    CALL SYSTEM_CLOCK(t2, cr, count_max)
    WRITE(*,'(/,1A12,1F7.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE write_sel_bader

!-----------------------------------------------------------------------------------!
! write_sum_bader: Write a CHGCAR file for the charge summed over the user specified 
! bader volumes.
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_sum_bader(bdr, opts, ions, chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp
    CHARACTER(LEN=128) :: atomfilename
    INTEGER, DIMENSION(bdr%nvols,2) :: volsig
    INTEGER :: n1, n2, n3, cr, count_max, t1, t2
    INTEGER :: tenths_done, i, bdimsig, bvolnum

    CALL SYSTEM_CLOCK(t1, cr, count_max)

    tmp = chg

! Correlate the number for each significant bader volume to its real number
    bdimsig = 0
    volsig = 0 

    DO i = 1,bdr%nvols
      IF (bdr%volchg(i) > bdr%tol) THEN
        bdimsig = bdimsig + 1
        volsig(bdimsig,1) = bdimsig
        volsig(bdimsig,2) = i
      END IF
    END DO
    
    WRITE(*,'(/,2x,A)') 'WRITING SELECTED BADER VOLUMES '
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    tenths_done = 0

    tmp%rho = 0._q2
    WRITE(atomfilename,'(A15)') "Bvol_summed.dat"

    DO i = 1,opts%sumbnum
      DO WHILE ((i*10/opts%sumbnum) > tenths_done)
        tenths_done = tenths_done + 1
        WRITE(*,'(A,$)') '**'
      END DO

      bvolnum = opts%sumbvol(i)
!      WHERE (bdr%volnum == volsig(bvolnum,2)) tmp%rho = chg%rho
      DO n1 = 1,chg%npts(1)
        DO n2 = 1,chg%npts(2)
          DO n3 = 1,chg%npts(3)
            IF (bdr%volnum(n1,n2,n3) == volsig(bvolnum,2))  &
            &  tmp%rho(n1,n2,n3) = chg%rho(n1,n2,n3)
          END DO
        END DO
      END DO
    END DO
      
    CALL write_charge(ions, tmp, opts, atomfilename)
    DEALLOCATE(tmp%rho)
    
    CALL SYSTEM_CLOCK(t2, cr, count_max)
    WRITE(*,'(/,1A12,1F7.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'
  
  RETURN
  END SUBROUTINE write_sum_bader


!-----------------------------------------------------------------------------------!
! write_bader_index: Write out a CHGCAR type file that the value at each point is
!                     Bader volume number for that grid point.
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_bader_index(bdr, opts, ions, chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp
    CHARACTER(LEN=128) :: filename
    REAL(q2):: vol
    INTEGER :: n1, n2, n3

    WRITE(*,'(/,2x,A)') 'WRITING BADER VOLUME INDEX TO BvIndex.dat'

    tmp = chg
!GH: to be compatible with ifort/x86_32
!    tmp%rho=bdr%volnum
    DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
          tmp%rho(n1,n2,n3) = bdr%volnum(n1,n2,n3)
        END DO
      END DO
    END DO

    filename = 'BvIndex.dat'
    IF(opts%out_opt == opts%out_cube) THEN
      vol = matrix_volume(ions%lattice)
!GH: to be compatible with ifort/x86_32
!      tmp%rho=tmp%rho*vol
      DO n1 = 1,chg%npts(1)
        DO n2 = 1,chg%npts(2)
          DO n3 = 1,chg%npts(3)
            tmp%rho(n1,n2,n3) = tmp%rho(n1,n2,n3)*vol
          END DO
        END DO
      END DO

    END IF

    CALL write_charge(ions, tmp, opts, filename)
    DEALLOCATE(tmp%rho)

  END SUBROUTINE write_bader_index

!-----------------------------------------------------------------------------------!
! write_atom_index: Write out a CHGCAR type file that the value at each point is
!                    associated atom number for that grid point.
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_atom_index(bdr, opts, ions, chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp
    CHARACTER(LEN=128) :: filename
    REAL(q2):: vol
    INTEGER :: n1, n2, n3

    WRITE(*,'(/,2x,A)') 'WRITING BADER ATOMIC INDEX TO AtIndex.dat'

    tmp = chg
    DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
          ! change this to accord with the vacuum volume calculation
          IF(bdr%volnum(n1,n2,n3) == bdr%bnum+1) THEN
            tmp%rho(n1,n2,n3) = ions%nions+1
          ELSE
            tmp%rho(n1,n2,n3) = bdr%nnion(bdr%volnum(n1,n2,n3))
          END IF
        END DO
      END DO
    END DO

    filename = 'AtIndex.dat'
    IF(opts%out_opt == opts%out_cube) THEN
      vol = matrix_volume(ions%lattice)
!      tmp%rho=tmp%rho*vol
      DO n1 = 1,chg%npts(1)
        DO n2 = 1,chg%npts(2)
          DO n3 = 1,chg%npts(3)
            tmp%rho(n1,n2,n3) = tmp%rho(n1,n2,n3)*vol
          END DO
        END DO
      END DO

    END IF

    CALL write_charge(ions, tmp, opts, filename)
    DEALLOCATE(tmp%rho)

  RETURN
  END SUBROUTINE write_atom_index

!-----------------------------------------------------------------------------------!
! write_atom_index: Write out a file that contains a list with grid indices that
!                   belong to a surface as well the atom indices.
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_surfaces_atoms(bdrscell, opts, ionsscell, chg, chgscell, sc_atom_map)
    IMPLICIT NONE

    TYPE(bader_obj) :: bdrscell
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ionsscell
    TYPE(charge_obj) :: chg
    TYPE(charge_obj) :: chgscell
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: sc_atom_map(:)
    INTEGER :: nx, ny, nz, atom, atmnbr, d1, d2, d3, nn, i
    INTEGER, DIMENSION(3) :: p, pt
    INTEGER, DIMENSION(27) :: neighbours
    LOGICAL :: is_present
 
    IF ( .NOT. ALLOCATED(sc_atom_map) ) THEN
      ALLOCATE(sc_atom_map(ionsscell%nions))
      DO i=1,ionsscell%nions
        sc_atom_map(i) = i
      ENDDO
    END IF

    WRITE(*,'(/,2x,A)') 'WRITING SURFACE TO FILE BaderSurfaces.dat'
    
    OPEN( 42,FILE="BaderSurfaces.dat",ACTION='write')
    WRITE( 42, *) "# List of point in grid coordinates ix, iy and iz marking the border between atoms"
    ! Even though obtained in the supercell, we only look at 
    ! surfaces in the primitive cell
    DO nz=1,chg%npts(3)
      DO ny=1,chg%npts(2)
        DO nx=1,chg%npts(1)
          p = (/ nx, ny, nz /)
          atom = bdrscell%nnion(bdrscell%volnum(p(1),p(2),p(3)))
          IF (bdrscell%volnum(p(1),p(2),p(3)) == bdrscell%bnum + 1) THEN
             CYCLE
          END IF
          nn = 1
          neighbours = 0
          DO d1 = -1,1
            DO d2 = -1,1
              DO d3 = -1,1
                pt = p + (/d1,d2,d3/)
                CALL pbc(pt,chgscell%npts)
                IF (bdrscell%volnum(pt(1),pt(2),pt(3)) == bdrscell%bnum + 1) THEN
                  ! Vaccum defines a surface, too, we encode it with -1
                  is_present = .FALSE.
                  DO i = 1,nn
                    IF ( neighbours(i) == -1 ) is_present = .TRUE.
                  END DO
                  IF ( .not. is_present ) THEN
                    neighbours(nn) = -1
                    nn = nn + 1
                  END IF
                END IF
                atmnbr = bdrscell%nnion(bdrscell%volnum(pt(1),pt(2),pt(3)))
                IF (atmnbr /= atom) THEN
                  is_present = .FALSE.
                  DO i = 1,nn
                    IF ( neighbours(i) == sc_atom_map(atmnbr) ) is_present = .TRUE.
                  END DO
                  IF ( .not. is_present ) THEN
                    neighbours(nn) = sc_atom_map(atmnbr)
                    nn = nn + 1
                  END IF
                END IF
              END DO
            END DO
          END DO
          IF ( ANY(neighbours(:) /= 0 ) ) THEN
            WRITE( 42, *) p(1),p(2),p(3),sc_atom_map(atom),neighbours(1:nn-1)
          ENDIF
        END DO
      END DO
    END DO
    CLOSE(42)
  END SUBROUTINE write_surfaces_atoms

!-----------------------------------------------------------------------------------!
! cal_atomic_vol: Integrate the atomic volume for each atom
!-----------------------------------------------------------------------------------!

  SUBROUTINE cal_atomic_vol(bdr,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2):: vol
    INTEGER :: n1, n2, n3, i, atom

    ALLOCATE(bdr%ionvol(ions%nions))
    bdr%ionvol = 0.0

    DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
          ! change to calculate vacuum volume
          IF (bdr%volnum(n1,n2,n3) /= bdr%nvols + 1) THEN
            atom = bdr%nnion(bdr%volnum(n1,n2,n3))
            bdr%ionvol(atom) = bdr%ionvol(atom) + 1
          END IF
        END DO
      END DO
    END DO

    vol = matrix_volume(ions%lattice)
    vol = vol/chg%nrho

    DO i = 1,ions%nions
      bdr%ionvol(i) = bdr%ionvol(i)*vol
    END DO

  RETURN
  END SUBROUTINE cal_atomic_vol

!-----------------------------------------------------------------------------------!
! bader_output: Write out a summary of the bader analysis.
!         AtomVolumes.dat: Stores the 'significant' Bader volumes associated with
!                          each atom.
!         ACF.dat        : Stores the main output to the screen.
!         BCF.dat        : Stores 'significant' Bader volumes, their coordinates and
!                          charge, atom associated and distance to it. 
!-----------------------------------------------------------------------------------!

  SUBROUTINE bader_output(bdr,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2) :: sum_ionchg, ne
    INTEGER :: i, bdimsig, mib, mab, cc, j, nmax
    INTEGER, DIMENSION(bdr%nvols) :: rck

    mab=MAXVAL(bdr%nnion)
    mib=MINVAL(bdr%nnion)
    OPEN(100,FILE='AVF.dat',STATUS='replace',ACTION='write')
    WRITE(100,'(A)') '   Atom                     Volume(s)'
    WRITE(100,'(A,A)') '  --------------------------------------------------------',&
    &                  '----------------'
    DO i = mib, mab
      cc = 0
      rck = 0
      nmax = 0
      DO j = 1,bdr%nvols
        IF (bdr%volchg(j) > bdr%tol) THEN
          nmax = nmax+1
          IF (bdr%nnion(j) == i) THEN
            cc = cc+1
            rck(cc) = nmax
          END IF
        END IF
      END DO 
      IF (cc == 0) CYCLE
      WRITE(100,'(2X,1I4,2X,A,2X,10000I5)') i,' ... ',rck(1:cc)
    END DO
    WRITE(100,'(A,A)') '  --------------------------------------------------------',&
    &                  '----------------'
    CLOSE(100)

    WRITE(*,'(/,A41)') 'WRITING BADER ATOMIC CHARGES TO ACF.dat'
    WRITE(*,'(A41,/)') 'WRITING BADER VOLUME CHARGES TO BCF.dat'

    OPEN(100,FILE='ACF.dat',STATUS='replace',ACTION='write')
    WRITE(100,555) '#','X','Y','Z','CHARGE','MIN DIST','ATOMIC VOL'
    555 FORMAT(4X,1A,9X,1A1,2(11X,1A1),7X,1A6,6X,1A8,3X,1A10)
    WRITE(100,'(A,A)') ' ----------------------------------------------------------------',&
    &                  '----------------'
    sum_ionchg = SUM(bdr%ionchg)
    ne = SUM(bdr%volchg(1:bdr%nvols)) + bdr%vacchg
    DO i = 1, ions%nions
!      WRITE(100,'(1I5,7F12.7)') i,ions%r_car(i,:),bdr%ionchg(i),bdr%minsurfdist(i),bdr%ionvol(i)
      WRITE(100,'(1I5,4F12.6,2F13.6)') i,ions%r_car(i,:),bdr%ionchg(i),bdr%minsurfdist(i),bdr%ionvol(i)
    END DO
    WRITE(100,'(A,A)') ' ----------------------------------------------------------------',&
    &                  '----------------'
    WRITE(100,'(4X,A,9X,1F12.4)') 'VACUUM CHARGE:',bdr%vacchg
    WRITE(100,'(4X,A,9X,1F12.4)') 'VACUUM VOLUME:',bdr%vacvol
    WRITE(100,'(4X,A,2X,1F12.4)')  'NUMBER OF ELECTRONS: ',ne
    CLOSE(100)

    bdimsig = 0
    OPEN(200,FILE='BCF.dat',STATUS='replace',ACTION='write')

    WRITE(200,556) '#','X','Y','Z','CHARGE','ATOM','DISTANCE'
    556 FORMAT(4X,1A1,9X,1A1,2(11X,1A1),8X,1A6,5X,1A4,4X,1A8)

    WRITE(200,'(A,A)') '  ----------------------------------------------------------',&
    &              '---------------'
    DO i = 1,bdr%nvols
        IF (bdr%volchg(i) > bdr%tol) THEN
           bdimsig = bdimsig+1
           WRITE(200,777) bdimsig,bdr%volpos_car(i,:),bdr%volchg(i), &
           &              bdr%nnion(i),bdr%iondist(i)
           777 FORMAT(1I5,4F12.4,3X,1I5,1F12.4)
        END IF
    END DO
    WRITE(200,'(A,A)') '  ----------------------------------------------------------',&
    &              '---------------'
    CLOSE(200)

    WRITE(*,'(2x,A,6X,1I8)')       'NUMBER OF BADER MAXIMA FOUND: ', bdr%nvols
    WRITE(*,'(2x,A,6X,1I8)')       '    SIGNIFICANT MAXIMA FOUND: ', bdimsig
    WRITE(*,'(2X,A,2X,1F12.4)')    '               VACUUM CHARGE: ', bdr%vacchg
    WRITE(*,'(2x,A,2X,1F12.5)')    '         NUMBER OF ELECTRONS: ', ne

  RETURN
  END SUBROUTINE bader_output

!-----------------------------------------------------------------------------------!
!  rho_val:  Return the density at the point (p1,p2,p3) taking into account the
!    boundary conditions.  This function is used to address points outside the
!    charge density array without a bunch of if statements at the place the value
!    is needed.
!-----------------------------------------------------------------------------------!

  FUNCTION volnum_val(bdr,chg,p1,p2,p3)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER, INTENT(IN) :: p1, p2, p3
    INTEGER :: i, volnum_val
    INTEGER, DIMENSION(3) :: p

    p = (/p1,p2,p3/)
    DO i = 1,3
      DO
        IF (p(i) >= 1) EXIT
        p(i) = p(i) + chg%npts(i)
      END DO
      DO
        IF (p(i) <= chg%npts(i)) EXIT
        p(i) = p(i) - chg%npts(i)
      END DO
    END DO

    volnum_val = bdr%volnum(p(1),p(2),p(3))

  RETURN
  END FUNCTION volnum_val

!-----------------------------------------------------------------------------------!
! known_volnum: return number of the associated bader volnum if all surrounding
!    grid points are known to be associated with the same bader volnum
!-----------------------------------------------------------------------------------!

  FUNCTION known_volnum(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    REAL(q2), DIMENSION(3), INTENT(IN) :: p
    INTEGER :: known_volnum

    INTEGER :: volnum, d1, d2, d3, p1, p2, p3
    LOGICAL :: first_flag

    known_volnum = 0
    first_flag = .TRUE.

    p1 = FLOOR(p(1))
    p2 = FLOOR(p(2))
    p3 = FLOOR(p(3))

    DO d1 = 0,1
      p1 = p(1)+d1
      DO d2 = 0,1
        p2 = p(2)+d2
        DO d3 = 0,1
          p3 = p(3)+d3
          IF (first_flag) THEN
            volnum = volnum_val(bdr,chg,p1,p2,p3)
            IF (volnum <= 0) RETURN
            first_flag = .FALSE.
          ELSE
            IF (volnum /= volnum_val(bdr,chg,p1,p2,p3)) RETURN
          END IF
        END DO
      END DO
    END DO
    known_volnum = volnum

  RETURN
  END FUNCTION known_volnum

!-----------------------------------------------------------------------------------!
! assign_surrounding_pts: check the surrounding points of p to see if their volnum
!                         is known
!-----------------------------------------------------------------------------------!
  SUBROUTINE assign_surrounding_pts(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER, DIMENSION(3), INTENT(IN) :: p
    INTEGER, DIMENSION(3) :: pt

    pt = p + (/1,0,0/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN 
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF
    pt = p + (/-1,0,0/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF
    pt = p + (/0,1,0/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF
    pt = p + (/0,-1,0/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF
    pt = p + (/0,0,1/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF
    pt = p + (/0,0,-1/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF

  RETURN
  END SUBROUTINE assign_surrounding_pts

!-----------------------------------------------------------------------------------!
! assign_surrounding_pts: check the surrounding points of p to see if their volnum
!                         is known
!-----------------------------------------------------------------------------------!

  SUBROUTINE assign_surrounding_pts2(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER, DIMENSION(3), INTENT(IN) :: p
    INTEGER, DIMENSION(3) :: pt
    INTEGER :: d1, d2, d3

    DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p + (/d1,d2,d3/)
          CALL pbc(pt,chg%npts)
          IF (bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
            CALL known_volnum_ongrid2(bdr,chg,pt)
          END IF
        END DO
      END DO
    END DO 

  RETURN
  END SUBROUTINE assign_surrounding_pts2

!-----------------------------------------------------------------------------------!
! known_volnum_ongrid: return number of the associated bader volnum if nearest
!    grid points are known to be associated with the same bader volnum
!-----------------------------------------------------------------------------------!

  SUBROUTINE known_volnum_ongrid(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER, DIMENSION(3), INTENT(IN) :: p
    INTEGER :: volnum, p1, p2, p3

    p1 = p(1)
    p2 = p(2)
    p3 = p(3)     

    volnum = volnum_val(bdr,chg,p1,p2,p3)
    IF(volnum <= 0) RETURN

    IF (volnum_val(bdr,chg,p1,p2,p3+1) /= volnum) RETURN
    IF (volnum_val(bdr,chg,p1,p2,p3-1) /= volnum) RETURN
    IF (volnum_val(bdr,chg,p1,p2+1,p3) /= volnum) RETURN
    IF (volnum_val(bdr,chg,p1,p2-1,p3) /= volnum) RETURN
    IF (volnum_val(bdr,chg,p1+1,p2,p3) /= volnum) RETURN
    IF (volnum_val(bdr,chg,p1-1,p2,p3) /= volnum) RETURN

    bdr%known(p1,p2,p3) = 2

  RETURN
  END SUBROUTINE known_volnum_ongrid

!-----------------------------------------------------------------------------------!
! known_volnum_ongrid: return number of the associated bader volnum if nearest
!    grid points are known to be associated with the same bader volnum
!-----------------------------------------------------------------------------------!

  SUBROUTINE known_volnum_ongrid2(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER, DIMENSION(3), INTENT(IN) :: p
    INTEGER, DIMENSION(3) :: pt
    INTEGER :: volnum, d1, d2, d3, p1, p2, p3

    p1 = p(1)
    p2 = p(2)
    p3 = p(3)

    volnum = volnum_val(bdr,chg,p1,p2,p3)
    IF (volnum <= 0) RETURN

    DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p + (/d1,d2,d3/)
          IF (volnum_val(bdr,chg,pt(1),pt(2),pt(3)) /= volnum) RETURN
        END DO
      END DO
    END DO
    bdr%known(p1,p2,p3) = 2

  RETURN
  END SUBROUTINE known_volnum_ongrid2

!-----------------------------------------------------------------------------------!
! reassign_volnum_ongrid: reassign the surrounding points of a edge point as unknown
!                         points
!-----------------------------------------------------------------------------------!

  SUBROUTINE reassign_volnum_ongrid(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER, DIMENSION(3), INTENT(IN) :: p
    INTEGER, DIMENSION(3) :: pt

    pt = (/p(1)+1,p(2),p(3)/)
    CALL pbc(pt,chg%npts)
    bdr%known(pt(1),pt(2),pt(3)) = 0
    pt = (/p(1)-1,p(2),p(3)/)
    CALL pbc(pt,chg%npts)
    bdr%known(pt(1),pt(2),pt(3)) = 0
    pt = (/p(1),p(2)+1,p(3)/)
    CALL pbc(pt,chg%npts)
    bdr%known(pt(1),pt(2),pt(3)) = 0
    pt = (/p(1),p(2)-1,p(3)/)
    CALL pbc(pt,chg%npts)
    bdr%known(pt(1),pt(2),pt(3)) = 0
    pt = (/p(1),p(2),p(3)+1/)
    CALL pbc(pt,chg%npts)
    bdr%known(pt(1),pt(2),pt(3)) = 0
    pt = (/p(1),p(2),p(3)-1/)
    CALL pbc(pt,chg%npts)
    bdr%known(pt(1),pt(2),pt(3)) = 0

  RETURN
  END SUBROUTINE reassign_volnum_ongrid

!-----------------------------------------------------------------------------------!
! reassign_volnum_ongrid: reassign the surrounding points of a edge point as unknown
!                         points
!-----------------------------------------------------------------------------------!

  SUBROUTINE reassign_volnum_ongrid2(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER, DIMENSION(3), INTENT(IN) :: p
    INTEGER, DIMENSION(3) :: pt
    INTEGER :: d1, d2, d3

    DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p + (/d1,d2,d3/)
          CALL pbc(pt,chg%npts)
          bdr%known(pt(1),pt(2),pt(3)) = 0
        END DO
      END DO
    END DO

  RETURN
  END SUBROUTINE reassign_volnum_ongrid2

!-----------------------------------------------------------------------------------!
! is_vol_edge: return .true. if the grid point is on the edge of a Bader volume.
!-----------------------------------------------------------------------------------!

  FUNCTION is_vol_edge(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    LOGICAL :: is_vol_edge

    INTEGER, DIMENSION(3), INTENT(IN) :: p
    INTEGER, DIMENSION(3) ::pt
    INTEGER :: d1,d2,d3,volnum,volnbr

    volnum = bdr%volnum(p(1),p(2),p(3))
    is_vol_edge = .FALSE.
    neighborloop: DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p + (/d1,d2,d3/)
          CALL pbc(pt,chg%npts)
          volnbr = bdr%volnum(pt(1),pt(2),pt(3))
          IF (ABS(volnbr) /= ABS(volnum)) THEN
            is_vol_edge = .TRUE.
            EXIT neighborloop  
          END IF
        END DO
      END DO
    END DO neighborloop

  RETURN
  END FUNCTION is_vol_edge

!-----------------------------------------------------------------------------------!
! is_atm_edge: return .true. if the grid point is on the edge of a Bader atom
!-----------------------------------------------------------------------------------!

  FUNCTION is_atm_edge(bdr,chg,p,atom)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    LOGICAL :: is_atm_edge

    INTEGER, DIMENSION(3), INTENT(IN) :: p
    INTEGER, DIMENSION(3) ::pt
    INTEGER :: d1, d2, d3, atmnbr
    INTEGER, INTENT(INOUT) ::atom 
    
    atom = bdr%nnion(bdr%volnum(p(1),p(2),p(3)))
    is_atm_edge = .FALSE.
    neighborloop: DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p + (/d1,d2,d3/)
          CALL pbc(pt,chg%npts)
          IF (bdr%volnum(pt(1),pt(2),pt(3)) == bdr%bnum + 1) CYCLE  ! vacuum point
          atmnbr = bdr%nnion(bdr%volnum(pt(1),pt(2),pt(3)))
          IF (atmnbr /= atom) THEN
            is_atm_edge = .TRUE.
            EXIT neighborloop
          END IF
        END DO
      END DO
    END DO neighborloop
    
    RETURN
    END FUNCTION is_atm_edge

!-----------------------------------------------------------------------------------!
! is_vol_neighbor: return .true. if the grid point neighboring the bader volume
!-----------------------------------------------------------------------------------!

  FUNCTION is_vol_neighbor(bdr,chg,p,vol)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    LOGICAL :: is_vol_neighbor

    INTEGER, DIMENSION(3), INTENT(IN) :: p
    INTEGER, DIMENSION(3) :: pt
    INTEGER :: d1, d2, d3, volneighbor
    INTEGER, INTENT(IN) :: vol

    is_vol_neighbor = .FALSE.

    ! only find neighbors, not points in the volume
    IF (bdr%volnum(pt(1),pt(2),pt(3)) == vol) RETURN

    neighborloop: DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p + (/d1,d2,d3/)
          CALL pbc(pt,chg%npts)
          volneighbor = bdr%volnum(pt(1),pt(2),pt(3))
          IF (volneighbor == vol) THEN
            is_vol_neighbor = .TRUE.
            EXIT neighborloop
          END IF
        END DO
      END DO
    END DO neighborloop

    RETURN
    END FUNCTION is_vol_neighbor

!-----------------------------------------------------------------------------------!
! reallocate_volpos: 
!-----------------------------------------------------------------------------------!

  SUBROUTINE reallocate_volpos(bdr,newsize)

    TYPE(bader_obj) :: bdr
    INTEGER :: newsize

    !STC: tmpvolpos was INTEGER should have been REAL(q2), this caused a crash
    !STC: when bader was run with -ref and compiled with floating point
    !STC: exception handling (-fpe0).
    REAL(q2), ALLOCATABLE, DIMENSION(:,:) :: tmpvolpos

    IF (newsize < bdr%bnum) write(*,*) 'Error: new volpos length too small'

    ALLOCATE(tmpvolpos(bdr%bdim,3))
    tmpvolpos = bdr%volpos_lat
    DEALLOCATE(bdr%volpos_lat)
    bdr%bdim = newsize
    ALLOCATE(bdr%volpos_lat(bdr%bdim,3))
    bdr%volpos_lat(1:bdr%bnum,:) = tmpvolpos(1:bdr%bnum,:)
    DEALLOCATE(tmpvolpos)

  END SUBROUTINE reallocate_volpos

!-----------------------------------------------------------------------------------!
! reallocate_path: 
!-----------------------------------------------------------------------------------!

  SUBROUTINE reallocate_path(bdr,newsize)

    TYPE(bader_obj) :: bdr
    INTEGER :: newsize

    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tmppath

    IF (newsize < bdr%pnum) write(*,*) 'Error: new path length too small'

    ALLOCATE(tmppath(bdr%pdim,3))
    tmppath = bdr%path
    DEALLOCATE(bdr%path)
    bdr%pdim = newsize
    ALLOCATE(bdr%path(bdr%pdim,3))
    bdr%path(1:bdr%pnum,:) = tmppath(1:bdr%pnum,:)
    DEALLOCATE(tmppath)

  END SUBROUTINE reallocate_path
!-----------------------------------------------------------------------------------!


!-----------------------------------------------------------------------------------!
! bader_check_partitioning
!-----------------------------------------------------------------------------------!

  SUBROUTINE bader_check_partitioning(bdr,chgval,scell_dir_out)

    IMPLICIT NONE
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgval
    LOGICAL, DIMENSION(3), intent(out) :: scell_dir_out

    INTEGER :: atom, n1, n2, n3
    INTEGER, DIMENSION(3) :: p
    INTEGER, ALLOCATABLE :: strings(:,:)
    ! Check if there is a path taking us from the x=0 to x=max plane without 
    ! crossing a surface in which case we need a supercell ... 
    ALLOCATE( strings(chgval%npts(2),chgval%npts(3)) )
    strings(:,:) = 1
    DO n2=1,chgval%npts(2)
      DO n3=1,chgval%npts(3)
        DO n1=1,chgval%npts(1)
          p = (/n1,n2,n3/)
          atom = bdr%nnion(bdr%volnum(n1,n2,n3))
          IF ( is_atm_edge(bdr,chgval,p,atom) ) THEN
            strings(n2,n3) = 0
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    scell_dir(1) = ANY(strings(:,:) == 1)
    DEALLOCATE(strings)

    ! Same thing for y and z ...
    ALLOCATE( strings(chgval%npts(1),chgval%npts(3)) )
    strings(:,:) = 1
    DO n1=1,chgval%npts(1)
      DO n3=1,chgval%npts(3)
        DO n2=1,chgval%npts(2)
          p = (/n1,n2,n3/)
          atom = bdr%nnion(bdr%volnum(n1,n2,n3))
          IF ( is_atm_edge(bdr,chgval,p,atom) ) THEN
            strings(n1,n3) = 0
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    scell_dir(2) = ANY(strings(:,:) == 1)
    DEALLOCATE(strings)

    ALLOCATE( strings(chgval%npts(1),chgval%npts(2)) )
    strings(:,:) = 1
    DO n1=1,chgval%npts(1)
      DO n2=1,chgval%npts(2)
        DO n3=1,chgval%npts(3)
          p = (/n1,n2,n3/)
          atom = bdr%nnion(bdr%volnum(n1,n2,n3))
          IF ( is_atm_edge(bdr,chgval,p,atom) ) THEN
            strings(n1,n2) = 0
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  
    scell_dir(3) = ANY(strings(:,:) == 1)
    DEALLOCATE(strings)

    scell_dir_out = scell_dir
    
  END SUBROUTINE bader_check_partitioning

END MODULE bader_mod
