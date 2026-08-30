  MODULE critpoint_mod
    USE kind_mod
    USE matrix_mod
    USE bader_mod
!    USE linsolve_mod
    USE chgcar_mod
    USE charge_mod 
    USE options_mod
    USE ions_mod
    USE io_mod
    USE ions_mod
    USE weight_mod
    USE dsyevj3_mod
    USE jacobi_mod
    IMPLICIT NONE

    PRIVATE 
    PUBLIC :: critpoint_find, StaticCheck

    TYPE static_cp_list
      CHARACTER(LEN=1) :: cp_type
      REAL(q2),DIMENSION(3) :: pos_direct
    END TYPE

    TYPE cpc ! stands for critical point candidate
      INTEGER, DIMENSION(3) :: ind  ! these are the indices of the cp
      REAL(q2), DIMENSION(3) :: trueind,truer,grad,eigvals,r
      REAL(q2), DIMENSION(3,3) :: hessianMatrix,eigvecs
      !REAL(q2), DIMENSION(3) :: tempcart, tempind
      REAL(q2) :: rho
      INTEGER, DIMENSION(2) :: connectedAtoms
      INTEGER :: negcount, weight
      LOGICAL :: hasProxy, isunique
    END TYPE
    REAL(q2) :: voxvol
    CONTAINS
  
    ! Counts the amount of negative modes in eigenvalues
    ! USED IN THIS MODULE
    FUNCTION CountNegModes(eigvals)
      REAL(q2),DIMENSION(3) :: eigvals
      INTEGER :: CountNegModes,i
      CountNegModes = 0
      DO i = 1,3
        IF (eigvals(i) < 0) CountNegModes = CountNegModes + 1
      END DO
      RETURN
    END FUNCTION 

!-----------------------------------------------------------------------------------!
!critpoint_find: find critical points on the edge of the Bader volumes
!NOTE: this subroutine should be called after refine_edge
!      in order to restrict the calculation to edge points
!-----------------------------------------------------------------------------------!

  SUBROUTINE GetCPCL(bdr,chg,cpl,cpcl,opts,cptnum)
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts
    TYPE(cpc), ALLOCATABLE, DIMENSION(:) :: cpl,cpcl,cpclt

    REAL(q2), DIMENSION(3,3) :: hessianMatrix
    REAL(q2), DIMENSION(3) :: tem,trueR,grad

    INTEGER, DIMENSION(3) :: p
    INTEGER :: n1,n2,n3,cptnum,i
 

    OUTER: DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
            ! check to see if this point is in the vacuum
          IF (bdr%volnum(n1,n2,n3) == bdr%bnum + 1) THEN
            CYCLE
          END IF
          ! We can only work with points that are on edge as defined by bader
          ! but that has been shown to be not reliable, as critical points are
          ! missed commonly. 
          p = (/n1,n2,n3/)
          trueR = (/REAL(n1,q2),REAL(n2,q2),REAL(n3,q2)/)
          tem = CalcTEMGrid(p,chg,grad,hessianMatrix)
          !IF ( (ABS(tem(1)) <= 1.5 + opts%par_tem .AND. &
          !     ABS(tem(2)) <= 1.5 + opts%par_tem .AND. &
          !     ABS(tem(3)) <= 1.5 + opts%par_tem )) THEN
          IF (ALL(tem <= 1.5 + opts%par_tem )) THEN
              ! ABS(tem(3)) <= 0.5 + opts%par_tem) .OR. &
              !(SUM(grad*grad) <= (0.1*opts%par_gradfloor)**2 )) THEN              
            ! finding proximity could potentially be costly
            IF (ProxyToCPCandidate(p,opts,cpcl,cptnum,chg)) THEN
              CYCLE
            END IF
            cptnum = cptnum + 1
            ! Check if the candidate list needs to be expanded.
            IF (cptnum < SIZE(cpcl) - 1 ) THEN
              cpcl(cptnum)%ind = (/n1,n2,n3/)
              cpcl(cptnum)%grad = grad
              cpcl(cptnum)%hasProxy = .FALSE.
              cpcl(cptnum)%r = tem
            ELSE 
              IF (cptnum > 100000) THEN
                PRINT *, "ERROR: Too many searches required. Aborting."
                EXIT OUTER
              END IF
              PRINT *, 'expanding cpcl size'
              ALLOCATE(cpclt(cptnum + 1))
              DO i = 1, cptnum - 1
                cpclt(i) = cpcl(i)
              END DO
              DEALLOCATE(cpcl)
              ALLOCATE(cpcl(cptnum*2))
              DO i = 1, cptnum - 1
                cpcl(i)=cpclt(i)
              END DO
              DEALLOCATE(cpclt)
              ALLOCATE(cpclt(cptnum + 1))
              DO i = 1, cptnum - 1
                cpclt(i) = cpl(i)
              END DO
              DEALLOCATE(cpl)
              ALLOCATE(cpl(cptnum*2))
              DO i = 1, cptnum - 1
                cpl(i)=cpclt(i)
              END DO
              DEALLOCATE(cpclt)
              cpcl(cptnum)%ind = (/n1,n2,n3/)
              cpcl(cptnum)%grad = grad
              cpcl(cptnum)%hasProxy = .FALSE.
              cpcl(cptnum)%r = tem
            END IF
          END IF
        END DO
      END DO
    END DO OUTER
  END SUBROUTINE GetCPCL

  SUBROUTINE SearchWithCPCL(bdr,chg,cpcl,cpl,cptnum,ucptnum,ucpCounts,opts)
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(cpc), ALLOCATABLE, DIMENSION(:) :: cpcl,cpl
    TYPE(options_obj) :: opts

    REAL(q2), DIMENSION(8,3,3) :: nnHes
    REAL(q2), DIMENSION(3,3) :: interpolHessian
    REAL(q2),DIMENSION(3) :: temcap,temscale,trueR,distance
    REAL(q2) :: temnormcap

    INTEGER, DIMENSION(4) :: ucpCounts
    INTEGER, DIMENSION(2) :: connectedAtoms
    INTEGER :: i,cptnum,ucptnum       
    DO i = 1, cptnum
      cpcl(i)%isunique = .FALSE.
      temcap = (/1.,1.,1./)
      temscale = (/1.,1.,1./)
      temnormcap = 1.
      IF (opts%gradMode) THEN
         !uses GradientDescend instead of NRTFGP          
         CALL GradientDescend(bdr,chg,opts,trueR,cpcl(i)%ind,&
         cpcl(i)%isUnique,3000)
      ELSE
         ! Begins newton raphson validation process
         CALL NRTFGP(bdr,chg,opts,trueR,&
         cpcl(i)%isUnique,cpcl(i)%r,cpcl(i)%ind,&
         1000)
      END IF
      IF (cpcl(i)%isUnique ) THEN
        !CALL MakeCPRoster(cpRoster,i,truer) !not sure what this does
        cpcl(i)%trueind = trueR
        interpolHessian = trilinear_interpol_hes(nnHes,distance)
        ucptnum = ucptnum + 1
        interpolHessian = CDHessianR(truer,chg)
        CALL RecordCPR(truer,chg,cpl,ucptnum,connectedAtoms, ucpCounts, &
          opts, interpolHessian, &
          cpcl(i)%ind)
        CYCLE
      ELSE
        CYCLE
      END IF
      CALL pbc_r_lat(trueR,chg%npts)
    ! moving on to the next critical pint candidate
    END DO
  END SUBROUTINE SearchWithCPCL

  ! This subroutine reads in a list of CPs and their types, runs it through ReduceCP and PHRuleExam
  SUBROUTINE StaticCheck(bdr,chg,opts,ions)
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cp_static
    TYPE(ions_obj) :: ions
    TYPE(options_obj) :: opts
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: atom_connectivity
    INTEGER,DIMENSION(4) :: ucpCounts
    INTEGER :: n, ucptnum, ij, it_num, rot_num
    LOGICAL :: phmrCompliant, isReduced

    isReduced = .FALSE.
    ALLOCATE(cp_static(ReadStaticSize()))
    !ALLOCATE(atom_connectivity(ions%nions,ions%nions))
    !atom_connectivity = 0
    CALL StaticCheckReadStatic(cp_static,chg)
    DO n = 1, SIZE(cp_static)
      cp_static(n)%truer(1) = cp_static(n)%truer(1) * chg%npts(1) + 1
      cp_static(n)%truer(2) = cp_static(n)%truer(2) * chg%npts(1) + 1
      cp_static(n)%truer(3) = cp_static(n)%truer(3) * chg%npts(1) + 1
    END DO
    ucptnum = SIZE(cp_static)
    ucpCounts = 0
    DO WHILE(.NOT. isReduced)
      CALL ReduceCPStatic(cp_static, ucptnum, ucpCounts, isReduced, opts)
    END DO
    PRINT *, "the cp set after reduction is "
    PRINT *, "cp number | signature | positions "
    DO n = 1, ucptnum
      PRINT *, n, (3 - cp_static(n)%negCount) - cp_static(n)%negCount , cp_static(n)%truer
    END DO
    PRINT *, "nuclear | bond | ring | cage cp numbers are"
    PRINT *, ucpCounts
    DO ij = 1,ucptnum
      !Saves pairs of connected atoms in connectedAtoms
      cp_static(ij)%hessianMatrix = CDHessianR(cp_static(ij)%trueR,chg)
      CALL jacobi_eigenvalue(3,cp_static(ij)%hessianMatrix,9999,&
        cp_static(ij)%eigvecs, cp_static(ij)%eigvals,&
        it_num, rot_num )
      IF (cp_static(ij)%negcount == 2) THEN
        cp_static(ij)%connectedAtoms = DoubleAscension(bdr,cp_static(ij),chg,ions,opts,1,cp_static,ucptnum)
        !atom_connectivity(cp_static(ij)%connectedAtoms(1),cp_static(ij)%connectedAtoms(2))=1
        !atom_connectivity(cp_static(ij)%connectedAtoms(2),cp_static(ij)%connectedAtoms(1))=1
      ELSEIF (cp_static(ij)%negcount == 1) THEN
        cp_static(ij)%connectedAtoms = DoubleAscension(bdr,cp_static(ij),chg,ions,opts,-1,cp_static,ucptnum)
      END IF
    END DO
    CALL PHRuleExam(ucpCounts,opts,ions,phmrCompliant)
    IF (phmrCompliant .AND. .NOT. CheckIsolatedAtom(cp_static,ucptnum)) THEN
      PRINT *, 'The CPs are self-consistent but isolated atoms are detected. &
        Declaring this a false positive.'
    END IF      
    IF (phmrCompliant .AND. .NOT. CheckIsolatedRing(cp_static,ucptnum)) THEN
      PRINT *, "A disconnected ring CP is found. Declaring this a false positive."
    END IF
    IF (phmrCompliant .AND. .NOT. CheckIsolatedCage(cp_static,ucptnum)) THEN
      PRINT *, "A disconnected cage CP is found. Declaring this a false positive."
    END IF
    DEALLOCATE(cp_static)
    !DEALLOCATE(atom_connectivity)

  END SUBROUTINE StaticCheck

  SUBROUTINE critpoint_find(bdr,chg,opts,ions,stat)
! These are for screening CP due to numerical error. 
    !TYPE(hessian) :: hes
    TYPE(bader_obj), INTENT(INOUT) :: bdr
    TYPE(charge_obj), INTENT(INOUT) :: chg
    TYPE(options_obj), INTENT(INOUT) :: opts
    TYPE(ions_obj), INTENT(INOUT) :: ions
    TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpcl, cpl, cpclt
    TYPE(static_cp_list), ALLOCATABLE, DIMENSION(:) :: cps_read

    ! The above three are CP candidate list, CP list and CP list temp
    ! copy
    ! for points, 1 and 2 are +1, -1

    REAL(q2), DIMENSION(:,:), ALLOCATABLE :: cpRoster, fullcpRoster, reducedcpRoster, &
      static_search
    REAL(q2), DIMENSION(8,3,3) :: nnhes !hessian of 8 nn
    REAL(q2), DIMENSION(10,3) :: rList,temList
    REAL(q2), DIMENSION(8,3) :: nngrad  
    REAL(q2), DIMENSION(3,3) :: hessianMatrix, eigvecs, interpolHessian, &
      ggrid
    REAL(q2), DIMENSION(3) :: tem, eigvals, truer, grad, prevgrad, temprealr, &
      distance, & ! vector to 000 in trilinear 
      finR, nexttem, previoustem, averager, temcap, temscale
    REAL(q2) :: temnormcap

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: descendPoints, ringPoints, nnind, &
      atom_connectivity, nucleiInd, ring_connectivity
    INTEGER, DIMENSION(:), ALLOCATABLE :: RingList
    INTEGER, DIMENSION(20,3) :: iniIList
    INTEGER, DIMENSION(4) :: cpCounts,ucpCounts ! the 4 elements are
    ! nuclear, bond, ring, cage CP coutns
    INTEGER, DIMENSION(3) :: p, tempr
    INTEGER, DIMENSION(2) :: connectedAtoms
    INTEGER :: n1, n2, n3, cptnum, ucptnum, i, j, k, debugnum,ij, &
      setcount, stat, axisnum,&
      avgMode, stepcount, nnlayers, averagecount

    CHARACTER(128) :: smoothenedCHGCAR

    LOGICAL :: HCF, LDM, LDM_DetectCircling, LDM_ReduceCP, LDM_DensityDescend,LDM_RecordCPRLight, &
      LDM_NRTFGP,LDM_CalcTEMLat,LDM_RecordCPR, LDM_GradMagGrad, LDM_RingAscend, LDM_Trajectories, &
      isUniqueTest, isReduced, phmrCompliant, isReducible

    IF (opts%static_check) THEN
      CALL StaticCheck(bdr,chg,opts,ions)
    ELSE 
      !CALL PrintFlavorText()
      ! below are variables for least sqaures gradient
      stat = 0 ! 0 means nothing


      !WRITE(*,'(A)')  'FINDING CRITICAL POINTS'
      CALL get_voxvol(chg,ions)
      IF (opts%enableCHGCARSmoothening) THEN
        ! First produce a smoothened CHGCAR
        CALL SmoothenCHGCAR(chg,avgMode,ions,opts)
        ! This will override chg object
        smoothenedCHGCAR = "smoothenedCHGCAR_sum"
        ! Deallocating objects allocated in read_charge_chgcar
        DEALLOCATE (ions%num_ion,ions%r_dir,ions%r_car,chg%rho,ions%r_lat)
        CALL read_charge_chgcar(ions,chg,smoothenedCHGCAR,opts)
      END IF
      HCF = .FALSE.
      LDM = .FALSE.
      IF (opts%debugMode) THEN
        !PRINT *, "Performing core functions check"
        !CALL CoreFunctionsCheck(chg)
        PRINT *, "Reading debug flags from debugConfig"
        CALL GetDebugFlags(opts,LDM,LDM_DetectCircling,&
          LDM_ReduceCP,LDM_DensityDescend,LDM_RecordCPRLight,&
          LDM_NRTFGP,LDM_CalcTEMLat,LDM_RecordCPR,LDM_GradMagGrad,LDM_RingAscend,LDM_Trajectories)
      ELSE
        LDM_DetectCircling = .FALSE.
        LDM_ReduceCP = .FALSE.
        LDM_DensityDescend = .FALSE.
        LDM_RecordCPRLight = .FALSE.
        LDM_NRTFGP = .FALSE.
        LDM_CalcTEMLat = .FALSE.
        LDM_RecordCPR = .FALSE.
        LDM_GradMagGrad = .FALSE.
        LDM_RingAscend = .FALSE.
        LDM_Trajectories = .FALSE.
      END IF
      ! get expected nucleus indices
!      WRITE (97,*) , 'expecting nuclei at:'
      ALLOCATE(nnInd(8,3))
!      ALLOCATE(nnInd((nnlayers*2+1)**3,3))
      ALLOCATE(nucleiInd(ions%nions,3))
      DO n1 = 1, ions%nions
        nucleiInd(n1,1) = NINT(ions%r_lat(n1,1))
        nucleiInd(n1,2) = NINT(ions%r_lat(n1,2))
        nucleiInd(n1,3) = NINT(ions%r_lat(n1,3))
!        WRITE (97,*) nucleiInd(d1,:)
      END DO
      setcount = 0
      cpCounts = 0
      ucpCounts = 0
      ucptnum = 0
      cptnum = 0
      debugnum = 0
      nnlayers = findnnlayers(ions)
      nnind = findnn(truer,nnlayers,chg,ions)
      PRINT *, '******************************************************' 
      ALLOCATE (cpl(10000)) ! start with 10000 capacity
      ALLOCATE (cpcl(10000))
      IF (opts%static_search) THEN
        ALLOCATE(cps_read(ReadStaticSize()))
        CALL ReadStatic(cps_read)
        DO i = 1, SIZE(cps_read)
          cpcl(i)%ind=INT(cps_read(i)%pos_direct(:)*chg%npts(:))+(/1,1,1/)
          PRINT *, "Starting search at assigned index: "
          PRINT *, cpcl(i)%ind
          CALL NRTFGP(bdr,chg,opts,trueR,cpcl(i)%isUnique,cpcl(i)%r,cpcl(i)%ind,&
          1000)
          PRINT *, "Search trajectory ended at "
          PRINT *, trueR
          cpcl(i)%trueind = trueR
          interpolHessian = trilinear_interpol_hes(nnHes,distance)
          interpolHessian = CDHessianR(truer,chg)
          ucptnum = ucptnum + 1
          CALL RecordCPR(trueR,chg,cpl,ucptnum,connectedAtoms, ucpCounts, &
            opts, interpolHessian, &
            cpcl(i)%ind)
          PRINT *, "cpl ucptnum trueR is ", cpl(ucptnum)%trueR(:)
          cpl(ucptnum)%isUnique = .TRUE.
          PRINT *, "cpl ucptnum negcount is ", cpl(ucptnum)%negCount
        END DO
        CALL VisAllCP(cpl,ucptnum,chg,ions,opts,ucpCounts)
        PRINT *, "called visallcp"
        STOP
      END IF
        ! Maxima should always be found first with gradient ascension
      !gradient descent unit test with arbitrary initial point
      !iniI and finR were defined in the main critpoint_fnd method for lack of a better place to define them.
      IF (LDM_GradMagGrad) THEN
        iniIList(6,:) = (/47,47,52/)
        axisnum = 3
        DO ij=1,11
            IF (axisnum == 1) THEN
               iniIList(ij,:) = (/iniIList(6,1)+ij-6,iniIList(6,2), iniIList(6,3)/)
            ELSE IF (axisnum == 2) THEN
               iniIList(ij,:) = (/iniIList(6,1),iniIList(6,2)+ij-6, iniIList(6,3)/)
            ELSE
               iniIList(ij,:) = (/iniIList(6,1),iniIList(6,2), iniIList(6,3)+ij-6/)
            END IF
        END DO
        iniIList(12,:) = (/43, 44, 49/)
        iniIList(13,:) = (/43, 45, 49/)
        DO ij=1,13
           ! PRINT *, "Segfault check loop",ij
            isUniqueTest = .TRUE.
            CALL GradientDescend(bdr,chg,opts,finR,iniIList(ij,:), isUniqueTest ,3000)
            PRINT *, finR
            PRINT *, " "
        END DO
      END IF

      IF (LDM) PRINT *, 'Finding maxima near atoms'
      DO n1 = 1, ions%nions
        tempr = NINT(ions%r_dir(n1,:) * chg%npts)
        temprealr = ascension(tempr,chg, &
                  ggrid,opts,ions)
        ucptnum = ucptnum + 1
        cpl(ucptnum)%trueind = temprealr
        ! these points are to stay in the list so cpl is used instead of cpcl
        p = (/NINT(cpl(ucptnum)%trueind(1)),NINT(cpl(ucptnum)%trueind(2)), &
          NINT(cpl(ucptnum)%trueind(3))/)
        hessianMatrix = CDHessian(p,chg)
        CALL RecordCPR(temprealr,chg,cpl,ucptnum,connectedAtoms, ucpCounts, &
          opts,hessianMatrix, p)
        IF (LDM) THEN
          PRINT *, "recording CP to cpl, this is number ", ucptnum
          PRINT *, "location is "
          PRINT *, temprealr
        END IF
      END DO
      ! ascension results may not give the best hessian. so the types are set
      ! manually.
      IF (ucpCounts(1) /= ucptnum) THEN
        !PRINT *, "ucpCounts(1) and ucptnum are ", ucpCounts(1), ucptnum
        PRINT *, 'WARNING: It was detected that the number of maxima found &
          does not equal to trials started. The found critical points are &
          being manually set as nuclear critical points'
        DO n1 = 1, ucptnum
          cpl(n1)%negcount = 3
        END DO
        ucpCounts(1) = ucptnum
      END IF

      IF (opts%static_search) THEN
        ALLOCATE(cps_read(ReadStaticSize()))
        CALL ReadStatic(cps_read)
        DO i = 1, SIZE(cps_read)
          cpcl(i)%ind=INT(cps_read(i)%pos_direct(:)*chg%npts(:))+(/1,1,1/)
          PRINT *, "Starting search at assigned index: "
          PRINT *, cpcl(i)%ind
          CALL NRTFGP(bdr,chg,opts,trueR,cpcl(i)%isUnique,cpcl(i)%r,cpcl(i)%ind,&
          1000)
          PRINT *, "Search trajectory ended at "
          PRINT *, trueR
          cpcl(i)%trueind = trueR
          interpolHessian = trilinear_interpol_hes(nnHes,distance)
          interpolHessian = CDHessianR(truer,chg)
          ucptnum = ucptnum + 1
          CALL RecordCPR(trueR,chg,cpl,ucptnum,connectedAtoms, ucpCounts, &
            opts, interpolHessian, &
            cpcl(i)%ind)
          PRINT *, "cpl ucptnum trueR is ", cpl(ucptnum)%trueR(:)
          cpl(ucptnum)%isUnique = .TRUE.
          PRINT *, "cpl ucptnum negcount is ", cpl(ucptnum)%negCount
        END DO
        CALL VisAllCP(cpl,ucptnum,chg,ions,opts,ucpCounts)
        PRINT *, "called visallcp"
        STOP
      ELSE 
        ! Loop through every grid point once and collect a list of points to start
        ! CP searching trajectories into cpcl, the CP candidate list.
        CALL GetCPCL(bdr,chg,cpl,cpcl,opts,cptnum)
        IF (cptnum > 100000) THEN
          stat = 0
        ELSE 
          PRINT *, "Number of Newton Rhapson trajectory needed: ", cptnum 
!!**      *****************************************************************
          ! To find critical points (unique), start with a cell that contains a
          ! critical point and its hessian and force. Use Newton's method to make a
          ! move. Interpolate the force inside the voxel. 
          ! Once moved, get the new force through trilinear interpolation, and
          ! get the new hessian which will be a matrix of constants, make moves until
          ! r is zero. get the coordinates of the new true critical point. If this
          ! point is within half lattice to another, do not record this new point.
          ALLOCATE(cpRoster(cptnum,3))
          IF (LDM_Trajectories) ALLOCATE(fullcpRoster(cptnum,3))
          CALL SearchWithCPCL(bdr,chg,cpcl,cpl,cptnum,ucptnum,ucpCounts,opts)

          PRINT *, 'Number of critical point count: ', ucptnum
          PRINT *, 'Number of nuclear, bond, ring and cage  critical point &
            counts : ', ucpCounts(:)
          ! remove duplicate CPs
          isReduced = .FALSE.
          isReducible = .TRUE.
          DO WHILE ( .NOT. isReduced .AND. isReducible)
            CALL ReduceCP(cpl,opts,ucptnum,chg,ucpCounts, &
              isReduced,LDM_RecordCPRLight,LDM_ReduceCP, isReducible)
          END DO
          
          ! After CP numbers are reduced, do density descend and reduce another
          ! round afterwards
          IF (opts%enableDensityDescend .AND. isReducible) THEN
            PRINT *, "starting density descend"
            ! First find the midpoint of all rings, and their mirror iMages across
            ! the closest pbc (to be implemented). 
            ALLOCATE(descendPoints(ucpCounts(3)*(ucpCounts(3) - 1 ),3))
            ALLOCATE(ringPoints(ucpCounts(3),3))
            CALL ResizeCPL(cpl,SIZE(cpl) + ucpCounts(3)*(ucpCounts(3) - 1))
            PRINT *, "finished resizecpl"
            j = 1
            DO i = 1, ucptnum
              IF (cpl(i)%negCount == 1 ) THEN
                ringPoints(j,:) = cpl(i)%truer 
                j = j + 1
              END IF
            END DO
            PRINT *, "finished getting ringPoints"
            k = 1
            DO i = 1, ucpCounts(3)
              DO j = i+1, ucpCounts(3)
                descendPoints(k,1) = INT((ringPoints(i,1) + ringPoints(j,1))/2)
                descendPoints(k,2) = INT((ringPoints(i,2) + ringPoints(j,2))/2)
                descendPoints(k,3) = INT((ringPoints(i,3) + ringPoints(j,3))/2)
                k = k + 1
              END DO
            END DO
            DO i = 1, ( ucpCounts(3) * ( ucpCounts(3) - 1) ) / 2
              CALL DensityDescendAndRecord(chg,bdr,opts,descendPoints(i,:),cpl,&
                ucptnum, ucpCounts, .FALSE.)
            END DO
            DEALLOCATE(descendPoints)
            DEALLOCATE(ringPoints)
            isReduced = .FALSE.
            isReducible = .TRUE.
            DO WHILE ( .NOT. isReduced .AND. isReducible)
              CALL ReduceCP(cpl,opts,ucptnum,chg,ucpCounts, &
                isReduced,LDM_RecordCPRLight, &
                LDM_ReduceCP,isReducible)
            END DO
          END IF
          ! This following debug line need to be togged on or off manually
          ! before compilling
          !ip = (/-0.554,1.953,1.985/)
          !CALL CPTracer(iP,chg,cpl,ucptnum)
          PRINT *, 'After a round of reduction'
          PRINT *, 'Number of atoms: ', ions%nions
          PRINT *, 'Number of critical point count: ', ucptnum
          PRINT *, 'Number of nuclear, bond, ring and cage  critical point &
            counts : ', ucpCounts(:)
          CALL PHRuleExam(ucpCounts,opts,ions,phmrCompliant)
          IF (phmrCompliant) THEN
            stat = 1
          ELSE
            stat = 0
          END IF
          ! stat = 1
          !Runs DoubleAscension and RingAscension on all detected critical points
          !ALLOCATE(atom_connectivity(ucptnum,ucptnum))
          !atom_connectivity = 0
          DO ij = 1,ucptnum
            !Saves pairs of connected atoms in connectedAtoms
            IF (cpl(ij)%negcount == 2) THEN
              cpl(ij)%connectedAtoms = DoubleAscension(bdr,cpl(ij),chg,ions,opts,1,cpl,ucptnum)
              !atom_connectivity(cpl(ij)%connectedAtoms(1),cpl(ij)%connectedAtoms(2))=1
              !atom_connectivity(cpl(ij)%connectedAtoms(2),cpl(ij)%connectedAtoms(1))=1
            ELSEIF (cpl(ij)%negcount == 1) THEN
              cpl(ij)%connectedAtoms = DoubleAscension(bdr,cpl(ij),chg,ions,opts,-1,cpl,ucptnum)
            END IF
           
            IF (LDM_RingAscend) THEN
               !Prints Rings
               CALL RingAscension(cpl(ij),chg,ions,RingList)
            END IF
          END DO
          ! output the cp to files
          CALL OutputCP(cpl,opts,ucptnum,chg,setcount, ucpCounts)
          !Output the found list of connectivity pairs to file 
          CALL OutputNetwork(cpl,ucptnum,setcount)
          IF (stat == 1) THEN
            IF ( .NOT. CheckIsolatedAtom(cpl,ucptnum)) THEN
              PRINT *, 'The CPs are self-consistent but isolated atoms are detected.'
              stat = 0
            END IF      
            IF (.NOT. CheckIsolatedRing(cpl,ucptnum)) THEN
              PRINT *, "A disconnected ring CP is found."
              stat = 0
            END IF
            IF (.NOT. CheckIsolatedCage(cpl,ucptnum)) THEN
              PRINT *, "A disconnected cage CP is found."
              stat = 0
            END IF
            IF (stat == 0) THEN 
              PRINT *, 'Declaring this a false positive.'
            END IF
          END IF


          IF (LDM_Trajectories) THEN
            !Performs statistical analysis (standard deviation) on the CP Roster 
            CALL CPRosterAnalysis(cpl,ions,fullcpRoster,chg)
            !unique_realcoords removes NaN/0 valued blank CPs, condenses the list to only the relevant ones inside reducedcpRoster
            CALL unique_realcoords(cpRoster,reducedcpRoster) 
            !Output the CP Roster to file
            CALL OutputCPRoster(fullcpRoster,setcount)
            DEALLOCATE(fullcpRoster)
            DEALLOCATE(reducedcpRoster)
          END IF
          DEALLOCATE(cpRoster)
        END IF
      END IF
      IF (cptnum > 100000) THEN
        PRINT *, ''//achar(27)//'[31m ERROR: FAILED Poincare Hopf Rule & 
          and Morse Relationship' //achar(27)//'[0m'      
      ELSE
        PRINT *, 'outputting debugging information to allcpPOSCAR'
        CALL VisAllCP(cpl,ucptnum,chg,ions,opts,ucpCounts)
      END IF
      DEALLOCATE(cpl)
      DEALLOCATE(cpcl)
      !PRINT *, "DEALLOCATING atom_connectivity"
      !DEALLOCATE(atom_connectivity)
      !DEALLOCATE(cpl,cpcl,atom_connectivity)
    END IF
    PRINT *, "end of critpoint_find"
  END SUBROUTINE critpoint_find

    ! this function determins when looking for nn, how many layers to search
    ! within. It looks for the smallest vector sum of lattice vectors, and the
    ! largest vector
    ! USED IN THIS MODULE
    FUNCTION findnnlayers(ions)
      INTEGER :: findnnlayers
      TYPE(ions_obj) :: ions
      REAL(q2), DIMENSION(4,3) :: latsums!lat12, lat13, lat23, lat123
      REAL(q2) :: latMag, minMag, maxMag
      INTEGER :: i
      latsums(1,:) = ions%lattice(1,:) + ions%lattice(2,:)
      latsums(2,:) = ions%lattice(1,:) + ions%lattice(3,:)
      latsums(3,:) = ions%lattice(2,:) + ions%lattice(3,:)
      latsums(4,:) = ions%lattice(1,:) + ions%lattice(2,:) + ions%lattice(3,:)
      minMag = Mag(latsums(1,:))
      maxMag = Mag(latsums(1,:))
      DO i = 2, 4
        latMag = Mag(latsums(i,:))
        IF (latMag >= maxMag) THEN
          maxMag = latMag
        END IF
        IF (latMag <= minMag) THEN
          minMag = latMag
        END IF
      END DO
      DO i = 1, 3
        latMag = Mag(ions%lattice(i,:))
        IF (latMag >= maxMag) THEN
          maxMag = latMag
        END IF
        IF (latMag <= minMag) THEN
          minMag = latMag
        END IF
      END DO
      findnnlayers = CEILING(maxMag/minMag)

    END FUNCTION findnnlayers
  
    ! This function gives the simple box for trilinear interpolation.
    ! USED IN THIS MODULE
    FUNCTION SimpleNN(p,chg)
      REAL(q2), DIMENSION(3) :: p
      TYPE(charge_obj) :: chg
      INTEGER, DIMENSION(8,3) :: SimpleNN
      SimpleNN(1,:) = (/FLOOR(p(1)),FLOOR(p(2)),FLOOR(p(3))/)
      SimpleNN(2,:) = (/CEILING(p(1)),FLOOR(p(2)),FLOOR(p(3))/)
      SimpleNN(3,:) = (/FLOOR(p(1)),CEILING(p(2)),FLOOR(p(3))/)
      SimpleNN(4,:) = (/CEILING(p(1)),CEILING(p(2)),FLOOR(p(3))/)
      SimpleNN(5,:) = (/floor(p(1)),FLOOR(p(2)),CEILING(p(3))/)
      SimpleNN(6,:) = (/ceiling(p(1)),FLOOR(p(2)),CEILING(p(3))/)
      SimpleNN(7,:) = (/floor(p(1)),CEILING(p(2)),CEILING(p(3))/)
      SimpleNN(8,:) = (/ceiling(p(1)),CEILING(p(2)),CEILING(p(3))/)
!      DO i = 1, 8
!        CALL pbc(SimpleNN(i,:),chg%npts)
!      END DO
    END FUNCTION


    ! find neares on grid points for a interpolated point
    ! to do trilinear interpolation
    ! p is the off grid point we want
    FUNCTION FindNN(truer,nnlayers,chg,ions) ! THIS FUNCTION IS NOT STABLE
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nnind
      !grid points
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      TYPE(weight_obj), ALLOCATABLE, DIMENSION(:) :: nndist 
      ! instead of storing weight, this is used to
      ! store distances so that they can be sorted
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: FindNN
      REAL(q2),DIMENSION(3) :: pcart
      REAL(q2), DIMENSION(3) :: truer,truecart
      INTEGER,DIMENSION(3) :: baseind
      INTEGER, ALLOCATABLE, DIMENSION(:) :: rank
      INTEGER :: i,j,k,counter, nnlayers
      ! since only closest neigherbos of current grid point
      ! will be used for interpolation
      counter = 0
      nnlayers = FindNNLayers(ions)
      ALLOCATE(nndist((nnlayers*2 + 1)**3 ))
      ALLOCATE(nnind((nnlayers*2 + 1)**3 ,3))
      ALLOCATE(rank((nnlayers*2 + 1)**3 ))
      ALLOCATE(FindNN((nnlayers*2+1)**3,3))
      !to calculate distance it is not necessary to run pbc
      !infact pbc should be avoided at this stage
      baseind = (/FLOOR(truer(1)),FLOOR(truer(2)),FLOOR(truer(3))/)
      truecart = MATMUL(chg%lat2car,truer)
      DO i = -nnlayers, nnlayers
        DO j = -nnlayers, nnlayers
          DO k = -nnlayers, nnlayers
            !IF (i == 0 .AND. j == 0 .AND. k == 0) THEN
            !  CYCLE
            !END IF
            counter = counter + 1
            nnind(counter,:) = (/i,j,k/) + baseind
            pcart = MATMUL(chg%lat2car,nnind(counter,:))
            nndist(counter)%rho = Mag(truecart - pcart)
            nndist(counter)%pos(:) = (/i,j,k/) + baseind
            ! remember, rho here is really the distance
            FindNN(counter,:) = (/i,j,k/) + baseind
          END DO
        END DO
      END DO
      DEALLOCATE(nnind)
      DEALLOCATE(nndist)
      RETURN
    END FUNCTION FindNN

    ! This function should take in a list of nearest neighbors predetermined.
    FUNCTION nn_grad(chg,r,rho,nn)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2),INTENT(OUT) :: rho
      REAL(q2),DIMENSION(3) :: nn_grad
      INTEGER :: p1, p2, p3
      REAL(q2),DIMENSION(3) :: rho_grad_lat
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110, rho111
      REAL(q2) :: rho00_, rho01_, rho10_, rho11_
      REAL(q2) :: rho0__, rho1__, rho_0_, rho_1_, rho__0, rho__1
      REAL(q2) :: rho_00, rho_01, rho_10, rho_11
      INTEGER,DIMENSION(8,3) :: nn
      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      rho000 = rho_val(chg,nn(1,1),nn(1,2),nn(1,3))
      rho001 = rho_val(chg,nn(2,1),nn(2,2),nn(2,3))
      rho010 = rho_val(chg,nn(3,1),nn(3,2),nn(3,3))
      rho100 = rho_val(chg,nn(4,1),nn(4,2),nn(4,3))
      rho011 = rho_val(chg,nn(5,1),nn(5,2),nn(5,3))
      rho101 = rho_val(chg,nn(6,1),nn(6,2),nn(6,3))
      rho110 = rho_val(chg,nn(7,1),nn(7,2),nn(7,3))
      rho111 = rho_val(chg,nn(8,1),nn(8,2),nn(8,3))
      rho00_ = rho000*g3 + rho001*f3
      rho01_ = rho010*g3 + rho011*f3
      rho10_ = rho100*g3 + rho101*f3
      rho11_ = rho110*g3 + rho111*f3
      rho0__ = rho00_*g2 + rho01_*f2
      rho1__ = rho10_*g2 + rho11_*f2
      rho = rho0__*g1 + rho1__*f1
  ! More work for gradients
      rho_0_ = rho00_*g1 + rho10_*f1
      rho_1_ = rho01_*g1 + rho11_*f1
      rho_00 = rho000*g1 + rho100*f1
      rho_01 = rho001*g1 + rho101*f1
      rho_10 = rho010*g1 + rho110*f1
      rho_11 = rho011*g1 + rho111*f1
      rho__0 = rho_00*g2 + rho_10*f2
      rho__1 = rho_01*g2 + rho_11*f2
      rho_grad_lat(1) = rho1__ - rho0__
      rho_grad_lat(2) = rho_1_ - rho_0_
      rho_grad_lat(3) = rho__1 - rho__0
  !   CALL vector_matrix(rho_grad_lat, chg%car2lat, rho_grad)
      nn_grad = MATMUL(chg%car2lat,rho_grad_lat)
    RETURN
    END FUNCTION nn_grad

    ! this funciton takes in 8 values, return a
    ! trilinear interpolated gradient of the values.
    ! the 8 value list value order is 
    ! 000 001 010 100 011 101 110 111
    ! Note 02042019: the above order is what I wrote previously
    ! Note 02042019: I believe the actuall order is the following
    ! 000 100 010 110 001 101 011 111
    !  1   2   3   4   5   6   7   8
    ! r is the indice of the predicted critical point
    ! The interpolation result is checked to be OK by mathematica
    ! USED IN THIS MODULE
    FUNCTION trilinear_interpol_grad(vals,r)
      ! varls come nngrad
      ! could r be the problem?
      ! r needs to be a vector where each component 
      ! number is between 0 and 1
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2),DIMENSION(3) :: trilinear_interpol_grad
      INTEGER :: p1, p2, p3
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2),DIMENSION(8,3) :: vals
      !REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110, rho111
      REAL(q2),DIMENSION(3) :: val00_, val01_, val10_, val11_
      REAL(q2),DIMENSION(3) :: val0__, val1__
 
      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      ! f1 f2 f3 are checked to be correct. 
      ! they should equal to tem for the first step
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      val00_ = vals(1,:)*g3 + vals(5,:)*f3
      val01_ = vals(3,:)*g3 + vals(7,:)*f3
      val10_ = vals(2,:)*g3 + vals(6,:)*f3
      val11_ = vals(4,:)*g3 + vals(8,:)*f3
      val0__ = val00_*g2 + val01_*f2
      val1__ = val10_*g2 + val11_*f2
      trilinear_interpol_grad = val0__*g1 + val1__*f1
    RETURN
    END FUNCTION trilinear_interpol_grad
 
    FUNCTION trilinear_interpol_rho(chg,r)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2) :: trilinear_interpol_rho
      INTEGER :: p1, p2, p3
      INTEGER, DIMENSION(3) :: tempp
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2),DIMENSION(8) :: vals
      !REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110, rho111
      REAL(q2) :: val01_, val10_, val11_
      REAL(q2) :: val0__, val1__
      REAL(q2) :: val00_
 
      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      tempp = (/p1,p2,p3/)
      CALL pbc(tempp,chg%npts)
      vals(1) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1+1,p2,p3/)
      CALL pbc(tempp,chg%npts)
      vals(2) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1,p2+1,p3/)
      CALL pbc(tempp,chg%npts)
      vals(3) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1+1,p2+1,p3/)
      CALL pbc(tempp,chg%npts)
      vals(4) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1,p2,p3+1/)
      CALL pbc(tempp,chg%npts)
      vals(5) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1+1,p2,p3+1/)
      CALL pbc(tempp,chg%npts)
      vals(6) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1,p2+1,p3+1/)
      CALL pbc(tempp,chg%npts)
      vals(7) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1+1,p2+1,p3+1/)
      CALL pbc(tempp,chg%npts)
      vals(8) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      ! f1 f2 f3 are checked to be correct. 
      ! they should equal to tem for the first step
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      val00_ = vals(1)*g3 + vals(5)*f3
      val01_ = vals(3)*g3 + vals(7)*f3
      val10_ = vals(2)*g3 + vals(6)*f3
      val11_ = vals(4)*g3 + vals(8)*f3
      val0__ = val00_*g2 + val01_*f2
      val1__ = val10_*g2 + val11_*f2
      trilinear_interpol_rho = val0__*g1 + val1__*f1
    RETURN
    END FUNCTION

    ! This function outputs hessian in direct coordinates
    ! USED IN THIS MODULE
    FUNCTION trilinear_interpol_hes(vals,r)
      ! varls come nnhes
      ! r needs to be a vector where each component
      ! number is between 0 and 1
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2),DIMENSION(3,3) :: trilinear_interpol_hes
      INTEGER :: p1, p2, p3
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2),DIMENSION(8,3,3) :: vals
      !REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110,
      !rho111
      REAL(q2),DIMENSION(3,3) :: val00_, val01_, val10_, val11_
      REAL(q2),DIMENSION(3,3) :: val0__, val1__

      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      ! f1 f2 f3 are checked to be correct.
      ! they should equal to tem for the first step
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      val00_(1,:) = vals(1,1,:)*g3 + vals(5,1,:)*f3
      val00_(2,:) = vals(1,2,:)*g3 + vals(5,2,:)*f3
      val00_(3,:) = vals(1,3,:)*g3 + vals(5,3,:)*f3
      val01_(1,:) = vals(3,1,:)*g3 + vals(7,1,:)*f3
      val01_(2,:) = vals(3,2,:)*g3 + vals(7,2,:)*f3
      val01_(3,:) = vals(3,3,:)*g3 + vals(7,3,:)*f3
      val10_(1,:) = vals(2,1,:)*g3 + vals(6,1,:)*f3
      val10_(2,:) = vals(2,2,:)*g3 + vals(6,2,:)*f3
      val10_(3,:) = vals(2,3,:)*g3 + vals(6,3,:)*f3
      val11_(1,:) = vals(4,1,:)*g3 + vals(8,1,:)*f3
      val11_(2,:) = vals(4,2,:)*g3 + vals(8,2,:)*f3
      val11_(3,:) = vals(4,3,:)*g3 + vals(8,3,:)*f3
      val0__ = val00_*g2 + val01_*f2
      val1__ = val10_*g2 + val11_*f2
      trilinear_interpol_hes = val0__*g1 + val1__*f1
    END FUNCTION


    
    ! The following subroutine gets gradient using central difference
    ! converts the gradient from lattice to cartesian by the end of it
    ! Note that the gradient is contravariant
    ! USED IN THIS MODULE
    FUNCTION CDGrad(p,chg)
      TYPE(charge_obj) :: chg
      INTEGER, DIMENSION(3) :: p
      INTEGER, DIMENSION(3) :: pzm,pzp,pxm,pxp,pym,pyp
      REAL(q2), DIMENSION(3) :: CDGrad
      pzm = p + (/0,0,-1/)
      pzp = p + (/0,0,1/)
      pxm = p + (/-1,0,0/)
      pxp = p + (/1,0,0/)
      pym = p + (/0,-1,0/)
      pyp = p + (/0,1,0/)
      CALL pbc(pxm,chg%npts)
      CALL pbc(pym,chg%npts)
      CALL pbc(pzm,chg%npts)
      CALL pbc(pxp,chg%npts)
      CALL pbc(pyp,chg%npts)
      CALL pbc(pzp,chg%npts)
      CDGrad(3) = 0.5*(rho_val(chg,pzp(1),pzp(2),pzp(3)) - &
                  rho_val(chg,pzm(1),pzm(2),pzm(3)))
      CDGrad(2) = 0.5*(rho_val(chg,pyp(1),pyp(2),pyp(3)) - &
                  rho_val(chg,pym(1),pym(2),pym(3)))
      CDGrad(1) = 0.5*(rho_val(chg,pxp(1),pxp(2),pxp(3)) - &
                  rho_val(chg,pxm(1),pxm(2),pxm(3)))
      CDGrad = MATMUL(CDGrad,chg%car2lat)
      RETURN
      ! now the gradient should be in cartesian
    END FUNCTION CDGrad
    
    ! the following subroutine gets hes and force in lattice units and converts
    ! to cartesian by the end
    ! USED IN THIS MODULE
    FUNCTION CDHessian(p,chg)
      REAL(q2),DIMENSION(3,3) :: CDHessian
      TYPE(charge_obj) :: chg
      INTEGER,DIMENSION(3) :: ptxy1, ptxy2, ptxy3, ptxy4, ptxz1, ptxz2, ptxz3 
      INTEGER,DIMENSION(3) :: ptxz4, ptyz1, ptyz2, ptyz3, ptyz4
      INTEGER,DIMENSION(3) :: p, ptx1, ptx2, pty1, pty2, ptz1, ptz2
      CALL pbc(p,chg%npts)
      ptx1 = p + (/1,0,0/)
      ptx2 = p + (/-1,0,0/)
      pty1 = p + (/0,1,0/)
      pty2 = p + (/0,-1,0/)
      ptz1 = p + (/0,0,1/)
      ptz2 = p + (/0,0,-1/)
      CALL pbc(ptx1,chg%npts)
      CALL pbc(ptx2,chg%npts)
      CALL pbc(pty1,chg%npts)
      CALL pbc(pty2,chg%npts)
      CALL pbc(ptz1,chg%npts)
      CALL pbc(ptz2,chg%npts)
      CDHessian(1,1) = &
         (rho_val(chg,ptx1(1),ptx1(2),ptx1(3)) - &
         rho_val(chg,p(1),p(2),p(3))) - &
         (rho_val(chg,p(1),p(2),p(3)) - &
         rho_val(chg,ptx2(1),ptx2(2),ptx2(3))) 
      CDHessian(2,2) = &
        (rho_val(chg,pty1(1),pty1(2),pty1(3)) - &
         rho_val(chg,p(1),p(2),p(3)))  - &
        (rho_val(chg,p(1),p(2),p(3)) - &
         rho_val(chg,pty2(1),pty2(2),pty2(3))) 
      CDHessian(3,3) = &
        (rho_val(chg,ptz1(1),ptz1(2),ptz1(3)) - &
         rho_val(chg,p(1),p(2),p(3)))  - &
        (rho_val(chg,p(1),p(2),p(3)) - &
         rho_val(chg,ptz2(1),ptz2(2),ptz2(3))) 
      ptxy1 = p + (/-1,-1,0/)
      ptxy2 = p + (/-1,+1,0/)
      ptxy3 = p + (/+1,+1,0/)
      ptxy4 = p + (/+1,-1,0/)
      ptxz1 = p + (/-1,0,-1/)
      ptxz2 = p + (/-1,0,+1/)
      ptxz3 = p + (/+1,0,+1/)
      ptxz4 = p + (/+1,0,-1/)
      ptyz1 = p + (/0,-1,-1/)
      ptyz2 = p + (/0,-1,+1/)
      ptyz3 = p + (/0,+1,+1/)
      ptyz4 = p + (/0,+1,-1/)
      CALL pbc(ptxy1,chg%npts)
      CALL pbc(ptxy2,chg%npts)
      CALL pbc(ptxy3,chg%npts)
      CALL pbc(ptxy4,chg%npts)
      CALL pbc(ptxz1,chg%npts)
      CALL pbc(ptxz2,chg%npts)
      CALL pbc(ptxz3,chg%npts)
      CALL pbc(ptxz4,chg%npts)
      CALL pbc(ptyz1,chg%npts)
      CALL pbc(ptyz2,chg%npts)
      CALL pbc(ptyz3,chg%npts)
      CALL pbc(ptyz4,chg%npts)
      CDHessian(1,2) = &
        ( &
        ! this is the backward dv
!        - 0.5_q2 / vMag * &
        - 0.25_q2 * & 
        ((rho_val(chg,ptxy2(1),ptxy2(2),ptxy2(3)) + &
            rho_val(chg,pty1(1),pty1(2),pty1(3)))  - &
          (rho_val(chg,ptxy1(1),ptxy1(2),ptxy1(3)) + &
            rho_val(chg,pty2(1),pty2(2),pty2(3)))  )  & 
        ! this is the forward dv
        + 0.25_q2 * & 
        ((rho_val(chg,ptxy3(1),ptxy3(2),ptxy3(3)) + &
            rho_val(chg,pty1(1),pty1(2),pty1(3)))  - &
          (rho_val(chg,ptxy4(1),ptxy4(2),ptxy4(3)) + &
            rho_val(chg,pty2(1),pty2(2),pty2(3)))  )  &
        )
      CDHessian(2,1) = CDHessian(1,2)
      CDHessian(1,3) = &
        ( &
        ! this is the bacward dw
        - 0.25_q2 * &
        ((rho_val(chg,ptxz2(1),ptxz2(2),ptxz2(3)) + &
            rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))   - &
          (rho_val(chg,ptxz1(1),ptxz1(2),ptxz1(3)) + &
            rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))   )   &
        ! this is the forward dw
        + 0.25_q2 * & 
        ((rho_val(chg,ptxz3(1),ptxz3(2),ptxz3(3)) + &
            rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))   - &
          (rho_val(chg,ptxz4(1),ptxz4(2),ptxz4(3)) + &
            rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))   )  &
        ) 
      CDHessian(3,1) = CDHessian(1,3)
      CDHessian(2,3) = &
        ( &
        ! this is the bacward dw
        - 0.25_q2 * &
        ((rho_val(chg,ptyz2(1),ptyz2(2),ptyz2(3)) + &
            rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))  - &
          (rho_val(chg,ptyz1(1),ptyz1(2),ptyz1(3)) + &
            rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))  )  & 
        ! this is the forward dw
        + 0.25_q2 * & 
        ((rho_val(chg,ptyz3(1),ptyz3(2),ptyz3(3)) + &
            rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))  - &
          (rho_val(chg,ptyz4(1),ptyz4(2),ptyz4(3)) + &
            rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))  )  &
        )
      CDHessian(3,2) = CDHessian(2,3)
      ! Convert the hessian, which is now in lattice coordinates, to cartesian
      !CDHessian = MATMUL(MATMUL(chg%car2lat,CDHessian),TRANSPOSE(chg%car2lat))
      CDHessian = MATMUL(TRANSPOSE(chg%car2lat),MATMUL(CDHessian,chg%car2lat))
      !CDHessian = MATMUL(CDHessian,chg%car2lat)
      
      RETURN
    END FUNCTION CDHessian
    

    ! This function takes in current position and grid point position in
    ! lattice, finds the distance this point is to the nearest cell, take half
    ! the distance as step size. 
    FUNCTION findstepsize(r)
      REAL(q2) :: findstepsize    
      REAL(q2),DIMENSION(3) :: r
      REAL(q2) :: f1,f2,f3
      f1 = MIN(ABS(r(1)),ABS(0.5-ABS(r(1)))) 
      f2 = MIN(ABS(r(2)),ABS(0.5-ABS(r(2))))
      f3 = MIN(ABS(r(3)),ABS(0.5-ABS(r(3))))
      findstepsize = MIN(MIN(f1,f2),MIN(f1,f3))/2
    RETURN
    END FUNCTION findstepsize

    ! this funciton finds hessian of a interpolated point using interpolated
    ! nearest neighbor gradients. Also gradients taken in here should be in
    ! lattice.
    FUNCTION inthessian(grad,stepsize)
      REAL(q2),DIMENSION(6,3) :: grad
      REAL(q2),DIMENSION(3,3) :: inthessian
      REAL(q2) :: stepsize
      ! again, intnngrad is following this order:
      ! +x -x +y -y +z -z
      inthessian(1,1) = (grad(1,1)-grad(2,1))*0.5_q2/stepsize
      inthessian(2,2) = (grad(3,2)-grad(4,2))*0.5_q2/stepsize
      inthessian(3,3) = (grad(5,3)-grad(6,3))*0.5_q2/stepsize
      inthessian(1,2) = (grad(3,1)-grad(4,1))*0.5_q2/stepsize
      inthessian(2,1) = inthessian(1,2)
      inthessian(1,3) = (grad(5,1)-grad(6,1))*0.5_q2/stepsize
      inthessian(3,1) = inthessian(1,3)
      inthessian(2,3) = (grad(5,2)-grad(6,2))*0.5_q2/stepsize
      inthessian(3,2) = inthessian(2,3)
    ! assuming that this function is fine
    RETURN
    END FUNCTION inthessian

   


    ! USED IN THIS MODULE
    FUNCTION makevi()
      INTEGER, DIMENSION(3,26) :: makevi
      makevi(:,1)=(/-1,-1,-1/)
      makevi(:,2)=(/-1,-1,0/)
      makevi(:,3)=(/-1,-1,1/)
      makevi(:,4)=(/-1,0,-1/)
      makevi(:,5)=(/-1,0,0/)
      makevi(:,6)=(/-1,0,1/)
      makevi(:,7)=(/-1,1,-1/)
      makevi(:,8)=(/-1,1,0/)
      makevi(:,9)=(/-1,1,1/)
      makevi(:,10)=(/0,-1,-1/)
      makevi(:,11)=(/0,-1,0/)
      makevi(:,12)=(/0,-1,1/)
      makevi(:,13)=(/0,0,-1/)
      makevi(:,14)=(/1,1,1/)
      makevi(:,15)=(/1,1,0/)
      makevi(:,16)=(/1,1,-1/)
      makevi(:,17)=(/1,0,1/)
      makevi(:,18)=(/1,0,0/)
      makevi(:,19)=(/1,0,-1/)
      makevi(:,20)=(/1,-1,1/)
      makevi(:,21)=(/1,-1,0/)
      makevi(:,22)=(/1,-1,-1/)
      makevi(:,23)=(/0,1,1/)
      makevi(:,24)=(/0,1,0/)
      makevi(:,25)=(/0,1,-1/)
      makevi(:,26)=(/0,0,1/)
      RETURN
    END FUNCTION

    ! USED IN THIS MODULE
    FUNCTION makeggrid(chg,ions)
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      REAL(q2),DIMENSION(3,3) :: makeggrid
      makeggrid(1,1) = DOT_PRODUCT(ions%lattice(1,:),ions%lattice(1,:))/ &
        (chg%npts(1)*chg%npts(1))
      makeggrid(1,2) = DOT_PRODUCT(ions%lattice(1,:),ions%lattice(2,:))/ &
        (chg%npts(1)*chg%npts(2))
      makeggrid(1,3) = DOT_PRODUCT(ions%lattice(1,:),ions%lattice(3,:))/ &
        (chg%npts(1)*chg%npts(3))
      makeggrid(2,1) = DOT_PRODUCT(ions%lattice(2,:),ions%lattice(1,:))/ &
        (chg%npts(2)*chg%npts(1))
      makeggrid(2,2) = DOT_PRODUCT(ions%lattice(2,:),ions%lattice(2,:))/ &
        (chg%npts(2)*chg%npts(2))
      makeggrid(2,3) = DOT_PRODUCT(ions%lattice(2,:),ions%lattice(3,:))/ &
        (chg%npts(2)*chg%npts(3))
      makeggrid(3,1) = DOT_PRODUCT(ions%lattice(3,:),ions%lattice(1,:))/ &
        (chg%npts(3)*chg%npts(1))
      makeggrid(3,2) = DOT_PRODUCT(ions%lattice(3,:),ions%lattice(2,:))/ &
        (chg%npts(3)*chg%npts(2))
      makeggrid(3,3) = DOT_PRODUCT(ions%lattice(3,:),ions%lattice(3,:))/ &
        (chg%npts(3)*chg%npts(3))
      RETURN
    END FUNCTION

    ! USED IN THIS MODULE
    FUNCTION ascension(ind,chg, &
                       ggrid,opts,ions)
      ! this function finds nucleus critical points. 
      INTEGER, DIMENSION(3) :: ind
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      TYPE(ions_obj) :: ions
      REAL(q2),DIMENSION(3) :: ascension, distance
      INTEGER, DIMENSION(:,:),ALLOCATABLE :: nnind
      REAL(q2), DIMENSION(8,3) :: nngrad
      INTEGER :: j, stepcount
      REAL(q2), DIMENSION(3) :: tempr, rn, rnm1 ! rn minus 1
      REAL(q2), DIMENSION(3) :: grad, stepsize, gradnm1
      REAL(q2), DIMENSION(3,3) :: ggrid
      
      stepcount = 0
      ALLOCATE(nnInd(8,3))
!      ALLOCATE(nnInd((nnlayers*2+1)**3,3))
      ! initialize the process
      CALL pbc(ind,chg%npts)
      rn(1) = REAL(ind(1),q2)
      rn(2) = REAL(ind(2),q2)
      rn(3) = REAL(ind(3),q2)
      stepsize(1) = 0.5
      stepsize(2) = 0.5
      stepsize(3) = 0.5
      grad = cdgrad(ind,chg)
      ! this gradient is in cartesian. convert it to lattice
      grad = MATMUL(chg%lat2car,grad)
      DO WHILE (stepsize(1) >= 0.01 .AND. &
                stepsize(2) >= 0.01 .AND. &
                stepsize(3) >= 0.01 )
        ! the gradient is in cartesian. 
        gradnm1 = grad
        rnm1 = rn
        ! determine where to go, a unit vector
        tempr = grad / SQRT(SUM(grad*grad))
        tempr(1) = stepsize(1) * tempr(1)
        tempr(2) = stepsize(2) * tempr(2)
        tempr(3) = stepsize(3) * tempr(3)
        rn = rn + tempr
        CALL pbc_r_lat(rn,chg%npts)
        nnind = SimpleNN(rn,chg)
        DO j = 1,8
          CALL pbc(nnind(j,:),chg%npts)
          nngrad(j,:) = cdgrad(nnind(j,:),chg)
        END DO
        distance = rn - nnind(1,:)
         !Row find the nearest neighbors at this new locaiton
         !First update critical point location
         !The next big step is to interpolate the force at predicted critical
         !point.
        grad = trilinear_interpol_grad(nngrad,distance) ! val r interpol
        !nnind = FindNN(rn,nnLayers,chg,ions)
        !grad = R2GradInterpol(nnind,rn,chg,nnLayers)
        ! this grad is in cartesian. convert it to lattice
        grad = MATMUL(grad,chg%lat2car)
        IF (ABS(grad(1)) < 0.001 &
            .AND. ABS(grad(2)) <= 0.001 &
            .AND. ABS(grad(3)) <= 0.001)  THEN
          ! we are at a critical point!
!          PRINT *, 'ascention gradient sufficiently small'
!          PRINT *, 'ascension rn is'
!          PRINT *, rn
          EXIT
        END IF
        ! if grad points backwards, reduce stepsize
        IF (DOT_PRODUCT(grad,gradnm1) <= 0) THEN
          stepsize = 0.5*stepsize
        END IF
        !IF (grad(1) * gradnm1(1) < 0) THEN
        !  stepsize(1) = stepsize(1)/2
        !END IF
        !IF (grad(2) * gradnm1(2) < 0) THEN
        !  stepsize(2) = stepsize(2)/2
        !END IF
        !IF (grad(3) * gradnm1(3) < 0) THEN
        !  stepsize(3) = stepsize(3)/2
        !END IF
        stepcount = stepcount + 1
!        PRINT *, rn
      END DO 
      ascension = rn
!      PRINT *, 'ascension step count is ', stepcount
      RETURN 
    END FUNCTION ascension
   

    ! USED IN THIS MODULE
    SUBROUTINE RingAscension(cp,chg,ions, UniqueCPs)
      !Starting from a ring critical point, locates the positions of associated nuclei critical points
      REAL(q2), DIMENSION(3) :: ind
      REAL(q2), DIMENSION(3) :: vec1, vec2,vecsum
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      TYPE(cpc) :: cp
      INTEGER, DIMENSION(:), ALLOCATABLE :: UniqueCPs
      INTEGER, DIMENSION(:), ALLOCATABLE :: CPsFound1
      INTEGER, DIMENSION(:), ALLOCATABLE :: CPsFound2
      INTEGER :: i, divs, divs2
      LOGICAL :: foundring
      
      ind = cp%truer
      !checks to make sure is a ring CP
      divs = 10
      divs2 = 30
      ALLOCATE(CPsFound1(divs))
      ALLOCATE(CPsFound2(divs2))
      CPsFound1 = CPsFound1*0
      CPsFound2 = CPsFound2*0
      IF (cp%eigvals(1) < 0 .AND. cp%eigvals(2) > 0 .AND. cp%eigvals(3) > 0 ) THEN
        vec1(1) = cp%eigvecs(1,2)
        vec1(2) = cp%eigvecs(2,2)
        vec1(3) = cp%eigvecs(3,2)
        
        vec2(1) = cp%eigvecs(1,3)
        vec2(2) = cp%eigvecs(2,3)
        vec2(3) = cp%eigvecs(3,3)
        foundring = .TRUE.
      ELSE IF (cp%eigvals(1) > 0 .AND. cp%eigvals(2) < 0 .AND. cp%eigvals(3) > 0 ) THEN
        vec1(1) = cp%eigvecs(1,1)
        vec1(2) = cp%eigvecs(2,1)
        vec1(3) = cp%eigvecs(3,1)

        vec2(1) = cp%eigvecs(1,3)
        vec2(2) = cp%eigvecs(2,3)
        vec2(3) = cp%eigvecs(3,3)
        foundring = .TRUE.
      ELSE IF (cp%eigvals(1) > 0 .AND. cp%eigvals(2) > 0 .AND. cp%eigvals(3) < 0 ) THEN
        vec1(1) = cp%eigvecs(1,1)
        vec1(2) = cp%eigvecs(2,1)
        vec1(3) = cp%eigvecs(3,1)
        
        vec2(1) = cp%eigvecs(1,2)
        vec2(2) = cp%eigvecs(2,2)
        vec2(3) = cp%eigvecs(3,2)
        foundring = .TRUE.
      ELSE
       ! PRINT *, "Not a ring CP"
        vec1 = (/0.0, 0.0, 0.0/)
        vec2 = (/0.0, 0.0, 0.0/)
        foundring = .FALSE.
      END IF
      vec1 = vec1/Mag(vec1)
      vec2 = vec2 - DOT_PRODUCT(vec1,vec2)/DOT_PRODUCT(vec1,vec1)*vec1
      vec2 = vec2/Mag(vec2)
      

      IF (foundring) THEN
        DO i=1,divs
          vecsum = vec1*COS(i*2*3.1415/divs) + vec2*SIN(i*2*3.1415/divs)
          CPsFound1(i) = ascension_new(ind+vecsum, chg, ions)
         ! PRINT *,"Nuclei found 1st run: ", CPsFound1(i)
        END DO

        vec1 = 2*vec1
        vec2 = 2*vec2
        PRINT *, " "
        DO i=1,divs2
          vecsum = vec1*COS(i*2*3.1415/divs2) + vec2*SIN(i*2*3.1415/divs2)
          CPsFound2(i) = ascension_new(ind+vecsum, chg, ions)
         ! PRINT *, "Nuclei found 2nd run: ", CPsFound2(i)

        END DO
       ! PRINT *,"CPsFound2: ",  CPsFound2
      END IF
      
      CALL unique(CPsFound1,UniqueCPs)
      IF (foundring) PRINT *, divs, " run: ", UniqueCPs
      CALL unique(CPsFound2,UniqueCPs)
      IF (foundring) PRINT *, divs2, " run: ", UniqueCPs 
      IF (foundring) PRINT *, " "
    END SUBROUTINE RingAscension
    
    !helper subroutine for RingAscension which takes an array of integers and returns an array of the unique elements of that array
    ! USED IN THIS MODULE
    SUBROUTINE unique(vec,vec_unique)
      IMPLICIT NONE
      INTEGER,DIMENSION(:),INTENT(in) :: vec
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: vec_unique
      INTEGER :: i, num
      LOGICAL, DIMENSION(size(vec)) :: mask

      mask = .FALSE.
      DO i=1,size(vec)
        num = count(vec(i)==vec)

        IF (num==1) THEN
          mask(i) = .TRUE.
        ELSE
          IF (.NOT. ANY(vec(i)==vec .AND. mask)) mask(i) = .TRUE.
        END IF
      END DO

      ALLOCATE(vec_unique(count(mask)) )
      vec_unique = pack(vec,mask)
    END SUBROUTINE unique

    ! removes the zero/NaN entries in a list of real-valued vectors 
    ! USED IN THIS MODULE
    SUBROUTINE unique_realcoords(vec,vec_unique)
      IMPLICIT NONE
      REAL(q2), DIMENSION(:,:), INTENT(in) :: vec
      REAL(q2), DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: vec_unique
      INTEGER :: i
      LOGICAL, DIMENSION(size(vec,1),size(vec,2)) :: mask
      mask = .TRUE.
      DO i=1,size(vec,1)
        IF (vec(i,1) .NE. vec(i,1) .OR. vec(i,2) .NE. vec(i,2) .OR. vec(i,3) .NE. vec(i,3) ) THEN
          mask(i,:) = .FALSE.
        ELSE IF (Mag(vec(i,:)) <1) THEN
          mask(i,:) = .FALSE.
        END IF
      END DO
      ALLOCATE(vec_unique(count(mask)/3, 3) )
      vec_unique =reshape( pack(vec,mask), (/count(mask)/3, 3 /) )
    END SUBROUTINE unique_realcoords

    ! USED IN THIS MODULE
    FUNCTION DoubleAscension(bdr,cp,chg,ions,opts,direction,cpl,ucptnum)
      !cp is a selected critical point candidate (cpc) object
      !performs ascension starting from an input bond CP to find the array indices of the two associated nuclei CPs.
      !will gradient ascend from two points shifted in opposite directions around it following the bond eigenvector. 
      !Returns the index numbers of the two nuclei found as an array of size 2
      REAL(q2), DIMENSION(3) :: ind
      REAL(q2), DIMENSION(3) ::vector
      TYPE(cpc), ALLOCATABLE, DIMENSION(:) :: cpl
      TYPE(bader_obj) :: bdr
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      TYPE(options_obj) :: opts
      REAL(q2), DIMENSION(3) :: ind_plus, ind_minus
      REAL(q2), DIMENSION(3) :: traj_end
      TYPE(cpc) :: cp
      INTEGER :: vecnum, nuc1num, nuc2num, direction, ucptnum
      INTEGER, DIMENSION(3) :: iinds
      INTEGER, DIMENSION(2) :: DoubleAscension
      !assigns ind to the starting critical point coordinate 
      ind = cp%truer
      IF (direction == 1) THEN
        !checks to make sure only one eigenvalue is positive (indicating bond CP), and finds which one it is.
        IF (cp%eigvals(1) < 0 .AND. cp%eigvals(2) < 0 .AND. cp%eigvals(3) > 0 ) THEN
          vecnum = 3
        ELSE IF (cp%eigvals(1) < 0 .AND. cp%eigvals(2) > 0 .AND. cp%eigvals(3) <0 ) THEN
          vecnum = 2
        ELSE IF (cp%eigvals(1) > 0 .AND. cp%eigvals(2) < 0 .AND. cp%eigvals(3) < 0 ) THEN
          vecnum = 1
        ELSE
          vecnum = 0
        END IF
      ELSEIF (direction == -1) THEN
        !checks to make sure only one eigenvalue is negative (indicating ring CP), and finds which one it is.
        IF (cp%eigvals(1) > 0 .AND. cp%eigvals(2) > 0 .AND. cp%eigvals(3) < 0 ) THEN
          vecnum = 3
        ELSE IF (cp%eigvals(1) > 0 .AND. cp%eigvals(2) < 0 .AND. cp%eigvals(3) > 0 ) THEN
          vecnum = 2
        ELSE IF (cp%eigvals(1) < 0 .AND. cp%eigvals(2) > 0 .AND. cp%eigvals(3) > 0 ) THEN
          vecnum = 1
        ELSE
          vecnum = 0
        END IF
      END IF
      nuc1num = 0
      nuc2num = 0
      IF (vecnum == 0) THEN
       ! PRINT *, "Not a valid bond critical point! Exactly one positive eigenvalue is required to proceed"
      ELSE
        vector = cp%eigvecs(:,vecnum)
        vector = vector/Mag(vector)
        ind_plus = ind + vector
        ind_minus = ind - vector
        CALL pbc_r_lat(ind_plus, chg%npts)
        CALL pbc_r_lat(ind_minus, chg%npts)
        IF (direction == 1) THEN
          nuc1num = ascension_new(ind_plus,chg,ions)
          nuc2num = ascension_new(ind_minus,chg,ions)
        ELSEIF (direction == -1) THEN
          iinds = NINT(ind + 2 * vector)
          CALL pbc(iinds,chg%npts)
          traj_end = DensityDescend(chg,bdr,opts,iinds)
          CALL pbc_r_lat(traj_end, chg%npts)
          nuc1num = CageSearch(traj_end,cpl,ucptnum,bdr,chg)
          iinds = NINT(ind - 2 * vector)
          CALL pbc(iinds,chg%npts)
          traj_end = DensityDescend(chg,bdr,opts,iinds)
          CALL pbc_r_lat(traj_end, chg%npts)
          nuc2num = CageSearch(traj_end,cpl,ucptnum,bdr,chg)
        END IF
      END IF
      DoubleAscension(1) = nuc1num
      DoubleAscension(2) = nuc2num
      RETURN
    END FUNCTION DoubleAscension

    ! search in the cp list and find a close by cage point 
    FUNCTION CageSearch(r,cpl,ucptnum,bdr,chg)
      TYPE(charge_obj) :: chg
      TYPE(bader_obj) :: bdr
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      REAL(q2),DIMENSION(3) :: r
      REAL(q2) :: min_distance, distance
      INTEGER, DIMENSION(3) :: p
      INTEGER :: CageSearch, ucptnum, nearest_cage, n
      p = NINT(r)
      CALL pbc(p,chg%npts)
      IF (bdr%volnum(p(1),p(2),p(3)) == bdr%bnum + 1) THEN
        CageSearch = -1 !connected to vacuum
      ELSE
        min_distance = 9999
        DO n=1, ucptnum
          IF (cpl(n)%negCount == 0) THEN
            distance = SQRT(SUM((cpl(n)%trueR(1:3) - r(1:3))**2) )
            IF (distance < min_distance) THEN
              min_distance = distance
              CageSearch = n
            END IF
          END IF        
        END DO
        IF (min_distance > 1) THEN
          CageSearch = -2 !no cage found
        END IF
      END IF
      RETURN
    END FUNCTION CageSearch


    ! USED IN THIS MODULE
    FUNCTION Ascension_New(ind, chg, ions)
    !used by DoubleAscension  
    !performs gradient ascent on a real coordinate, returns list index of the first nuclei found within 0.5 units of position .  
      REAL(q2), DIMENSION(3) ::startpos, ind
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      
      REAL(q2), DIMENSION(3) :: grad,oldgrad, distance
      REAL(q2) :: maxstepsize
      INTEGER,DIMENSION(8,3) :: nnInd
      REAL(q2), DIMENSION(8,3) :: nngrad
      INTEGER :: i, j, loopcount, stepMax
      INTEGER :: ascension_new
      LOGICAL :: foundAtom
      CALL pbc_r_lat(ind, chg%npts)
      startpos = ind
      nnInd = SimpleNN(ind, chg)
      DO j=1,8
        CALL pbc(nnInd(j,:),chg%npts)
        nngrad(j,:) = CDGrad(nnInd(j,:),chg)
      END DO
      distance = ind - nnInd(1,:)
      grad = trilinear_interpol_grad(nnGrad,distance)
      grad = MATMUL(grad,chg%lat2car)
      foundAtom = .FALSE.
      ascension_new = 0
      maxstepsize = 0.5
      loopcount = 0
      stepMax = 1000
      DO WHILE(loopcount < stepMax)
        loopcount = loopcount + 1
       ! PRINT *, "ind", ind
        DO i=1,ions%nions
          ! PRINT *, ions%r_lat(i,:)
          IF (Mag(MATMUL(ions%r_lat(i,:),chg%lat2car)-MATMUL(ind,chg%lat2car)) < 0.15 ) THEN
           ! PRINT *, "Found an atom within 0.5" 
            ascension_new = i
            foundAtom = .TRUE.
            !PRINT *, "Cartesian Distance: ", Mag(MATMUL(ind,chg%lat2car)-MATMUL(startpos,chg%lat2car))
            EXIT
          END IF
        END DO
        IF (foundAtom) THEN
          EXIT
        END IF
        grad = MIN(Mag(grad),maxstepsize) * grad/Mag(grad)
        ind = ind + grad
        oldgrad = grad
        CALL pbc_r_lat(ind, chg%npts)
        nnInd = SimpleNN(ind,chg)
        DO j=1,8
          CALL pbc(nnInd(j,:),chg%npts)
          nngrad(j,:) = CDGrad(nnInd(j,:),chg)
        END DO
        distance = ind - nnInd(1,:)
        grad = trilinear_interpol_grad(nnGrad,distance)
        grad = MATMUL(grad,chg%lat2car)
        IF (SUM(grad*oldgrad)<0) THEN
          maxstepsize = 1*maxstepsize
        END IF
      END DO
      IF (.NOT. foundAtom) THEN
       ! PRINT *, "No atom was found"
      END IF
      RETURN
    END FUNCTION ascension_new



    ! USED IN THIS MODULE
    FUNCTION Descension(ind,chg, &
                        opts,ions)
      ! this function finds nucleus critical points. 
      INTEGER, DIMENSION(3) :: ind
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      TYPE(ions_obj) :: ions
      REAL(q2),DIMENSION(3) :: descension, distance
      INTEGER, DIMENSION(:,:),ALLOCATABLE :: nnind
      REAL(q2), DIMENSION(8,3) :: nngrad
      INTEGER :: j, stepcount
      REAL(q2), DIMENSION(3) :: tempr, rn, rnm1 ! rn minus 1
      REAL(q2), DIMENSION(3) :: grad, stepsize, gradnm1
      stepcount = 0
      ALLOCATE(nnInd(8,3))
      ! initialize the process
      CALL pbc(ind,chg%npts)
      rn(1) = REAL(ind(1),q2)
      rn(2) = REAL(ind(2),q2)
      rn(3) = REAL(ind(3),q2)
      stepsize(1) = 0.5
      stepsize(2) = 0.5
      stepsize(3) = 0.5
      grad = cdgrad(ind,chg)
      ! this gradient is in cartesian. convert it to lattice
      grad = MATMUL(chg%lat2car,grad)
      DO WHILE (stepsize(1) >= 0.1 .AND. &
                stepsize(2) >= 0.1 .AND. &
                stepsize(3) >= 0.1 )
        ! the gradient is in cartesian. 
        IF (SUM(grad*grad) <= (0.1*opts%par_gradfloor)**2) THEN
          PRINT *, 'descension grad small enough'
          EXIT
        END IF
        gradnm1 = grad
        rnm1 = rn
        ! determine where to go, a unit vector
        tempr = grad / SQRT(SUM(grad*grad))
        
        tempr(1) = stepsize(1) * tempr(1)
        tempr(2) = stepsize(2) * tempr(2)
        tempr(3) = stepsize(3) * tempr(3)
        rn = rn - tempr
        CALL pbc_r_lat(rn,chg%npts)
        nnind = SimpleNN(rn,chg)
        DO j = 1,8
          CALL pbc(nnind(j,:),chg%npts)
          nngrad(j,:) = cdgrad(nnind(j,:),chg)
        END DO
        distance = rn - nnind(1,:)
        grad = trilinear_interpol_grad(nngrad,distance) ! val r interpol
        ! this grad is in cartesian. convert it to lattice
        grad = MATMUL(chg%lat2car,grad)
        IF (SUM(grad*grad)<=(0.1*opts%par_gradfloor)**2) THEN
          ! we are at a critical point!
          EXIT
        END IF
        ! if grad points backwards, reduce stepsize
        IF (grad(1) * gradnm1(1) < 0) THEN
          stepsize(1) = stepsize(1)/2
        END IF
        IF (grad(2) * gradnm1(2) < 0) THEN
          stepsize(2) = stepsize(2)/2
        END IF
        IF (grad(3) * gradnm1(3) < 0) THEN
          stepsize(3) = stepsize(3)/2
        END IF
        stepcount = stepcount + 1
      END DO 
      descension = rn
      RETURN 
    END FUNCTION 

  
    ! this function takes in the current point, see if a near by point has
    ! already been marked as a critical point. if yes, skip it. The criteria is
    ! set to be based on halfo the search radius, which is 1 + knob_tem
    ! THIS FUNCTION AS OF 20191221 WONT WORK WELL AROUND PBC
    ! As a loose first round checking, it is OK that a few points near the PBC
    ! are permitted into the candidacy. 
    ! Just like the Democratic 2020 primary, not everyone has to be super
    ! qualified to enter.
    ! USED IN THIS MODULE
    !FUNCTION ProxyToCPCandidate(p,opts,cpl,cptnum,chg,nnLayers)
    FUNCTION ProxyToCPCandidate(p,opts,cpl,cptnum,chg)
      LOGICAL :: ProxyToCPCandidate
      INTEGER :: i
      INTEGER, DIMENSION(3) :: p
      TYPE(options_obj) :: opts
      TYPE(cpc), ALLOCATABLE, DIMENSION(:) :: cpl
      TYPE(charge_obj) :: chg
      INTEGER :: cptnum
      !INTEGER :: nnLayers
      ProxyToCPCandidate = .FALSE.
      ! These are lattice based and do not take PBC into consideration
      DO i = 1, cptnum
        IF (ABS(cpl(i)%ind(1) - p(1)) <= opts%cp_search_radius .AND. &
            ABS(cpl(i)%ind(2) - p(2)) <= opts%cp_search_radius .AND. &
            ABS(cpl(i)%ind(3) - p(3)) <= opts%cp_search_radius ) THEN
          ProxyToCPCandidate = .TRUE.
        END IF
      END DO  
      RETURN
    END FUNCTION ProxyToCPCandidate

    ! this function should always run. it first checks if things are stuck
    ! this function is for when newton raphson hops between two points
    FUNCTION unstuck( nexttem, previoustem, temscale, temnormcap)
      REAL(q2), DIMENSION(3) :: nexttem, previoustem, temscale, unstuck
      REAL(q2) :: temnormcap
      LOGICAL :: stuck
      LOGICAL :: capprev,capnext ! see if tem was or will be capped
      stuck = .FALSE. ! we check for stuckness
      capprev = .FALSE.
      capnext = .FALSE.
      !temnormcap is all rounded cap designed to prevent overstepping.
      !temnormcap should not change at all
      !temscaleois direction sensitive. It gradually narrows in if tem changes
      !sign a lot.
      IF (Mag(previoustem) == temnormcap) THEN
        capprev = .TRUE.
      END IF
      IF (Mag(nexttem) == temnormcap ) THEN
        capnext = .TRUE.
      END IF
      ! if nexttem and prev tem are complete opposites, it's stuck for sure
      IF (nexttem(1) == -previoustem(1) .AND. &
          nexttem(2) == -previoustem(2) .AND. &
          nexttem(3) == -previoustem(3) ) THEN
        stuck = .TRUE.
      END IF
      IF ( stuck ) THEN
        nexttem = rngkick(nexttem)
      END IF
      unstuck = nexttem
      ! stuck is defined as nexttem and previous tem being equal but complete
      ! opposite

      RETURN
    END FUNCTION

    ! this function takes in a vector, give it a random scale within 10%
    FUNCTION rngkick(nexttem)
      REAL(q2), DIMENSION(3) :: nexttem, rngkick
      REAL(q2) :: ran
      INTEGER :: i
      DO i = 1,3
        CALL RANDOM_NUMBER(ran)
        nexttem(1) = ((0.5-ran)*0.2 + 1)*nexttem(1)
      END DO
      rngkick = nexttem
      RETURN
    END FUNCTION

    ! give the Magnitude of a 3d vector
    ! USED IN THIS MODULE
    FUNCTION Mag(vec3d)
      REAL(q2), DIMENSION(3) :: vec3d
      REAL(q2) :: Mag
      Mag = SQRT(SUM(vec3d*vec3d))
      RETURN
    END FUNCTION

    ! this function takes in the movement factor, modifies it as necessary
    ! USED IN THIS MODULE
    FUNCTION TemMods(nexttem,temscale,temnormcap)
      REAL(q2), DIMENSION(3) :: TemMods, nexttem, temscale
      REAL(q2) :: temnormcap
      INTEGER :: i
      TemMods = nexttem
      DO i = 1 , 3
        TemMods(i) = TemMods(i) * temscale(i)
      END DO
      IF (Mag(TemMods) > temNormCap) THEN
        TemMods = TemMods/Mag(TemMods) * temScale
      END IF
      RETURN
    END FUNCTION
   
    ! this function inspects if it is beneficial to reduce temscale
    ! right now it does nothing because I'm not sure if limiting it helps in any
    ! way at all
    ! USED IN THIS MODULE
    FUNCTION scaleinspector( nexttem, previoustem, temScale)
      REAL(q2), DIMENSION(3) :: nexttem, previoustem, temScale, scaleinspector
      !DO i = 1, 3
      !  IF (nexttem(i)*previoustem(i) <= 0) THEN
      !    temScale(i) = temScale(i) * 0.5
      !  END IF
      !END DO
      !IF (nexttem(1)*previoustem(1)<=0 .AND. &
      !    nexttem(1)*previoustem(1)<=0 .AND. &
      !    nexttem(1)*previoustem(1)<=0 ) THEN
      !    temScale = temScale * 0.5
      !END IF
      scaleinspector = temScale
      RETURN
    END FUNCTION

    ! if it is detected that the hessian 
    FUNCTION unZeroHessian(truer,chg,ggrid)
      TYPE(charge_obj) :: chg
      REAL(q2), DIMENSION(3,3) :: unZeroHessian
      INTEGER, DIMENSION(3) :: truer
      REAL(q2), DIMENSION(3,3) :: ggrid
      REAL(q2) :: ran
      INTEGER :: i
      
      ! first step is to move a little bit
      DO i = 1 , 3
        CALL RANDOM_NUMBER(ran)
        truer(i) = truer(i) + (ran - 0.5)*0.000001
      END DO
      ! second step is to find the new nearest neighbors
      PRINT *, truer
      PRINT *,"broken code executed"
      STOP
      unZeroHessian = 0
      RETURN
    END FUNCTION 

    
    ! This subroutine checks if PH rule is satisfied given crystal/molecule
    ! inport or not
    ! USED IN THIS MODULE
    SUBROUTINE PHRuleExam(ucpCounts,opts,ions,phmrCompliant)
      TYPE(options_obj) :: opts
      TYPE(ions_obj) :: ions
      INTEGER, DIMENSION(4) :: ucpCounts
      INTEGER :: phSum,iphsum
      LOGICAL :: phmrCompliant
      phSum = ucpCounts(1) - ucpCounts(2) + ucpCounts(3) - ucpCounts(4)
      IF (ucpCounts(1) /= ions%nions) PRINT *, "WARNING: There is more nuclear CP &
        than number of atoms!" 
      iphSum = ions%nions - ucpCounts(2) + ucpCounts(3) - ucpCounts(4) 
      !phSum = iphSum ! Using atom count to override
      phmrCompliant = .FALSE.
      IF (opts%isCrystal) THEN
        PRINT *, 'The system is assigned as a Crystal'
        IF (.NOT. (ucpCounts(1)>=1 .AND. ucpCounts(2)>=3 .AND. ucpCounts(3)>=3 .AND. ucpCounts(4)>=1)) THEN
          PRINT *,''//achar(27)//'[31m ERROR: FAILED Morse relationship. Number&
            of CPs found is not enough.' &
            //achar(27)//'[0m'
        ELSE IF (phSum == 0) THEN
          PRINT *, ''//achar(27)//'[32m Satisfies the Morse Relationship' &
            //achar(27)//'[0m'
          phmrCompliant = .TRUE.
        ELSE IF (phSum == 1 .AND. iphSum == phSum ) THEN
          PRINT *, ''//achar(27)//'[31m ERROR: The result satisfies the  & 
            Poincare Hopf Rule for a & molecule, not a crystal.' &
            //achar(27)//'[0m'
        ELSE IF (iphSum == 0 .AND. iphSum /= phSum) THEN
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          PRINT *, ''//achar(27)//'[32m Satisfies the Morse Relationship' &
            //achar(27)//'[0m'
          phmrCompliant = .TRUE.
        ELSE IF (iphSum == 1 .AND. iphSum /= phSum) THEN
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          PRINT *, ''//achar(27)//'[31m ERROR: The result satisfies the  & 
            Poincare Hopf Rule for a & molecule, not a crystal.' &
            //achar(27)//'[0m'
        ELSE 
          PRINT *, ''//achar(27)//'[31m ERROR: FAILED Morse relationship' &
            //achar(27)//'[0m'
        END IF
      ELSE IF (opts%isMolecule) THEN
        PRINT *, 'The system is assigned as a Molecule'
        IF (phSum == 0) THEN
          PRINT *, ''//achar(27)//'[31m ERROR: The result satisfies the Morse Relationship &
            for a crystal, not a molecule.'//char(27)//'[0m'
        ELSE IF (phSum == 1) THEN
          phmrCompliant = .TRUE.
          PRINT *, ''//achar(27)//'[32m Satisfies the Poincare Hopf Rule' &
            //achar(27)//'[0m'
        ELSE IF (iphSum == 0 .AND. iphSum /= phSum) THEN
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          PRINT *, ''//achar(27)//'[31m ERROR: The result satisfies the Morse Relationship &
            for a crystal, not a molecule.'//char(27)//'[0m'
        ELSE IF (iphSum == 1 .AND. iphSum /= phSum) THEN
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          phmrCompliant = .TRUE.
          PRINT *, ''//achar(27)//'[32m Satisfies the Poincare Hopf Rule' &
            //achar(27)//'[0m'
        ELSE 
          PRINT *, ''//achar(27)//'[31m ERROR: FAILED Poincare Hopf Rule' &
            //achar(27)//'[0m'
        END IF
      ELSE 
        IF (phSum == 0) THEN
          IF (.NOT. (ucpCounts(1)>=1 .AND. ucpCounts(2)>=3 .AND. ucpCounts(3)>=3 .AND. ucpCounts(4)>=1)) THEN
            PRINT *,''//achar(27)//'[31m ERROR: FAILED Morse relationship. Number&
              of CPs found is not enough and FAILED Poincare Hopf Rule.' &
            //achar(27)//'[0m'
            phmrCompliant = .FALSE.
          ELSE
            phmrCompliant = .TRUE.
            PRINT *, ''//achar(27)//'[32m This system has not been designated & 
              as a molecule or crystal but the Morse relationship for & 
              crystals are satisfied.' //achar(27)//'[0m'
          END IF
        ELSE IF (phSum == 1) THEN
          phmrCompliant = .TRUE.
          PRINT *, ''//achar(27)//'[32m This system has not been designated & 
            as a molecule or crystal but the Poincare-Hopf rule for &
            molecules are satisfied.' //achar(27)//'[0m'
        END IF
        IF (iphSum == 0 .AND. iphSum /= phSum) THEN
          phmrCompliant = .TRUE.
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          PRINT *, ''//achar(27)//'[32m This system has not been designated & 
            as a molecule or crystal but the Morse relationship for & 
            crystals are satisfied.' //achar(27)//'[0m'
        ELSE IF (iphSum == 1 .AND. iphSum /= phSum) THEN
          phmrCompliant = .TRUE.
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          PRINT *, ''//achar(27)//'[32m This system has not been designated & 
            as a molecule or crystal but the Poincare-Hopf rule for &
            molecules are satisfied.' //achar(27)//'[0m'
        END IF
        IF (.NOT. phmrCompliant) THEN
          PRINT *, ''//achar(27)//'[31m ERROR: FAILED Poincare Hopf Rule & 
            and Morse Relationship' //achar(27)//'[0m'
        END IF
      END IF
    END SUBROUTINE PHRuleExam

    ! count the number of negative eigenvalues to characterize a critical point
    ! the version of the above subroutine where p is real not integer
    ! USED IN THIS MODULE
    SUBROUTINE RecordCPR(p,chg,cpl,ucptnum,connectedAtoms, ucpCounts, &
      opts,hessianMatrix,ind)
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      REAL(q2), DIMENSION(3) :: eigvals, grad
      REAL(q2), DIMENSION(3,3) :: eigvecs, hessianMatrix
      INTEGER, DIMENSION(4) :: ucpCounts
      INTEGER, DIMENSION(2) :: connectedAtoms
      REAL(q2), DIMENSION(3) :: p,realdump
      REAL(q2) :: rho
      INTEGER,DIMENSION(3) :: ind
      INTEGER :: ucptnum, negCount
      INTEGER :: it_num, rot_num
      LOGICAL :: LDM
      cpl(ucptnum)%hessianMatrix = hessianMatrix
      CALL jacobi_eigenvalue(3,hessianMatrix,9999,eigvecs,eigvals,&
        it_num, rot_num )
      negCount = CountNegModes(eigvals)
      IF (LDM) THEN
        PRINT *, "CountNegModes counted negCount as "
        PRINT *, negCount
      END IF
      CALL UpDateCounts(negCount,ucpCounts)
      cpl(ucptnum)%truer = p
      !cpl(ucptnum)%grad = grad !where is this getting value from?
      cpl(ucptnum)%eigvecs = eigvecs
      cpl(ucptnum)%eigvals = eigvals
      cpl(ucptnum)%negcount = negcount
      cpl(ucptnum)%hasProxy = .FALSE.
      cpl(ucptnum)%ind = ind
      cpl(ucptnum)%connectedAtoms = connectedAtoms
      realDump = rho_grad(chg,p,rho) 
      rho = rho 
      cpl(ucptnum)%rho = rho
      IF (LDM) THEN
        PRINT *, "RecordCPR recorded CP number ", ucptnum
        PRINT *, "location is "
        PRINT *, p
        PRINT *, "eigvals is "
        PRINT *, eigvals
        PRINT *, "negcount is "
        PRINT *, cpl(ucptnum)%negCount
      END IF
    END SUBROUTINE RecordCPR
   
    ! This is for recording cage points found through density descend. 
    ! These cage points are on grid points. They can have a gradient 
    ! and a hessian but the hessian may not have 3 positive eigenvalues.
    SUBROUTINE RecordDensityDescendCage(p,cpl,ucptnum,ucpCounts,chg)
      TYPE(charge_obj) :: chg
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      INTEGER, DIMENSION(3) :: p
      INTEGER, DIMENSION(4) :: ucpCounts
      INTEGER :: ucptnum, it_num, rot_num
      cpl(ucptnum)%grad = CDGrad(p,chg)
      cpl(ucptnum)%hessianMatrix = CDHessian(p,chg)
      CALL jacobi_eigenvalue(3,cpl(ucptnum)%hessianMatrix, &
        9999,cpl(ucptnum)%eigvecs,cpl(ucptnum)%eigvals,it_num,rot_num) 
      cpl(ucptnum)%negCount = 0
      CALL UpDateCounts(0,ucpCounts)
      cpl(ucptnum)%trueR = p
      cpl(ucptnum)%rho = chg%rho(p(1),p(2),p(3))
      cpl(ucptnum)%hasProxy = .FALSE.
    END SUBROUTINE RecordDensityDescendCage

    ! USED IN THIS MODULE
    SUBROUTINE RecordCPRLight(p,chg,cpl,ucptnum, ucpCounts, &
      ind,LDM)
      TYPE(charge_obj) :: chg
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      REAL(q2), DIMENSION(3) :: eigvals, grad
      REAL(q2), DIMENSION(3,3) :: eigvecs, hessianMatrix
      REAL(q2), DIMENSION(3) :: p, realDump
      REAL(q2) :: rho
      INTEGER, DIMENSION(4) :: ucpCounts
      INTEGER,DIMENSION(3) :: ind
      INTEGER :: ucptnum, negCount
      INTEGER :: it_num, rot_num
      LOGICAL :: LDM
      hessianMatrix = CDHessianR(p,chg)
      grad = CDGradR(p,chg)
      cpl(ucptnum)%hessianMatrix = hessianMatrix
      !CALL DSYEVJ3(hessianMatrix,eigvecs,eigvals)
      CALL jacobi_eigenvalue(3,hessianMatrix,9999,eigvecs,eigvals,&
        it_num, rot_num )
      negCount = CountNegModes(eigvals)
      CALL UpDateCounts(negCount,ucpCounts)
      cpl(ucptnum)%truer = p
      cpl(ucptnum)%grad = grad
      cpl(ucptnum)%eigvecs = eigvecs
      cpl(ucptnum)%eigvals = eigvals
      cpl(ucptnum)%negCount = negCount
      cpl(ucptnum)%hasProxy = .FALSE.
      cpl(ucptnum)%ind = ind
      realDump = rho_grad(chg,p,rho)
      cpl(ucptnum)%rho = rho
      IF (LDM) THEN
        PRINT *, "RecordCPRLight recorded CP number ",ucptnum
        PRINT *, "truer is"
        PRINT *, p
        PRINT *, "eigvals is"
        PRINT *, eigvals
        PRINT *, "negCount is"
        PRINT *, negCount
        PRINT *, "Hessian matrix is "
        PRINT *, hessianMatrix(1,:)
        PRINT *, hessianMatrix(2,:)
        PRINT *, hessianMatrix(3,:)
      END IF
    END SUBROUTINE RecordCPRLight

    SUBROUTINE RecordCPStatic(i,cp_static, reduced_cp_static, reduced_ucptnum, reduced_ucpCounts)
      TYPE(cpc), ALLOCATABLE, DIMENSION(:) :: cp_static, reduced_cp_static
      INTEGER, DIMENSION(4) :: reduced_ucpCounts
      INTEGER :: reduced_ucptnum, i
     
      reduced_cp_static(reduced_ucptnum) = cp_static(i)
      CALL UpDateCounts(reduced_cp_static(i)%negCount,reduced_ucpCounts)
    END SUBROUTINE RecordCPStatic

    ! USED IN THIS MODULE
    SUBROUTINE  OutputCP(cpl,opts,ucptnum,chg,setcount, ucpCounts)
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      TYPE(options_obj) :: opts
      TYPE(charge_obj) :: chg
      INTEGER, DIMENSION(4) :: ucpCounts
      INTEGER :: ucptnum, i, setcount
      ! TYPE(ions_obj) :: ions
      ! vol was used to calculate charge density with cpl(i)%rho/vol
      ! REAL(q2) :: vol
      ! vol = matrix_volume(ions%lattice)
      CHARACTER(10) :: fileName
      PRINT *, 'Writting critical point output files'
      WRITE(fileName,fmt='(a,i2.2,a)') TRIM('CPFU'), setcount,TRIM('.dat')
      PRINT *, 'Critical point information are written in file: ', filename
      OPEN(98,FILE=filename,STATUS='REPLACE',ACTION='WRITE')
      WRITE(98,*) 'Number of critical points written: ', ucptnum
      WRITE(98,*) 'Number of bond, ring and cage critical point written: ', ucpCounts(2:4)
      WRITE(98,*) 'Nuclear critical points are omitted. Use atomic positions.'
      DO i = 1, ucptnum
        !IF (cpl(i)%negCount == 3 ) CYCLE
        WRITE (98,*) '_______________________________________'
        WRITE (98,*) 'Unique CP # : ', i
        WRITE (98,*) 'Critical point is found at indices'
        WRITE (98,*) cpl(i)%truer
        WRITE (98,*) 'Coordinates in cartesian are'
        WRITE (98,*) MATMUL(chg%lat2car,cpl(i)%truer)
        WRITE (98,*) 'Direct coordinates are'
        WRITE (98,*) cpl(i)%truer(1)/chg%npts(1), &
          cpl(i)%truer(2)/chg%npts(2), &
          cpl(i)%truer(3)/chg%npts(3)
        WRITE (98,*) "Charge density is"
        ! It appears rho_val values are densities indeed. Volume is already
        ! factored in.
        WRITE (98,*) cpl(i)%rho
        WRITE (98,*) 'Gradient is'
        WRITE (98,*) cpl(i)%grad
        WRITE (98,*) 'Hessian is'
        WRITE (98,*) cpl(i)%hessianMatrix(1,:)
        WRITE (98,*) cpl(i)%hessianMatrix(2,:)
        WRITE (98,*) cpl(i)%hessianMatrix(3,:)
        WRITE (98,*) 'Laplacian is'
        WRITE (98,*) cpl(i)%hessianMatrix(1,1)**2 &
                     + cpl(i)%hessianMatrix(2,2)**2 &
                     + cpl(i)%hessianMatrix(3,3)**2
        WRITE (98,*) 'Eigenvalues are'
        WRITE (98,*) cpl(i)%eigvals
        WRITE (98,*) 'Eigenvectors are'
        WRITE (98,*) cpl(i)%eigvecs(1,:)
        WRITE (98,*) cpl(i)%eigvecs(2,:)
        WRITE (98,*) cpl(i)%eigvecs(3,:)

        IF (cpl(i)%negcount == 2) THEN
          WRITE(98,*) 'Connected Atoms are'
          WRITE(98,*) cpl(i)%connectedAtoms
        END IF
        IF (cpl(i)%negcount == 0) THEN
          WRITE(98,*) 'This is a cage critical point'
          WRITE(98,*) ' '
        END IF
        IF (cpl(i)%negcount == 2) THEN
          WRITE(98,*) 'This is a bond critical point'
          WRITE(98,*) ' '
        ELSEIF(cpl(i)%negcount == 1) THEN
          WRITE(98,*) 'This is a ring critical point'
          WRITE(98,*) ' '
        ELSEIF(cpl(i)%negcount == 3) THEN
          WRITE(98,*) 'This is a nuclear critical point'
          WRITE(98,*) ' ' 
        END IF
        WRITE(98,*) '_________________________________________'
      END DO
    WRITE(98,*) ''
    CLOSE(98)
    END SUBROUTINE OutputCP
    
    !Outputs a list of connections in a file
    ! USED IN THIS MODULE
    SUBROUTINE OutputNetwork(cpl,ucptnum,setcount)
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      INTEGER :: ucptnum, i, setcount
      CHARACTER(10) :: fileName
      INTEGER, DIMENSION(2) :: pair
      WRITE(fileName,fmt='(a,i2.2,a)') TRIM('LINK'),setcount, TRIM('.dat')
      PRINT *, "Network Graph is written in file: ", fileName
      OPEN(98,FILE=fileName,STATUS='REPLACE',ACTION='WRITE')
      !WRITE(98,*) size(ions%atomic_num)
      !DO i=1,size(ions%atomic_num)
      !  WRITE (98,*) ions%atomic_num(i)  
      !END DO
      WRITE(98,*) "Connections between bond CP and atoms"
      WRITE(98,*) "Connection 1 , ", " CP # , "," Connection 2"
      DO i=1,ucptnum
        IF (cpl(i)%negCount /= 2) THEN
          CYCLE
        END IF
        pair = cpl(i)%connectedAtoms
        IF (pair(1) == 0 .AND. pair(2) == 0 ) THEN
          !PRINT *, "Invalid bond CP"
        ELSE
          WRITE(98,*) pair(1),' ', i, ' ', pair(2)
        END IF
      END DO
      WRITE(98,*) "Connections between ring CP and cage CP"
      WRITE(98,*) "connection to -1 means the connected cage point is located in vacuum"
      WRITE(98,*) "connection to -2 means the connected cage point is not found"
      WRITE(98,*) "Connection 1 , ", " CP # , "," Connection 2"
      DO i=1,ucptnum
        IF (cpl(i)%negCount /= 1) THEN
          CYCLE
        END IF
        pair = cpl(i)%connectedAtoms
        IF (pair(1) == 0 .AND. pair(2) == 0 ) THEN
          !PRINT *, "Invalid bond CP"
        ELSE
          WRITE(98,*) pair(1),' ', i, ' ', pair(2)
        END IF
      END DO
      ! Not sure if we need the connectivity matrix
      !WRITE(98,*) ""
      !WRITE(98,*) "The connectivity matrix is"
      !DO i = 1,size(atom_connectivity,1)
      !  WRITE(98,*) atom_connectivity(i,:)
      !END DO
      CLOSE(98)
    END SUBROUTINE OutputNetwork

    FUNCTION CheckIsolatedAtom(cpl,ucptnum)
      TYPE(cpc),DIMENSION(:),ALLOCATABLE :: cpl
      INTEGER :: i,j, ucptnum, connection_count
      LOGICAL :: CheckIsolatedAtom
      
      CheckIsolatedAtom = .TRUE.
      DO i = 1, ucptnum
        connection_count = 0
        IF (cpl(i)%negCount == 3) THEN
          DO j = 1, ucptnum
            IF (cpl(j)%negCount == 2) THEN
              IF (ANY(cpl(j)%connectedAtoms == i)) connection_count = connection_count + 1 
            END IF
          END DO
          IF (connection_count < 1) THEN
            CheckIsolatedAtom = .FALSE.
            PRINT *, "ERROR: cp # ", i, "has no connection to any bond CP"
          END IF
        END IF
      END DO

      RETURN 
    END FUNCTION CheckIsolatedAtom

    ! TRUE means there is no problem.
    FUNCTION CheckIsolatedRing(cpl,ucptnum)
      TYPE(cpc),DIMENSION(:),ALLOCATABLE :: cpl
      INTEGER :: n, ucptnum
      LOGICAL :: CheckIsolatedRing, voidDweller
      CheckIsolatedRing = .FALSE.  
      voidDweller = .FALSE. !if there are rings connecting to vacuum
      DO n=1,ucptnum
        IF (cpl(n)%negCount == 1) THEN
          IF (CheckIsolatedRing) THEN
            CYCLE
          ELSE
            CheckIsolatedRing = ALL(cpl(n)%connectedAtoms /= -2 ) .AND. ALL(cpl(n)%connectedAtoms /= 0) 
          END IF
          IF (voidDweller) THEN
            CYCLE
          ELSE
            voidDweller = voidDweller .OR. ANY(cpl(n)%connectedAtoms == -1)
          END IF
        END IF
      END DO
      IF (voidDweller) PRINT *, "WARNING: at least one ring point is expected in the vacuum. &
        Examine the results carefully. This may be OK for non-crystaline system."
      RETURN
    END FUNCTION CheckIsolatedRing
   
    ! TRUE means there is no problem.
    FUNCTION CheckIsolatedCage(cpl,ucptnum)
      TYPE(cpc),DIMENSION(:),ALLOCATABLE :: cpl
      INTEGER :: i, j, ucptnum, connection_count
      LOGICAL :: CheckIsolatedCage
      CheckIsolatedCage = .TRUE.  
      DO i = 1, ucptnum
        connection_count = 0
        IF (cpl(i)%negCount == 0) THEN
          DO j = 1, ucptnum
            IF (cpl(j)%negCount == 1) THEN
              IF (ANY(cpl(j)%connectedAtoms == i)) connection_count = connection_count + 1 
            END IF
          END DO
          IF (connection_count < 1) THEN
            CheckIsolatedCage = .FALSE.
            PRINT *, "ERROR: cp # ", i, "has no connection to any ring CPs" ! todo: check if it's connecting to a periodic self.
          END IF
        END IF
      END DO
      RETURN
    END FUNCTION CheckIsolatedCage


    ! USED IN THIS MODULE
    SUBROUTINE OutputCPRoster(cpRoster,setcount)
      REAL(q2), DIMENSION(:,:), ALLOCATABLE :: cpRoster
      CHARACTER(19) :: fileName
      INTEGER :: i,setcount
      WRITE(fileName,fmt='(a,i2.2,a)') TRIM ('TrajEndPoints'), setcount, TRIM('.dat')
      PRINT *, "Trajectory End Point Roster is written in file: ", fileName 
      OPEN(98,FILE=fileName,STATUS='REPLACE',ACTION='WRITE')

      DO i=1,size(cpRoster,1)
        WRITE(98,*) cpRoster(i,:)
      END DO
      CLOSE(98)
    END SUBROUTINE OutputCPRoster
   
    ! USED IN THIS MODULE
    SUBROUTINE CPRosterAnalysis(cpl,ions, cproster, chg)
      TYPE(cpc), ALLOCATABLE,DIMENSION(:) :: cpl
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      REAL(q2), DIMENSION(3) :: r, cp, nnp 
      REAL(q2), DIMENSION(:,:), ALLOCATABLE :: cproster
      REAL(q2), DIMENSION(:), ALLOCATABLE :: nearby
      LOGICAL, DIMENSION(size(cproster,1)) :: mask
      REAL(q2), DIMENSION(size(cproster,1)) :: distances
      INTEGER :: i,trajcount, ij,i1,i2,i3
      REAL(q2) :: val, sumsq, StdDev, truedist
      REAL(q2), DIMENSION(size(cpl)) :: StdDevList
      
      trajcount = 0
      DO ij = 1,size(cpl)
        cp = cpl(ij)%truer
        mask = .FALSE.
        DO i=1,size(cproster,1)
          r = cproster(i,:)
          val = Mag(r-cp)
          IF (val < 1.0) THEN
            mask(i) = .TRUE.
          END IF
          distances(i) = val
        END DO

        nearby = PACK(distances,mask)
        sumsq= 0
        trajcount = trajcount + size(nearby)
        DO i =1, size(nearby)
          sumsq = sumsq + nearby(i)**2
        END DO
        IF (size(nearby) == 0 ) THEN
          StdDev = 0.00
        ELSE
          StdDev = SQRT(sumsq/(size(nearby)) )
        END IF
        
        StdDevList(ij) = StdDev

        PRINT "(I5,A,es10.3,A,i5,A)", ij, ": StdDev is ", StdDev, " with ", size(nearby), " trajectories"
      !   PRINT *, cp
        IF (StdDev > 0.01) THEN
          PRINT *, "High Standard Deviation detected."
          PRINT *, " " 
        END IF
        DEALLOCATE(nearby)
      END DO
      PRINT *, "Total Trajectories Used: ", trajcount
      
      DO ij = 1,size(cpl)
         cp = cpl(ij)%truer 
         DO i=ij,size(cpl)
           IF (i /= ij) THEN
             truedist = Mag(cpl(i)%truer - cp)
             DO i1 = -1,1
               DO i2 = -1,1
                 DO i3 = -1,1
                   IF (i1 == 0 .AND. i2 == 0 .AND. i3 == 0) CYCLE
                     nnp = cp + (/i1 * chg%npts(1), i2 * chg%npts(2), i3 * chg%npts(3)/)
                     truedist = MIN(truedist, Mag(cpl(i)%truer - nnp) )
                 END DO
               END DO
             END DO
             
             
             
             IF (truedist < 5) THEN
               PRINT "(A,i5,A,i5,A,es10.3)", "Distance between", ij, " and ", i, " is: ", truedist
                
               IF (cpl(ij)%negcount /= cpl(i)%negcount) THEN
                 IF (cpl(ij)%negcount == 0) THEN
                   PRINT "(i5,A)", ij," is of type: nuclear"
                 ELSE IF (cpl(ij)%negcount == 1) THEN
                   PRINT "(i5,A)", ij," is of type: bond"
                 ELSE IF (cpl(ij)%negcount == 2) THEN
                   PRINT "(i5,A)", ij, " is of type: ring"
                 ELSE
                   PRINT "(i5,A)", ij," is of type: cage"
                 END IF

                 IF (cpl(i)%negcount == 0) THEN
                    PRINT "(i5,A)", i," is of type: nuclear"
                 ELSE IF (cpl(i)%negcount == 1) THEN
                    PRINT "(i5,A)", i," is of type: bond"
                 ELSE IF (cpl(i)%negcount == 2) THEN
                   PRINT "(i5,A)", i, " is of type: ring"
                 ELSE
                   PRINT "(i5,A)", i," is of type: cage"
                 END IF
                 PRINT *, " "
               END IF
             END IF
           END IF
         END DO

      END DO

    END SUBROUTINE CPRosterAnalysis
    
    
    ! get the distance between two points, considering that these two points can
    ! be neighbors across the periodic boundary
    ! The output distance is in cartesian
    FUNCTION GetPointDistance(p, pt, chg, nnlayer)
      REAL(q2) :: GetPointDistance, distance
      ! these are indices
      INTEGER, DIMENSION(3) :: p, pt
      INTEGER, DIMENSION( 3) :: nnp
      ! these are cartesian coordinates
      REAL(q2), DIMENSION(3) :: r, rt
      REAL(q2), DIMENSION(3) :: nnr
      TYPE(charge_obj) :: chg
      INTEGER :: nnlayer, i, j, k
      ! look for all equivalents of p within nnlayer, find the smallest distance
      rt = MATMUL(chg%lat2car,pt)
      r = MATMUL(chg%lat2car,p)
      GetPointDistance = Mag( r - rt)
      DO i = -nnlayer, nnlayer
        DO j = -nnlayer, nnlayer
          DO k = -nnlayer, nnlayer
            IF (i == 0.AND.j == 0.AND.k == 0) CYCLE
            nnp = p + (/i * chg%npts(1), j * chg%npts(2), k * chg%npts(3)/)
            nnr = MATMUL(chg%lat2car,nnp)
            distance = Mag( nnr - rt)
            GetPointDistance = MIN(GetPointDistance, distance)
          END DO
        END DO
      END DO
      RETURN
    END FUNCTION GetPointDistance

    ! This function takes two fractional lattice coordinate sets and produces
    ! the distance in cartesian
    FUNCTION GetPointDistanceR(p, pt, chg, nnlayer)
      REAL(q2) :: GetPointDistanceR, distance
      ! these are indices
      REAL(q2), DIMENSION(3) :: p, pt
      INTEGER, DIMENSION( 3) :: nnp
      ! these are cartesian coordinates
      REAL(q2), DIMENSION(3) :: r, rt
      REAL(q2), DIMENSION(3) :: nnr
      TYPE(charge_obj) :: chg
      INTEGER :: nnlayer, i, j, k
      ! look for all equivalents of p within nnlayer, find the smallest distance
      rt = MATMUL(chg%lat2car,pt)
      r = MATMUL(chg%lat2car,p)
      GetPointDistanceR = Mag( r - rt)
      DO i = -nnlayer, nnlayer
        DO j = -nnlayer, nnlayer
          DO k = -nnlayer, nnlayer
            IF (i == 0.AND.j == 0.AND.k == 0) CYCLE
            nnp = p + (/REAL(i * chg%npts(1),q2), REAL(j * chg%npts(2),q2), &
                  REAL(k * chg%npts(3),q2)/)
            nnr = MATMUL(chg%lat2car,nnp)
            distance = Mag( nnr - rt)
            GetPointDistanceR = MIN(GetPointDistanceR, distance)
          END DO
        END DO
      END DO
      RETURN
    END FUNCTION GetPointDistanceR


    ! USED IN THIS MODULE
    SUBROUTINE MakeCPRoster(cpr,cptnum,r)
      REAL(q2), DIMENSION(:,:), ALLOCATABLE :: cpr
      INTEGER :: cptnum
      REAL(q2), DIMENSION(3) :: r
      cpr(cptnum,:) = r
    END SUBROUTINE MakeCPRoster
    
    !Stores every coordinate converged to, even if nonunique
    ! USED IN THIS MODULE
    SUBROUTINE MakeFullCPRoster(cpr,cptnum,r)
      REAL(q2), DIMENSION(:,:), ALLOCATABLE :: cpr
      INTEGER :: cptnum
      REAL(q2), DIMENSION(3) :: r
      cpr(cptnum, :) = r
    END SUBROUTINE MakeFullCPRoster

    ! this function interpolates the gradient of a point. weight towards each nn
    ! is determined by its distance to that neighbor. 
    FUNCTION R2GradInterpol(nnInd,r,chg,nnlayers)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3) :: r, R2GradInterpol
      REAL(q2),DIMENSION(:),ALLOCATABLE :: weight
      REAL(q2) :: normalizer,distance
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: nnInd
      INTEGER :: nnlayers, i
      LOGICAL :: onGrid
      normalizer = 0
      onGrid = .FALSE.
      ALLOCATE(weight(SIZE(nnInd)/3))
      DO i = 1,SIZE(nnInd)/3
        distance = GetPointDistanceR( &
          MATMUL(chg%lat2car,r),MATMUL(chg%lat2car,nnInd(i,:)) &
          ,chg,nnlayers)
        IF (distance == 0) THEN
          onGrid = .TRUE.
          EXIT
        END IF
        weight(i) = (1/distance )**2
        normalizer = normalizer + weight(i)
      END DO   
      IF (ongrid) THEN
        R2GradInterpol = CDGrad(nnind(i,:),chg)
      ELSE
        weight = weight / normalizer
        R2GradInterpol = 0
        DO i = 1,SIZE(nnInd)/3
          R2GradInterpol = R2GradInterpol + CDGrad(nnInd(i,:),chg) * weight(i) 
        END DO
      END IF
      DEALLOCATE(weight)
      RETURN
    END FUNCTION R2GradInterpol
   
    !FUNCTION R2HesInterpol(nnInd,r,chg,nnLayers)
    FUNCTION R2HesInterpol(nnInd,r,chg)  
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3,3) :: R2HesInterpol
      REAL(q2),DIMENSION(:),ALLOCATABLE :: weight
      REAL(q2),DIMENSION(3) :: r
      REAL :: normalizer, distance
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: nnInd
      INTEGER :: nnlayers, i
      LOGICAL :: onGrid
      normalizer = 0
      onGrid = .FALSE.
      ALLOCATE(weight(SIZE(nnInd)/3))
      DO i = 1,SIZE(nnInd)/3
        distance = GetPointDistanceR( &
          MATMUL(chg%lat2car,r),MATMUL(chg%lat2car,nnInd(i,:)) &
          ,chg,nnlayers)
        IF (distance == 0) THEN
          onGrid = .TRUE.
          EXIT
        END IF
        weight(i) =  (1/distance)**2
        normalizer = normalizer + weight(i)
      END DO
      IF (onGrid) THEN
        R2HesInterpol = CDHessian(nnInd(i,:),chg)
      ELSE
        weight = weight / normalizer
        R2HesInterpol = 0
        DO i = 1,SIZE(nnInd)/3
          R2HesInterpol = R2HesInterpol + CDHessian(nnInd(i,:),chg) * weight(i)
        END DO
      END IF
      DEALLOCATE(weight)
      RETURN
    END FUNCTION R2HesInterpol

    !FUNCTION R2RhoInterpol(nnInd,r,chg,nnLayers)
    FUNCTION R2RhoInterpol(nnInd,r,chg)  
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(:),ALLOCATABLE :: weight
      REAL(q2),DIMENSION(3) :: r
      REAL :: normalizer, R2RhoInterpol, distance
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: nnInd
      INTEGER :: nnlayers, i
      LOGICAL :: onGrid
      normalizer = 0
      onGrid = .FALSE.
      ALLOCATE(weight(SIZE(nnInd)/3))
      DO i = 1,SIZE(nnInd)/3
        distance = GetPointDistanceR( &
          MATMUL(chg%lat2car,r),MATMUL(chg%lat2car,nnInd(i,:)) &
          ,chg,nnlayers)
        IF (distance == 0) THEN
          onGrid = .TRUE.
          EXIT
        END IF
        weight(i) = (1/distance)**2
        normalizer = normalizer + weight(i)
      END DO
      IF (onGrid) THEN
        R2RhoInterpol = rho_val(chg,nnind(i,1),nnind(i,2),nnind(i,3))
      ELSE
        weight = weight / normalizer
        R2RhoInterpol = 0
        DO i = 1,SIZE(nnInd)/3
          R2RhoInterpol = R2RhoInterpol + rho_val(chg,nnind(i,1),nnind(i,2),&
            nnind(i,3)) * weight(i)
        END DO
      END IF
      DEALLOCATE(weight)
      RETURN
    END FUNCTION R2RhoInterpol
 
    ! updates the count on all types of CPs
    ! USED IN THIS MODULE
    SUBROUTINE UpDateCounts(negCount,ucpCounts)
      INTEGER, DIMENSION(4) :: ucpCounts
      INTEGER :: negCount
      IF (negCount == 3) THEN
        ucpCounts(1) = ucpCounts(1) + 1
      ELSE IF (negCount == 2) THEN
        ucpCounts(2) = ucpCounts(2) + 1
      ELSE IF (negCount == 1) THEN
        ucpCounts(3) = ucpCounts(3) + 1
      ELSE 
        ucpCounts(4) = ucpCounts(4) + 1
      END IF
    END SUBROUTINE UpDateCounts

    ! This is a subroutine for debugging purpose. It takes all converged
    ! locations of all CP candidates, and outputs them. Bond CP are marked as He,
    ! Ring CP are marked as Ne, cage CP are marked as Ar. Nucleus are written as
    ! normal
    ! USED IN THIS MODULE
    SUBROUTINE VisAllCP(cpcl,cptnum,chg,ions,opts,&
      ucpCounts)
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      TYPE(options_obj) :: opts
      TYPE(cpc),DIMENSION(:),ALLOCATABLE :: cpcl
      !REAL(q2),DIMENSION(8,3) :: nnGrad
      REAL(q2),DIMENSION(3) :: r
      INTEGER, DIMENSION(4) :: ucpCounts
      INTEGER :: cptnum,n1,j
      CHARACTER(LEN=128) :: atoms, natoms
      OPEN(100,FILE=opts%chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
      IF(opts%in_opt == opts%in_chgcar5) THEN
        DO n1 = 1, 5
          READ(100,'(a)')
        END DO
        READ(100,'(a)') atoms
        READ(100,'(a)') natoms
      ELSE
        ! have to hope that somewhere in the CHGCAR elements are specified.
        READ(100,'(a)') atoms
        DO n1 = 1,4
          READ(100,'(a)')
          ! skipping 1 line of scaling factor and 3 lines of lattice coordinates
        END DO
        READ(100,'(a)') natoms
      END IF
      CLOSE(100)
      ! Write header of allcpPOSCAR
      IF (opts%gradMode) THEN
         OPEN(11,FILE='allcpPOSCAR_GD',STATUS='REPLACE',ACTION='WRITE')
      ELSE
         OPEN(11,FILE='allcpPOSCAR_NM',STATUS='REPLACE',ACTION='WRITE')
      END IF
      IF (ucpCounts(4)>0) WRITE(11,'(a)',ADVANCE='NO') 'Ar '
      IF (ucpCounts(3)>0) WRITE(11,'(a)',ADVANCE='NO') 'Ne '
      IF (ucpCounts(2)>0) WRITE(11,'(a)',ADVANCE='NO') 'He '
      WRITE(11,'(a)') TRIM(atoms)
      WRITE(11,*) ions%scalefactor
      WRITE(11,*) '     ',ions%lattice(1,:)
      WRITE(11,*) '     ',ions%lattice(2,:)
      WRITE(11,*) '     ',ions%lattice(3,:)
      IF (ucpCounts(4)>0) WRITE(11,'(a)',ADVANCE='NO') 'Ar '
      IF (ucpCounts(3)>0) WRITE(11,'(a)',ADVANCE='NO') 'Ne '
      IF (ucpCounts(2)>0) WRITE(11,'(a)',ADVANCE='NO') 'He '
      WRITE(11,'(a)') TRIM(atoms)
      IF (ucpCounts(4)>0) WRITE(11,'(I4)',ADVANCE='NO') ucpCounts(4)
      IF (ucpCounts(3)>0) WRITE(11,'(I4)',ADVANCE='NO') ucpCounts(3)
      IF (ucpCounts(2)>0) WRITE(11,'(I4)',ADVANCE='NO') ucpCounts(2)
      WRITE(11,'(a)') TRIM(natoms)
      WRITE(11,'(a)') 'Cartesian'
      ! Loop three times, each time writes out one type of CP
      DO j = 0, 2
        DO n1 = 1, cptnum
          IF (.NOT.cpcl(n1)%isunique) CYCLE
          IF (cpcl(n1)%negCount == j) THEN
            ! Because lattice starts at 1 1 1 but cartesian starts at 0 0 0
            r(1) = cpcl(n1)%trueR(1) - 1
            r(2) = cpcl(n1)%trueR(2) - 1
            r(3) = cpcl(n1)%trueR(3) - 1
            WRITE(11,*) MATMUL(chg%lat2car,r)
           ! WRITE(11,*) MATMUL(chg%lat2car,cpcl(n1)%trueR)
           ! dir(1) = cpcl(n1)%trueR(1)/chg%npts(1)
           ! dir(2) = cpcl(n1)%trueR(2)/chg%npts(2)
           ! dir(3) = cpcl(n1)%trueR(3)/chg%npts(3)
          END IF
        END DO
      END DO
      PRINT *, 'writting atomic locations'
      DO j = 1, ions%nions
        WRITE(11,*) ions%r_car(j,:)
      END DO
      CLOSE(11)
    END SUBROUTINE VisAllCP
  
    ! This subroutine tracks steps in up to the past 10 steps. If the next step is
    ! identical to one taken before, it gives the location when repeat is
    ! detected 
    ! it could also give averaged location of the past
    ! 10 steps in ther future, given treatments to PBC.
    ! USED IN THIS MODULE
    SUBROUTINE DetectCircling(stepCount,rList,temList,trueR,nextTem,averageR,LDM,ind)
      REAL(q2),DIMENSION(10,3) :: rList,temList
      REAL(q2),DIMENSION(3) :: trueR,nextTem,averageR
      INTEGER,DIMENSION(3) :: ind
      INTEGER :: stepCount,i,j
      LOGICAL :: isRunningCircles,LDM
      ! establish lists if stepCount is low
      isRunningCircles = .FALSE.
      IF ( stepCount <= 10 ) THEN
        rList(stepCount,:) = trueR
        temList(stepCount,:) = nextTem
      ELSE 
        ! update lists
        DO i = 1, 9
          rList(i,:) = rList(i+1,:)
          temList(i,:) = temList(i+1,:)
        END DO
        rList(10,:) = trueR
        temList(10,:) = nextTem
        ! Detect if tem has repeated after 10 steps
        outer: DO i = 1, 10
          DO j = 1, 10
            IF ( j == i ) CYCLE
            IF (ALL(temList(i,:) == temList(j,:),1) ) THEN 
              IF (Mag(temList(j,:))== 1.) CYCLE
              isRunningCircles = .TRUE.
              EXIT outer
            END IF
          END DO
        END DO outer
        IF (isRunningCircles) THEN
          IF (LDM) THEN
            PRINT *, "De Bugger: Circling Detected"
            PRINT *, 'This trajetory initiated at ', ind
            PRINT *, "The past ten tems are"
            DO j = 1, 10
              PRINT *, temList(j,:)
            END DO
            PRINT *, "The past ten trueR are"
            DO j = 1, 10
              PRINT *, rList(j,:) 
            END DO
          END IF
          averageR = 0.
          averageR = trueR
! need to figure out how to deal with PBC first
!          DO i = 1,10
!            averageR = averageR + rList(i,:)
!          END DO
!          averageR = averageR/10
        END IF
      END IF
    END SUBROUTINE DetectCircling

    ! This subroutine looks for critical points in the list that is too close to
    ! another, and averages the same types into one to remove duplicate critical
    ! points.
    ! USED IN THIS MODULE
    SUBROUTINE ReduceCP(cpl,opts,ucptnum,chg,ucpCounts, &
      isReduced,LDM_RecordCPRLight,LDM,isReducible)
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl,rcpl
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      REAL(q2),DIMENSION(3) :: avgR
      INTEGER, DIMENSION(4) :: ucpCounts,rcpCounts
      INTEGER :: uCPTNum
      INTEGER :: i,j, nUCPTNum, weight, dupCount
      LOGICAL :: isReduced,LDM_RecordCPRLight,LDM, switch, isReducible
      IF (LDM) PRINT *, 'Checking for duplicate CP'
      isReduced = .TRUE.
      dupCount = 0
      rcpCounts = 0
      ! Give the reduced list same length as before, it's ok if a little goes to
      ! waste
      ALLOCATE(rcpl(SIZE(cpl)))
      ! start from a critical point. loop through the entire list, look for
      ! another entry 
      nUCPTNum = 0
      OUTER: DO i = 1, ucptnum
        ! check if this point is already determined as a proxy to some other
        ! point before
        weight = 1
        avgR = cpl(i)%truer
        switch = .TRUE.
        IF (cpl(i)%hasProxy ) THEN 
          CYCLE
        END IF
        DO j = i, ucptnum
          ! Periodic boundary condition will come to haunt me, periodically. 
          IF (j == i) CYCLE
          IF ( Mag(cpl(i)%truer - cpl(j)%truer) .LE. opts%cp_min_dinstance ) THEN
            cpl(j)%hasProxy = .TRUE.
            avgR = (avgR * weight + cpl(j)%truer )/(weight + 1)
            weight = weight + 1
            dupCount = dupCount + 1
            ! The two CP should be the same type!
            IF (cpl(i)%negCount /= cpl(j)%negCount) THEN
              IF (cpl(i)%negCount == 3) CYCLE !Skip validating maxima.
              IF (cpl(j)%negCount == 3) CYCLE !ascension is bugged. 
              ! it some times go to bond CP
              IF (switch) THEN
                PRINT *,'ERROR: TWO TYPES OF CP ARE TOO CLOSE TO EACH OTHER.'
                PRINT *, 'The CPs ',i, j, 'have number of negative eigenvalues: ',cpl(i)%negCount, cpl(j)%negCount
                PRINT *, ''//achar(27)//'[31m STOPPING. NO MORE CP WILL BE OUTPUT.'//achar(27)//'[0m'
                switch = .FALSE.
                isReducible = .FALSE.
                isReduced = .FALSE.
                IF (.NOT. opts%ignore_cp_conflict) THEN
                  EXIT OUTER
                END IF
              END IF
            END IF 
            isReduced = .FALSE.
          END IF
        END DO
        nUCPTNum = nUCPTNum + 1
        ! record the reduced CP
        CALL RecordCPRLight(avgR,chg,rcpl,nUCPTnum, rcpCounts, &
          cpl(i)%ind,.FALSE.)
      END DO OUTER
      CALL ReplaceCPL(cpl,rcpl)
      DO i = 1, SIZE(cpl)
        cpl(i)%isUnique = .TRUE.
      END DO
      ucpCounts = rcpCounts
      ucptnum = nUCPTnum
      IF (LDM) PRINT *, 'The number of duplicate CP found is', dupcount
      DEALLOCATE(rcpl)
    END SUBROUTINE ReduceCP

    SUBROUTINE ReduceCPStatic(cp_static,ucptnum,ucpCounts,isReduced,opts)
      TYPE(cpc), ALLOCATABLE, DIMENSION(:) :: cp_static, reduced_cp_static
      TYPE(options_obj) :: opts
      INTEGER, DIMENSION(4) :: ucpCounts, reduced_ucpCounts
      INTEGER :: ucptnum, reduced_ucptnum, i, j, weight, dupCount
      REAL(q2), DIMENSION(3) :: avg_r
      LOGICAL :: isReduced

      isReduced = .TRUE.
      dupCount = 0
      reduced_ucpCounts = 0
      ! Give the reduced list same length as before, it's ok if a little goes to
      ! waste
      ALLOCATE(reduced_cp_static(ucptnum))
      ! start from a critical point. loop through the entire list, look for
      ! another entry 
      reduced_ucptnum = 0
      DO i = 1, ucptnum
        ! check if this point is already determined as a proxy to some other
        ! point before
        cp_static(i)%weight = 1
        avg_r = cp_static(i)%truer
        IF (cp_static(i)%hasProxy ) THEN 
          CYCLE
        END IF
        DO j = i, ucptnum
          ! Periodic boundary condition will come to haunt me, periodically. 
          IF (j == i) CYCLE
          IF ( Mag(cp_static(i)%truer - cp_static(j)%truer) .LE. opts%cp_min_dinstance ) THEN
            cp_static(j)%hasProxy = .TRUE.
            avg_r = (avg_r * cp_static(i)%weight + cp_static(j)%truer )/(cp_static(i)%weight + 1)
            cp_static(i)%weight = cp_static(i)%weight + 1
            dupCount = dupCount + 1
            ! The two CP should be the same type!
            IF (cp_static(i)%negCount /= cp_static(j)%negCount) THEN
              IF (cp_static(i)%negCount == 3) CYCLE !Skip validating maxima.
              IF (cp_static(j)%negCount == 3) CYCLE !ascension is bugged. 
              ! it some times go to bond CP
              PRINT *,'ERROR: TWO TYPES OF CP ARE TOO CLOSE TO EACH OTHER'
            END IF 
            isReduced = .FALSE.
          END IF
        END DO
        reduced_ucptnum = reduced_ucptnum + 1
        ! record the reduced CP
        CALL RecordCPStatic(i, cp_static, reduced_cp_static, reduced_ucptnum, reduced_ucpCounts)
      END DO
      CALL ReplaceCPL(cp_static,reduced_cp_static)
      DO i = 1, reduced_ucptnum
        cp_static(i)%isUnique = .TRUE.
      END DO
      ucpCounts = reduced_ucpCounts
      ucptnum = reduced_ucptnum
      DEALLOCATE(reduced_cp_static)

    END SUBROUTINE ReduceCPStatic
  
    ! USED IN THIS MODULE
    SUBROUTINE ReplaceCPL(replacee,replacer)
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: replacee,replacer
      INTEGER :: i
      DEALLOCATE (replacee)
      ALLOCATE (replacee(SIZE(replacer)))
      DO i = 1, SIZE(replacer)
        replacee(i) = replacer(i)
      END DO
    END SUBROUTINE ReplaceCPL

    ! USED IN THIS MODULE
    SUBROUTINE ResizeCPL(cpl,newSize)
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) ::cpl,newcpl
      INTEGER :: newSize, i, oldSize
      ALLOCATE (newcpl(newSize))
      oldSize = SIZE(cpl)
      DO i = 1, oldSize
        newcpl(i) = cpl(i)
      END DO
      DEALLOCATE (cpl)
      ALLOCATE (cpl(newSize))
      DO i = 1, oldSize
        cpl(i) = newcpl(i)
      END DO
      DEALLOCATE (newCPL)
    END SUBROUTINE ResizeCPL
 
    ! takes in coordinates, gives out interpolated Hessian using central
    ! difference
    ! USED IN THIS MODULE
    FUNCTION CDHessianR(r,chg)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(8,3,3) :: nnhes
      REAL(q2),DIMENSION(3,3) :: CDHessianR
      REAL(q2),DIMENSION(3) :: r, distance
      INTEGER,DIMENSION(8,3) :: nnind
      INTEGER :: i
      nnind = SimpleNN(r,chg)
      distance = r - nnind(1,:)
      DO i = 1,8
        nnHes(i,:,:) = CDHessian(nnind(i,:),chg)
      END DO
       !row find the nearest neighbors at this new locaiton
       !first update critical point location
       !the next big step is to interpolate the force at predicted critical
       !point.
      CDHessianR = trilinear_interpol_hes(nnHes,distance)
      IF ( GetHessianMag(CDHessianR) < 0.000000001) THEN
        PRINT *, "WARNING, zero Hessian matrix at ", r
      END IF
      RETURN
    END FUNCTION CDHessianR

    ! USED IN THIS MODULE
    FUNCTION CDGradR(r,chg)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(8,3) :: nnGrad
      REAL(q2),DIMENSION(3) :: CDGradR
      REAL(q2),DIMENSION(3) :: r, distance
      INTEGER,DIMENSION(8,3) :: nnind
      INTEGER :: i
      nnind = SimpleNN(r,chg)
      distance = r - nnind(1,:)
      DO i = 1,8
        nnGrad(i,:) = CDGrad(nnind(i,:),chg)
      END DO
       !row find the nearest neighbors at this new locaiton
       !first update critical point location
       !the next big step is to interpolate the force at predicted critical
       !point.
      CDGradR = trilinear_interpol_grad(nnGrad,distance)
      RETURN
    END FUNCTION CDGradR


    ! This subroutine takes in two cartesian coordinates, draw a line in between
    ! with given interval, output charge density, gradient, hessian, tem into
    ! seperate debug files, and terminates program at the end of this function.
    SUBROUTINE DebugLine(iP,fP,iS,chg)
    ! iP is initial Point, fP is final Point, iS is interval Size
    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(8,3,3) :: nnHes
    REAL(q2),DIMENSION(8,3) :: nnGrad
    REAL(q2),DIMENSION(3,3) :: hes
    REAL(q2),DIMENSION(3) :: p, iP, fP, sZ ! sZ is the acutal step size
    REAL(q2),DIMENSION(3) :: grad, distance
    ! sZ should be slightly different from iS due to rounding 
    REAL(q2) :: iS, rho
    INTEGER,DIMENSION(8,3) :: nnind
    INTEGER :: sN ! step number
    INTEGER :: i, j
    OPEN(50,FILE='debug_line_charge.dat',STATUS='REPLACE',ACTION='WRITE')
    OPEN(51,FILE='debug_line_gradient_cart.dat',STATUS='REPLACE',ACTION='WRITE')
    OPEN(52,FILE='debug_line_hessian_cart.dat',STATUS='REPLACE',ACTION='WRITE')
    OPEN(53,FILE='debug_line_rhograd_cart.dat',STATUS='REPLACE',ACTION='WRITE')
    sN = CEILING( Mag( fP - iP ) / iS )
    sZ = ( fP - iP ) / sN
    DO i = 0, sN
      p = iP + i * sZ ! This is in cartesian
      PRINT *, 'Position in cartesian is'
      PRINT *, p
      p = MATMUL(chg%car2lat,p) ! Now it's in lattice
      PRINT *, 'Position in lattie is'
      PRINT *, p
      CALL pbc_r_lat(p,chg%npts)
      grad = rho_grad(chg,p,rho)
      rho = rho
      WRITE (53,*) grad
      PRINT *, "cartesian grad from rho_grad is "
      PRINT *, grad 
      nnind = SimpleNN(p,chg)
      distance = p - nnind(1,:)
      DO j = 1,8
        nngrad(j,:) = CDGrad(nnind(j,:),chg)
        nnhes(j,:,:) = CDHessian(nnind(j,:),chg)
      END DO
      grad = trilinear_interpol_grad(nnGrad,distance) ! val r interpol
      PRINT *, 'cartesian grad from this mod is'
      PRINT *, grad
      hes = trilinear_interpol_hes(nnHes,distance)
      WRITE (50,*) rho
      PRINT *, 'rho is ', rho
      WRITE (51,*) grad
      grad = MATMUL(grad,chg%lat2car)
      PRINT *, 'lattice grad from this mod is '
      PRINT *, grad
      WRITE (52,*) hes(1,:)
      WRITE (52,*) hes(2,:)
      WRITE (52,*) hes(3,:)
      PRINT *, 'hes is'
      PRINT *, hes(1,:)
      PRINT *, hes(2,:)
      PRINT *, hes(3,:)
    END DO
    CLOSE(50)
    CLOSE(51)
    CLOSE(52)
    CLOSE(53)
    END SUBROUTINE

    ! This subroutines takes lattice ranges and output on lattice rho, gradient
    ! and hessian.
    SUBROUTINE DebugRGHLat(xmin,xmax,ymin,ymax,zmin,zmax,chg)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3,3) :: hes
      REAL(q2),DIMENSION(3) :: grad
      INTEGER,DIMENSION(3) :: lat
      INTEGER :: xmin,xmax,ymin,ymax,zmin,zmax,i,j,k
      DO i = xmin,xmax
        DO j = ymin,ymax
          DO k = zmin,zmax
            lat=(/i,j,k/)
            PRINT *, 'lattice point is '
            PRINT *, lat
            PRINT *, 'cartesian point is'
            PRINT *, MATMUL(chg%lat2car,lat)
            PRINT *, 'gradient in cartesian units is'
            grad = CDGrad(lat,chg)
            hes = CDHessian(lat,chg)
            PRINT *, grad
            PRINT *, 'gradient in cartesian is'
            PRINT *, MATMUL(chg%lat2car,grad)
            PRINT *, 'hessian in cartesian is'
            PRINT *, hes(1,:)
            PRINT *, hes(2,:)
            PRINT *, hes(3,:)
            PRINT *, 'hessian in lattice units is'
            hes = MATMUL(TRANSPOSE(chg%lat2car),MATMUL(hes,chg%lat2car))
            PRINT *, hes(1,:)
            PRINT *, hes(2,:)
            PRINT *, hes(3,:)
          END DO
        END DO
      END DO
    END SUBROUTINE
  
    ! USED IN THIS MODULE
    SUBROUTINE DiagonalOnlyHes(hes)
      REAL(q2), DIMENSION(3,3) :: hes
      hes(1,2) = 0
      hes(1,3) = 0
      hes(2,3) = 0
      hes(2,1) = 0
      hes(3,2) = 0
      hes(3,1) = 0
    END SUBROUTINE

    ! This subroutine finds out which point leads to finding of a critical point


    ! This function calculates TEM for a grid point
    ! USED IN THIS MODULE
    FUNCTION CalcTEMGrid(p,chg,grad,hessianMatrix)
      TYPE(charge_obj) :: chg
      INTEGER,DIMENSION(3) :: p
      REAL(q2),DIMENSION(3,3) :: hessianMatrix
      REAL(q2),DIMENSION(3) :: grad,CalcTEMGrid
        grad = CDGrad(p,chg)
        hessianMatrix = CDHessian(p,chg)
        CalcTEMGrid = - MATMUL(INVERSE(hessianMatrix),grad)
        CalcTEMGrid = MATMUL(chg%car2lat,CalcTEMGrid)
        RETURN
    END FUNCTION CalcTEMGrid
 
    ! This function will not work if calculating TEM at a grid point.
    ! USED IN THIS MODULE
    FUNCTION CalcTEMLat(trueR,chg,temScale,previousTEM,grad, &
      temNormCap,LDM)
      TYPE(charge_obj) :: chg
      INTEGER,DIMENSION(8,3) :: nnind
      INTEGER :: j
      REAL(q2),DIMENSION(8,3,3) :: nnHes
      REAL(q2),DIMENSION(8,3) :: nnGrad
      REAL(q2),DIMENSION(3,3) :: hessianMatrix
      REAL(q2),DIMENSION(3) :: temScale,previousTEM,CalcTEMLat,grad&
        ,distance,trueR
      REAL(q2) :: temNormCap
      LOGICAL :: LDM
      nnInd = SimpleNN(trueR,chg)
      IF (LDM) THEN
        PRINT *, "Inside CalcTEMLat"
        PRINT *, "At point "
        PRINT *, trueR
        PRINT *, "nnInd are"
        DO j = 1,8
          PRINT *, nnInd(j,:)
        END DO
      END IF
      distance = truer - nnind(1,:)
      DO j = 1,8
        nngrad(j,:) = CDGrad(nnind(j,:),chg)
        hessianMatrix = CDHessian(nnind(j,:),chg)
        nnHes(j,:,:) = hessianMatrix
      END DO
      grad = trilinear_interpol_grad(nnGrad,distance) ! val r interpol
      hessianMatrix = trilinear_interpol_hes(nnHes,distance)
      IF (LDM) THEN
        PRINT *, "CalcTEMLat Calculated grad is"
        PRINT *, grad
        PRINT *, "CalcTEMLat HessianMatrix is"
        PRINT *, hessianMatrix(1,:)
        PRINT *, hessianMatrix(2,:)
        PRINT *, hessianMatrix(3,:)
        !PRINT *, "All components in nnHes are"
        !DO j=1,8
        !  PRINT *, j
        !  PRINT *, nnHes(j,1,:)
        !  PRINT *, nnHes(j,2,:)
        !  PRINT *, nnHes(j,3,:)
        !END DO
        !PRINT *, "Distance is"
        !PRINT *, distance
      END IF
      
      
      CalcTEMLat = - MATMUL(inverse(hessianMatrix),grad)
      IF (LDM) THEN
        PRINT *, "CalcTEMLat Preconversion TEM is"
        PRINT *, CalcTEMLat
        PRINT *, - MATMUL(inverse(hessianMatrix),grad)
        PRINT *, "Grad used here is"
        PRINT *, grad
        PRINT *, "Hessian used here is"
        PRINT *, hessianMatrix
        PRINT *, "for TemMods, temScale, temNormCap are"
        PRINT *, temScale
        PRINT *, temNormCap
      END IF
      CalcTEMLat  = MATMUL(chg%car2lat,CalcTEMLat)
      CalcTEMLat  = TemMods(CalcTEMLat,temScale,temNormCap)
      temScale = scaleinspector(CalcTEMLat, previousTEM, temScale)
      IF (LDM) THEN
        PRINT *, "CalcTEMLat after car2lat is "
        PRINT *, MATMUL(chg%car2lat,- MATMUL(inverse(hessianMatrix),grad))
        PRINT *, "CalcTEMLat Postconversion TEM is"
        PRINT *, CalcTEMLat
        PRINT *, "Leaving CalcTEMLat"
      END IF
      RETURN
    END FUNCTION CalcTEMLat

    ! USED IN THIS MODULE
    SUBROUTINE GetDebugFlags(opts,LDM,LDM_DetectCircling,&
      LDM_ReduceCP,LDM_DensityDescend,LDM_RecordCPRLight,&
      LDM_NRTFGP,LDM_CalcTEMLat,LDM_RecordCPR,LDM_GradMagGrad,LDM_RingAscend,LDM_Trajectories)
      TYPE(options_obj) :: opts
      CHARACTER(128) :: debugFlags
      INTEGER :: ios
      LOGICAL :: HCF,LDM ! has config file, local debug mode
      LOGICAL :: LDM_RecordCPRLight, LDM_NRTFGP
      LOGICAL :: LDM_DetectCircling, LDM_ReduceCP, LDM_DensityDescend
      LOGICAL :: LDM_CalcTEMLat,LDM_RecordCPR, LDM_GradMagGrad,LDM_RingAscend,LDM_Trajectories
      LDM = .FALSE.
      LDM_DetectCircling = .FALSE.
      LDM_ReduceCP = .FALSE.
      LDM_DensityDescend = .FALSE.
      LDM_RecordCPRLight = .FALSE.
      LDM_GradMagGrad = .FALSE.
      LDM_RingAscend = .FALSE.
      LDM_Trajectories = .FALSE.
      LDM_RecordCPR = .FALSE.
      LDM_CalcTEMLat = .FALSE.
      LDM_NRTFGP = .FALSE.
      
      INQUIRE(FILE="debugConfig",EXIST=HCF)
      IF (HCF) THEN
        OPEN(60,FILE="debugConfig",STATUS='old',ACTION='read',BLANK='null',PAD='yes')
        DO
          READ(60, '(a)',IOSTAT=ios) debugFlags
          IF (debugFlags == "critpoint_find") THEN
            LDM = .TRUE.
            PRINT *, "De Bugger: Debugging SUBROUTINE critpoint_find"
          END IF
          IF (debugFlags == "DetectCircling") THEN
            LDM_DetectCircling = .TRUE.
            PRINT *, "De Bugger: Debugging SUBROUTINE DetectCircling"
          END IF
          IF (debugFlags == "ReduceCP") THEN
            LDM_ReduceCP = .TRUE.
            PRINT *, "De Bugger: Debugging SUBROUTINE ReduceCP"
          END IF
          IF (debugFlags == "DensityDescend") THEN
            LDM_DensityDescend = .TRUE.
            PRINT *, "De Bugger: Debugging SUBROUTINE DensityDescend"
          END IF
          IF (debugFlags == "RecordCPRLight") THEN
            LDM_RecordCPRLight = .TRUE.
            PRINT *, "De Bugger: Debugging SUBROUTINE RecordCPRLight"
          END IF
          IF (debugFlags == "NRTFGP") THEN
            LDM_NRTFGP = .TRUE.
            PRINT *, "De Bugger: Debugging SUBROUTINE NRTFGP"
          END IF
          IF (debugFlags == "CalcTEMLat") THEN
            LDM_CalcTEMLat = .TRUE.
            PRINT *, "De Bugger: Debugging SUBROUTINE CalcTEMLat"
          END IF
          IF (debugFlags == "RecordCPR") THEN
            LDM_RecordCPR = .TRUE.
            PRINT *, "De Bugger: Debugging SUBROUTINE RecordCPR"
          END IF
          IF (debugFlags == "GradMagGrad") THEN
            LDM_GradMagGrad = .TRUE.
            PRINT *, "De Bugger: Debugging SUBROUTINE GradientDescend"
          END IF
          IF (debugFlags == "RingAscend") THEN
            LDM_RingAscend = .TRUE.
            PRINT *, "De Bugger: Debugging SUBROUTINE RingAscension"
          END IF
          IF (debugFlags == "Trajectories") THEN
            LDM_Trajectories = .TRUE.
            PRINT *, "De Bugger: Debugging Trajectories (FUNCTIONs CPRosterAnalysis, OutputCPRoster)"
          END IF
          IF (ios/=0) EXIT
        END DO
        CLOSE(60)
      ELSE
        PRINT *, "ERROR : in debug mode but debugConfig file was not found!"
      END IF
    END SUBROUTINE GetDebugFlags
 
    ! Finishes a Newton Rhapson trajectory from point ind
    ! trueR is the output converged point
    ! isUnique shows if the trajectory converged
    ! USED IN THIS MODULE
    SUBROUTINE NRTFGP(bdr,chg,opts,trueR,&
      isUnique,r,ind,stepMax)
      ! is r used?
      TYPE(bader_obj) :: bdr
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      INTEGER,DIMENSION(8,3) :: nnInd
      INTEGER,DIMENSION(3) :: tempR,ind
      INTEGER :: stepCount,AverageCount,j,stepMax
      REAL(q2),DIMENSION(10,3) :: temList,rList
      REAL(q2),DIMENSION(8,3) :: nnGrad
      REAL(q2),DIMENSION(3) :: grad,averageR,trueR,nexttem,previoustem,&
        prevGrad,distance,temScale,temCap,r,indR
      REAL(q2) :: temNormCap
      LOGICAL :: isUnique
      LOGICAL :: LDM,LDM_detectCircling
      LDM = .FALSE.
      LDM_detectCircling = .FALSE.
      isUnique = .FALSE.
      temcap = (/1.,1.,1./)
      temScale = (/1.,1.,1./)
      temNormCap = 1.
      stepCount = 1
      averageCount = 0
      averageR = (/-1.,-1.,-1./)
      indR(1) = REAL(ind(1),q2)
      indR(2) = REAL(ind(2),q2)
      indR(3) = REAL(ind(3),q2)
      nnInd = SimpleNN(indR,chg)
      ! get gradients " of all nearest neighbors
      DO j = 1, 8
        nngrad(j,:) = CDGrad(nnInd,chg)
      END DO
      previoustem = CalcTEMLat(trueR,chg,temScale,previousTEM,grad,temNormCap,&
        LDM)
      ! Now start newton method iterations
      ! First step
      trueR =  ind + previoustem
      CALL pbc_r_lat(trueR,chg%npts)
      ! All the rest of the steps
      DO stepcount = 1,stepMax
        IF (LDM) PRINT *, "This is step ",stepCount
        IF (stepcount >= 1) THEN
          prevgrad = grad
        END IF
        CALL pbc_r_lat(truer,chg%npts)
        nnind = SimpleNN(truer,chg)
        distance = truer - nnind(1,:)
        nexttem = CalcTEMLat(trueR,chg,temScale,previousTEM,grad,temNormCap,&
          LDM)
        !grad = R2GradInterpol(nnind,truer,chg,nnLayers)
        IF (ABS(grad(1)) <= 0.1*opts%par_gradfloor .AND. &
            ABS(grad(2)) <= 0.1*opts%par_gradfloor .AND. &
            ABS(grad(3)) <= 0.1*opts%par_gradfloor) THEN
          isunique = .TRUE.
          EXIT
        END IF
        CALL DetectCircling(stepCount,rList,temList,trueR,nextTem,averageR, &
          LDM_DetectCircling,ind)
        IF (ALL(averageR /= -1.,1)) THEN
          !cpcl(i)%isUnique = .TRUE.
          trueR = averageR
          ! the following code is temporary. it disables averaging.
          ! upon seeing averaging, this trajectory is marked unusable.
          isUnique = .FALSE.
          EXIT
        END IF
        previoustem = nexttem
        tempr(1) = NINT(truer(1))
        tempr(2) = NINT(truer(2))
        tempr(3) = NINT(truer(3))
        CALL pbc(tempr,chg%npts)
        IF (bdr%volnum(tempr(1), &
            tempr(2),tempr(3)) == bdr%bnum + 1) THEN
          ! We are heading into the vacuum space, cosmonaughts! 
          isunique = .FALSE.
          EXIT
        END IF
        IF ( ABS(nexttem(1)) .LE. 0.1*opts%par_newtonr .AND. &
             ABS(nexttem(2)) .LE. 0.1*opts%par_newtonr .AND. &
             ABS(nexttem(3)) .LE. 0.1*opts%par_newtonr ) THEN
          isUnique = .TRUE.
          EXIT
        END IF
        truer = truer + nexttem

      END DO
      truer = truer + nexttem ! this keeps track the total movement
      CALL pbc_r_lat(truer,chg%npts)
    END SUBROUTINE NRTFGP 

    ! USED IN THIS MODULE
    ! Follow the charge density on grid points down to local minimums
    FUNCTION DensityDescend(chg,bdr,opts,p)
      TYPE(bader_obj) :: bdr
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      INTEGER, DIMENSION(3) :: p,pn
      INTEGER :: n1,n2,n3
      REAL(q2),DIMENSION(3,3) :: hessianMatrix 
      REAL(q2),DIMENSION(3) :: r,grad,DensityDescend
      LOGICAL :: isUnique,minimized
      minimized = .FALSE.
      DO WHILE (.NOT. minimized) 
        outer: DO n1 = -1 , 1
          DO n2 = -1 , 1
            DO n3 = -1 , 1
              pn = p + (/n1,n2,n3/)
              IF ( rho_val(chg,pn(1),pn(2),pn(3)) < rho_val(chg,p(1),p(2),p(3)) ) THEN
                minimized = .FALSE.
                p = pn
                CALL pbc(p,chg%npts)
                EXIT outer
              ELSE 
                minimized = .TRUE.      
              END IF
            END DO
          END DO
        END DO outer
      END DO
      ! Need to obtain initial grad at the grid point and tem
      r = CalcTEMGrid(p,chg,grad,hessianMatrix)
      CALL NRTFGP(bdr,chg,opts,DensityDescend,&
        isUnique,r,p,1000)
      ! note: what if NRTFGP gets carried away? 
      RETURN
    END FUNCTION DensityDescend


    SUBROUTINE DensityDescendAndRecord(chg,bdr,opts,p,cpl,UCPTnum, ucpCounts, &
      noRecording )
      TYPE(bader_obj) :: bdr
      TYPE(charge_obj) :: chg
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      TYPE(options_obj) :: opts
      INTEGER, DIMENSION(4) :: ucpCounts
      INTEGER,DIMENSION(3) :: p,pn,trueInd
      INTEGER :: n1,n2,n3
      INTEGER :: UCPTnum, stepCount
      INTEGER :: uCCbefore,uCCafter !cage count before and after
      REAL(q2),DIMENSION(3,3) :: hessianMatrix 
      REAL(q2),DIMENSION(3) :: trueR,r,grad 
      REAL(q2) :: current_rho, test_rho
      LOGICAL :: isUnique, minimized, noRecording
      minimized = .FALSE.
      stepCount = 0
      current_rho = chg%rho(p(1),p(2),p(3))
      DO WHILE (.NOT. minimized) 
        OUTER: DO n1 = -1 , 1
          DO n2 = -1 , 1
            DO n3 = -1 , 1
              stepCount = stepCount + 1
              pn = p + (/n1,n2,n3/)
              CALL pbc(pn,chg%npts)
              test_rho = chg%rho(pn(1),pn(2),pn(3))
              IF (test_rho < current_rho) THEN
              !IF ( rho_val(chg,pn(1),pn(2),pn(3)) < rho_val(chg,p(1),p(2),p(3)) ) THEN
                current_rho = test_rho
                p = pn
                EXIT OUTER
              ELSE
                minimized = .TRUE.
              END IF
            END DO
          END DO
        END DO OUTER
      END DO
      ! Need to obtain initial grad at the grid point and tem
      !r = CalcTEMGrid(p,chg,grad,hessianMatrix)
      !CALL NRTFGP(bdr,chg,opts,trueR,&
      !  isUnique,r,p,&
      !  1000)
      !! note: what if NRTFGP gets carried away? 
      !IF (ABS(trueR(1) - p(1)) > 1 .OR. ABS(trueR(2) - p(2)) > 1 .OR. &
      !    ABS(trueR(3) - p(3)) > 1) THEN
      !  trueR = p
      !  !PRINT *, "Please report bug that DensityDescend wondered off. The distance is "
      !  !PRINT *, SQRT((trueR(1) - p(1))**2 + (trueR(2) - p(2))**2 + (trueR(3) - p(3))**2)
      !END IF

      !trueInd(:) = NINT(trueR(:))
      trueR = p
      IF (noRecording) THEN
        ! trueR is what we are after
      ELSEIF (bdr%volnum(p(1),p(2),p(3) ) == bdr%bnum + 1) THEN
        ! Density descend went into the vacuum. DO NOTHING.
      ELSE
        UCPTnum = UCPTnum + 1
        uCCbefore = ucpCounts(4)
        CALL RecordDensityDescendCage(p,cpl,ucptnum,ucpCounts,chg)
        !CALL RecordCPRLight(trueR,chg,cpl,UCPTnum, ucpCounts, &
        !  trueInd, .FALSE.)
        uCCafter = ucpCounts(4)
      END IF
      IF (uCCbefore == uCCafter .AND. .FALSE. ) THEN
        PRINT *, "ERROR :: Density descend found a cage point that did not &
          produce three positive eigenvalues!"
        PRINT *, "The found cage point is at ", trueR
      END IF
      ! To catch cages that does not produce three positive eigenvalues:
    END SUBROUTINE DensityDescendAndRecord

    ! USED IN THIS MODULE
    SUBROUTINE PrintNeighborCharges(p,chg)
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3) :: p
    INTEGER :: i,j,k
    PRINT *, "Printing charges of all neighbors at "
    PRINT *, p, rho_val(chg,p(1),p(2),p(3))
    DO i = 1, 3
      DO j = 1,3
        DO k = 1,3
          IF ( i == 0 .AND. j == 0 .AND. k == 0) THEN
            CYCLE
          END IF
          PRINT *, rho_val(chg,p(1)+i,p(1)+j,p(1)+k)
        END DO
      END DO
    END DO
    END SUBROUTINE PrintNeighborCharges



    FUNCTION trace(mat3x3)
      REAL(q2) :: trace
      REAL(q2),DIMENSION(3,3) :: mat3x3
      trace = mat3x3(1,1) + mat3x3(2,2) + mat3x3(3,3)
    END FUNCTION 

    ! Find the eigenvalues by finding roots to the characteristic polynomial
    SUBROUTINE EigvalCharPoly(hessianMatrix,eigvals,eigvecs)
      ! characteristic polynomial of a 3x3 matrix is 
      ! https://mathworld.wolfram.com/CharacteristicPolynomial.html
      ! P3(x) = 1/6 * (trace(A)**3 + 2 trace(A**3) - 3 trace (A) trace(A**2) -
      !         1/2 * (trace(A)**2 - trace(A**2))x + trace(A)x**2 - x**3
      REAL(q2),DIMENSION(3,3) :: hessianMatrix, eigvecs, hessianMatrix2
      REAL(q2),DIMENSION(3) :: eigvals
      REAL(q2) :: a,b,c,d !coefficients for the characteristic polynomial
      REAL(q2) :: b23c ! -b**2 - 3*c
      COMPLEX :: cabbage ! or junk or crap etc.
      COMPLEX :: cabbage1, cabbage2,cabbage3
      COMPLEX :: top1,bot1,top2,bot2
      COMPLEX :: eigval1,eigval2,eigval3
      a = -1
      b = trace(hessianMatrix)
      hessianMatrix2 = MATMUL(hessianMatrix,hessianMatrix)
      c = -(0.5 * (trace(hessianMatrix)**2 - &
        trace(hessianMatrix2)))
      d = (1./6.) * (trace(hessianMatrix)**3 + &
          2. * trace(MATMUL(MATMUL(hessianMatrix,hessianMatrix),hessianMatrix))&
          - 3. * trace(hessianMatrix) * trace(MATMUL(hessianMatrix,hessianMatrix)) )
      ! since a is -1, a is already subsituted in for roots
      b23c = -b**2. - 3.*c
      !PRINT *, "b23c is ", b23c
      cabbage1 = -2.*(b**3) - (9.*b*c) - (27.*d)
      cabbage2 = -b**2. * c**2. - 4.*c**3. +4.*b**3 * d + &
        18.*b*c*d + 27.*d**2
      !PRINT *, "cabbage1 is ", cabbage1
      !PRINT *, "cabbage2 is ", cabbage2
      cabbage3 = 3. * SQRT(3.) * SQRT(cabbage2)
      !PRINT *, "cabbage3 is ", cabbage3
      cabbage = (cabbage1 + cabbage3)**(1./3.)
      !PRINT *, "cabbage is ", cabbage
      top1 = 2.**(1./3.)*b23c
      bot1 = 3*cabbage
      top2 = cabbage
      bot2 = 3.*2.**(1./3.)
      eigval1 = b/3. + top1/bot1 - top2/bot2
      !PRINT *, "eigval1 is ", eigval1
      top1 = (1 + CMPLX(0,SQRT(3.)))*(b23c)
      bot1 = 3.*2.**(2./3.) * cabbage
      top2 = (1-CMPLX(0,SQRT(3.)))*cabbage
      bot2 = 6.*2.**(1./3.)
      eigval2 = b/3. - top1/bot1 + top2/bot2
      !PRINT *, "eigval2 is ", eigval2
      top1 = (1 - CMPLX(0,SQRT(3.)))*b23c
      bot1 = 3.*2.**(2./3.)*cabbage
      top2 = (1 + CMPLX(0,SQRT(3.)))*cabbage
      bot2 = 6.*2.**(1./3.)
      eigval3 = b/3. - top1/bot1 + top2/bot2
      !PRINT *, "eigval3 is ", eigval3
      eigvals(1) = REAL(eigval1)
      eigvals(2) = REAL(eigval2)
      eigvals(3) = REAL(eigval3)
      !PRINT *, "eigvals are"
      !PRINT *, eigvals
    END SUBROUTINE

    ! This subroutine aims at reducing the bips and bumps in a CHGCAR by
    ! averaging it.
    ! USED IN THIS MODULE
    SUBROUTINE SmoothenCHGCAR(chg,avgMode,ions,opts)
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      TYPE(options_obj) :: opts
      REAL(q2),DIMENSION(3) :: rcar
      REAL(q2) :: sumChg, r, rmin, weight
      INTEGER,DIMENSION(3) :: p,pp,rlat
      INTEGER :: n1,n2,n3,i,j,k
      INTEGER :: avgMode, ios
      INTEGER :: counter
      CHARACTER(128) :: line

      ! These code will be in debug mode in the future
      ! first write the headder of CHGCAR from CHGCAR  
      !OPEN(55,FILE='aflow.CHGCAR_sum',STATUS='old',ACTION='read',BLANK='null',PAD='yes')
      OPEN(55,FILE=opts%chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
      OPEN(56,FILE='smoothenedCHGCAR_sum',STATUS='REPLACE',ACTION='WRITE')
      counter = 0
      ! need to differentiate format. 
      DO
        IF(opts%in_opt == opts%in_chgcar5) THEN
          IF (counter == 10 + ions%nions) EXIT
        ELSE
          IF (counter == 9 + ions%nions) EXIT
        END IF
        READ(55, '(A)',IOSTAT=ios) line
        WRITE(56,'(A)') TRIM(line)
        ! CHGCAR has 7 + #atoms + 2 lines before getting into charge densities
        IF (ios/=0) EXIT
        counter = counter + 1
      END DO
      CLOSE(55)
      ! mode:
      ! 1: simple average of 26 neighbors
      ! 2: weighted average. weight is 1/r where r is distance in cartesian.
      !    the cell itself has weight equal to 4/r of the smallest r.
      ! 3: simple average of all neighbors -2 to +2 on all axis
      ! 99: debug mode. simply write out chgcar again
      avgMode = 1
      !newchg = chg
      counter = 0
      DO n1 = 1, chg%npts(1)
        DO n2 = 1, chg%npts(2)
          DO n3 = 1, chg%npts(3)
            p(1) = n1
            p(2) = n2
            p(3) = n3
            IF (avgMode == 2) THEN
              sumChg = 0
              weight = 0
            ELSE
              sumChg = rho_val(chg,n1,n2,n3)
            END IF
            rmin = 99999.
            IF (avgMode == 1 .OR. avgMode == 2) THEN
              DO i = -1,1
                DO j = -1,1
                  DO k = -1,1
                    IF (i == 0 .AND. j == 0 .AND. k == 0) CYCLE
                    pp(1) = n1 + i
                    pp(2) = n2 + j
                    pp(3) = n3 + k
                    CALL pbc(pp,chg%npts)
                    IF (avgMode == 1) THEN
                      sumChg = sumChg + rho_val(chg,pp(1),pp(2),pp(3))
                    ELSE IF (avgMode == 2) THEN
                      rlat(1) = i
                      rlat(2) = j
                      rlat(3) = k
                      rcar = MATMUL(rlat,ions%lattice)
                      r = SQRT(rcar(1)**2 + rcar(2)**2 + rcar(3)**2)
                      IF (rmin < r ) rmin = r
                      sumChg = sumChg + 1/r * rho_val(chg,pp(1),pp(2),pp(3))
                      weight = weight + 1/r
                    END IF
                  END DO
                END DO
              END DO
            ELSE IF (avgMode == 3) THEN
              DO i = -2,2
                DO j = -2,2
                  DO k = -2,2
                    IF (i == 0 .AND. j == 0 .AND. k == 0) CYCLE
                    pp(1) = n1 + i
                    pp(2) = n2 + j
                    pp(3) = n3 + k
                    CALL pbc(pp,chg%npts)
                    sumChg = sumChg + rho_val(chg,pp(1),pp(2),pp(3))
                  END DO
                END DO
              END DO
            ELSE IF (avgMode == 99)  THEN
              ! Does nothing
            END IF
            IF (avgMode == 1) THEN
              sumChg = sumChg/27.
            ELSE IF (avgMode == 2) THEN
              sumChg = sumChg + 4/rmin 
              weight = weight + 4/rmin
              sumChg = sumChg / weight
            ELSE IF (avgMode == 3) THEN
              sumChg = sumChg/126.
            END IF
            IF (counter == 4) THEN
              WRITE(56,'(A)',advance='no') '  '
              WRITE(56,'(ES18.11E2)') sumChg 
              counter = 0
            ELSE
              WRITE(56,'(A)',advance='no') '  '
              WRITE(56,'(ES18.11E2)',advance='no') sumChg
              counter = counter + 1
            END IF
           
          END DO
        END DO
      END DO
      CLOSE(56)
    END SUBROUTINE SmoothenCHGCAR


    ! This subroutine is used when NR trajectories fail to converge.
    ! It calculates the modulus of the charge density gradient
    ! And calculate the gradient of this modulus
    ! And descend the gradient of the modulus.
    ! It takes in only the starting point index, 
    ! and returns the converged point in lattice coordiantes.
    ! USED IN THIS MODULE
    SUBROUTINE GradientDescend(bdr, chg, opts, rn, iniI,isUnique,stepMax)
      INTEGER,DIMENSION(8,3) :: nnInd
      INTEGER,DIMENSION(3) :: iniI, tempint
      REAL(q2),DIMENSION(3) :: tempr
      REAL(q2),DIMENSION(3):: distance
      REAL(q2),DIMENSION(3) :: gMCDG, oldgMCDG 
      ! gradient of modulus of charge density gradient
      REAL(q2),DIMENSION(3) :: maxstepsize
      REAL(q2) :: maxstepsize_mag
      !magMode uses the overall magnitude of the gradient vector for criteria. When off, it checks separately along each axis.
      LOGICAL :: magMode 
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      TYPE(bader_obj) :: bdr
      LOGICAL :: isUnique
      INTEGER ::j, stepcount, loopcount, stepMax
      REAL(q2), DIMENSION(3)::rn 
      REAL(q2), DIMENSION(8,3) :: nngrad
      magMode = opts%GD_magMode
      !maximum step size
      maxstepsize =(/0.5,0.5,0.5/)
      maxstepsize_mag = 0.5
      gMCDG = CDGMCDG(iniI,chg)
      rn(1) = REAL(iniI(1),q2)
      rn(2) = REAL(iniI(2),q2)
      rn(3) = REAL(iniI(3),q2)
      oldgMCDG(1) = 0
      oldgMCDG(2) = 0
      oldgMCDG(3) = 0
      stepcount=0
      loopcount=0
      !The main gradient descent iteration loop
      DO WHILE (loopcount<stepMax)
         loopcount = loopcount + 1
         !checks if the gradient is close enough to zero which implies it is a critical point
         !par_GDgradfloor
         IF (Mag(gMCDG)<= 0.00001 )  THEN
            ! we are at a critical point !
            isUnique= .TRUE.
            EXIT
         END IF
         !checks if maxstepsize is too small
         IF (magMode) THEN
            IF (maxstepsize_mag <= 0.001) THEN
              isUnique = .TRUE.
              EXIT
            END IF
         ELSE
            IF (Mag(maxstepsize) <= 0.001) THEN
              !critical point found by step size becoming too small
              isUnique = .TRUE.
              EXIT
            END IF
         END IF
         !tempr is the gradient vector with minimum step size condition included.
         tempr = gMCDG
         !adjusts the magnitude of each tempr coordinate to be the minimum of its normed magnitude and the max stepsize value.
         !ensures while the magnitudes are constrained, the directions are still correct.
         IF (magMode) THEN
           tempr = MIN(Mag(tempr),maxstepsize_mag) * tempr/Mag(tempr)
         ELSE
           tempr(1) =  MIN(ABS(tempr(1)), maxstepsize(1))* tempr(1)/ABS(tempr(1))
           tempr(2) =  MIN(ABS(tempr(2)), maxstepsize(2))* tempr(2)/ABS(tempr(2))
           tempr(3) =  MIN(ABS(tempr(3)), maxstepsize(3))* tempr(3)/ABS(tempr(3))
         END IF
         !gradient position update
         IF (tempr(1) /= tempr(1) ) tempr(1) = 0
         IF (tempr(2) /= tempr(2) ) tempr(2) = 0
         IF (tempr(3) /= tempr(3) ) tempr(3) = 0
         rn = rn - tempr
         CALL pbc_r_lat(rn,chg%npts)
         nnInd = SimpleNN(rn,chg)
         DO j = 1,8
            CALL pbc(nnind(j,:),chg%npts)
            nngrad(j,:) = CDGMCDG(nnInd(j,:),chg)
         END DO
         distance = rn - nnind(1,:)
          
         tempint(1) = NINT(rn(1))
         tempint(2) = NINT(rn(2))
         tempint(3) = NINT(rn(3))
          
         CALL pbc(tempint,chg%npts) 

         IF (bdr%volnum(tempint(1),tempint(2),tempint(3)) == bdr%bnum + 1) THEN
            isUnique = .FALSE.
            EXIT
         END IF


         oldgMCDG(1) = tempr(1)
         oldgMCDG(2) = tempr(2)
         oldgMCDG(3) = tempr(3)

         gMCDG = trilinear_interpol_grad(nngrad,distance)
        
        !method of looking at each coordinate separately when restricting stepsize
        !If any of the step directions are reversed from the previous step, decrease the overall step size.
        ! IF (gMCDG(1) * oldgMCDG(1) <= 0 .OR.  gMCDG(2) * oldgMCDG(2) <= 0 .OR.  gMCDG(3) * oldgMCDG(3) <= 0 ) THEN
        !     stepsize = stepsize*0.5
        ! END IF
        IF (magMode) THEN
           IF(SUM(gMCDG*oldgMCDG) < 0) THEN
              maxstepsize_mag = 0.5* maxstepsize_mag
           END IF

        ELSE
        
          IF ( gMCDG(1) * oldgMCDG(1) <= 0 ) THEN
            maxstepsize(1) = maxstepsize(1)*0.5
          END IF
          IF (gMCDG(2) * oldgMCDG(2) <= 0 ) THEN
            maxstepsize(2) = maxstepsize(2)*0.5
          END IF
          IF (gMCDG(3) * oldgMCDG(3) <= 0) THEN
            maxstepsize(3) = maxstepsize(3)*0.5
          END IF
      
        END IF
        !If the Dot Product is negative (new gradient in reverse direction), the maximum stepsize is decreased 
        ! IF (SUM(gMCDG*oldgMCDG) .LT. 0) THEN
        !    maxstepsize = 0.5 *maxstepsize
        ! END IF
     END DO     

    END SUBROUTINE GradientDescend

    ! Central Difference Gradient of Modulus of Charge Density Gradient
    ! Operates on grid points
    ! USED IN THIS MODULE
    FUNCTION CDGMCDG(p,chg)
      TYPE(charge_obj) :: chg
      INTEGER,DIMENSION(3) :: p
      INTEGER, DIMENSION(3) :: pzm,pzp,pxm,pxp,pym,pyp
      REAL(q2) :: gradz,gradx,grady
      REAL(q2), DIMENSION(3) :: CDGMCDG
      REAL(q2):: mCDGpzm,mCDGpzp,mCDGpxm,&
        mCDGpxp,mCDGpym,mCDGpyp
      ! First get get the Magnitude of all neary by charge density gradient
      pzm = p + (/0,0,-1/)
      pzp = p + (/0,0,1/)
      pxm = p + (/-1,0,0/)
      pxp = p + (/1,0,0/)
      pym = p + (/0,-1,0/)
      pyp = p + (/0,1,0/)
      CALL pbc(pxm,chg%npts)
      CALL pbc(pym,chg%npts)
      CALL pbc(pzm,chg%npts)
      CALL pbc(pxp,chg%npts)
      CALL pbc(pyp,chg%npts)
      CALL pbc(pzp,chg%npts)
      mCDGpzm = Mag(MATMUL(chg%lat2car,CDGrad(pzm,chg)))
      mCDGpzp = Mag(MATMUL(chg%lat2car,CDGrad(pzp,chg)))
      mCDGpxm = Mag(MATMUL(chg%lat2car,CDGrad(pxm,chg)))
      mCDGpxp = Mag(MATMUL(chg%lat2car,CDGrad(pxp,chg)))
      mCDGpym = Mag(MATMUL(chg%lat2car,CDGrad(pym,chg)))
      mCDGpyp = Mag(MATMUL(chg%lat2car,CDGrad(pyp,chg)))
      gradz = 0.5*(MCDGpzp-MCDGpzm)
      gradx = 0.5*(MCDGpxp-MCDGpxm)
      grady = 0.5*(MCDGpyp-MCDGpym)
      CDGMCDG(1) = REAL(gradx,q2)
      CDGMCDG(2) = REAL(grady,q2)
      CDGMCDG(3) = REAL(gradz,q2)
      RETURN
    END FUNCTION CDGMCDG

    FUNCTION ReadStaticSize()
      INTEGER :: ReadStaticSize
      OPEN(13,FILE="static_search",STATUS='old',ACTION="read")
      READ(13,*) ReadStaticSize
      CLOSE(13)
      RETURN
    END FUNCTION ReadStaticSize

    SUBROUTINE ReadStatic(cps_read)
      TYPE(static_cp_list), ALLOCATABLE, DIMENSION(:) :: cps_read
      CHARACTER(LEN=1) :: cp_type
      REAL(q2), DIMENSION(3) :: pos_direct
      INTEGER :: natoms, ios, n
      
      OPEN(13,FILE="static_search",IOSTAT=ios,STATUS='old',ACTION="read")
      READ(13,*) natoms
      DO n = 1, natoms
        READ(13,*,IOSTAT=ios) cp_type, pos_direct
        cps_read(n)%cp_type = cp_type
        cps_read(n)%pos_direct(:) = pos_direct(:)
      END DO
      CLOSE(13)
      
    END SUBROUTINE ReadStatic

    ! This subroutine is the same as the other version except for it returns a cpc type.
    SUBROUTINE StaticCheckReadStatic(cp_static,chg)
      TYPE(charge_obj) :: chg
      TYPE(cpc), ALLOCATABLE, DIMENSION(:) :: cp_static
      CHARACTER(LEN=1) :: cp_type
      REAL(q2), DIMENSION(3) :: pos_direct
      INTEGER :: natoms, ios, n, it_num, rot_num

      ! only negcount and truer in cp_static needs to be set

      OPEN(13,FILE="static_search",IOSTAT=ios,STATUS='old',ACTION="read")
      READ(13,*) natoms
      DO n = 1, natoms
        READ(13,*,IOSTAT=ios) cp_type, pos_direct
        IF (cp_type == 'n') THEN
          cp_static(n)%negcount = 3
        ELSE IF (cp_type == 'b') THEN        
          cp_static(n)%negcount = 2
        ELSE IF (cp_type == 'r') THEN        
          cp_static(n)%negcount = 1
        ELSE IF (cp_type == 'c') THEN        
          cp_static(n)%negcount = 0
        END IF
        cp_static(n)%truer(:) = pos_direct(:)
        cp_static(n)%hessianMatrix = CDHessianR(pos_direct(:),chg)
        CALL jacobi_eigenvalue(3,cp_static(n)%hessianMatrix,9999,&
          cp_static(n)%eigvecs, cp_static(n)%eigvals,&
          it_num, rot_num )
        cp_static(n)%hasProxy = .FALSE.
      END DO
      CLOSE(13)
    END SUBROUTINE StaticCheckReadStatic


    ! Returns the number of symmetry operations in the AFLOW file.
    FUNCTION ReadSymOpSize()
      INTEGER :: ReadSymOpSize

      RETURN
    END FUNCTION ReadSymOpSize

    ! USED IN THIS MODULE
    SUBROUTINE get_voxvol(chg,ions)
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      voxvol = matrix_volume(ions%lattice)/ chg%nrho
    END SUBROUTINE
  
    ! USED IN THIS MODULE
    FUNCTION GetHessianMag(hes)
      REAL(q2),DIMENSION(3,3) :: hes
      REAL(q2) :: GetHessianMag
      INTEGER :: i,j
      GetHessianMag = 0
      DO i = 1, 3
        DO j = 1,3
          GetHessianMag = GetHessianMag + hes(i,j)**2 
        END DO
      END DO
    END FUNCTION GetHessianMag

    SUBROUTINE PrintFlavorText()
      PRINT *, ''//achar(27)//'[31m Finding Critical points'//achar(27)//'[0m'
      PRINT *, ''//achar(27)//'[91m Interrogation of the soul:'//achar(27)//'[0m'
      PRINT *, ''//achar(27)//'[33m Did I turn on vacuum ?'//achar(27)//'[0m'
      PRINT *, ''//achar(27)//'[92m Did I tell if this is a crystall or molecule?'//achar(27)//'[0m'
      PRINT *, ''//achar(27)//'[36m Did I use the CHGCAR_sum ?'//achar(27)//'[0m'
      PRINT *, ''//achar(27)//'[34m Did I use reasonable values for parameters ? '//achar(27)//'[0m'
      PRINT *, ''//achar(27)//'[95m Critical point is like a box of chocolates. &
        You never know what you are gonna get.'//achar(27)//'[0m'

    END SUBROUTINE

    SUBROUTINE DebugGetRhoAround(p,chg,ions)
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      INTEGER,DIMENSION(3) :: p
      PRINT *, "printing rho_val around ", p       
      PRINT *, "x-1 to x+1 are ", rho_val(chg,p(1)-1,p(2),p(3)),rho_val(chg,p(1),p(2),p(3)),rho_val(chg,p(1)+1,p(2),p(3))
      PRINT *, "y-1 to y+1 are ", rho_val(chg,p(1),p(2)-1,p(3)),rho_val(chg,p(1),p(2),p(3)),rho_val(chg,p(1),p(2)+1,p(3))
      PRINT *, "z-1 to z+1 are ", rho_val(chg,p(1),p(2),p(3)-1),rho_val(chg,p(1),p(2),p(3)),rho_val(chg,p(1),p(2),p(3)+1)
      PRINT *, "voxel length on 3 directions are"
      PRINT *, "x: ", SQRT(ions%lattice(1,1)**2 +ions%lattice(1,2)**2 +ions%lattice(1,3)**2) &
        / chg%npts(1)
      PRINT *, "y: ", SQRT(ions%lattice(2,1)**2 +ions%lattice(2,2)**2 +ions%lattice(2,3)**2) &
        / chg%npts(2)
      PRINT *, "z: ", SQRT(ions%lattice(3,1)**2 +ions%lattice(3,2)**2 +ions%lattice(3,3)**2) &
        / chg%npts(3)
    END SUBROUTINE DebugGetRhoAround

    ! This function should print the density, gradient and curvature of a grid point.
    SUBROUTINE DebugGetInfo(grid,chg)
      TYPE(charge_obj) :: chg
      INTEGER, DIMENSION(3) :: grid
      PRINT *, 'grid point is'
      PRINT *, grid
      PRINT *, 'mag(grad:'
      PRINT *, Mag(CDGrad(grid,chg))
    END SUBROUTINE DebugGetInfo

    SUBROUTINE DebugCPTracer(rt,chg,cpl,ucptnum)
      TYPE(charge_obj) :: chg
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      REAL(q2),DIMENSION(3) :: rt ! This is the target cart location.
      REAL(q2),DIMENSION(3) :: rc ! This is current cart location.
      INTEGER :: i, ucptnum
      DO i = 1, ucptnum
        rc = MATMUL(chg%lat2car,cpl(i)%truer)
        IF (ABS(rc(1) - rt(1)) .le. 0.001 .AND. &
            ABS(rc(2) - rt(2)) .le. 0.001 .AND. &
            ABS(rc(3) - rt(3)) .le. 0.001 ) THEN
          PRINT *, "De Bugger: The given CP is found by trajectory starting at"
          PRINT *, cpl(i)%ind
        END IF
      END DO
    END SUBROUTINE DebugCPTracer

    ! Below is a check on all core functions
    ! Functions being checked: eigenvalues and eigenvectors
    SUBROUTINE DebugCoreFunctionsCheck(chg)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3,3) :: hessianMatrix,eigvecs
      REAL(q2),DIMENSION(3) :: eigvals
      INTEGER :: it_num, rot_num 
      !This following test is not working well
      !hessianMatrix(1,1) = 0.3987508399038720
      !hessianMatrix(1,2) = 0.4997246160218287 
      !hessianMatrix(1,3) = 0.7752726241455388
      !hessianMatrix(2,1) = 0.5993296708147727
      !hessianMatrix(2,2) = 0.6510026529797813
      !hessianMatrix(2,3) = 0.7907512182895509
      !hessianMatrix(3,1) = 0.2018206933392253
      !hessianMatrix(3,2) = 0.7068772764652435
      !hessianMatrix(3,3) = 0.05810136673621315
      !---------------------------------------
      hessianMatrix(1,1) = 4570.188002120697
      hessianMatrix(1,2) = 852.1787784831047
      hessianMatrix(1,3) = 835.9831203358196
      hessianMatrix(2,1) = 852.1787784831047
      hessianMatrix(2,2) = 7811.922523438933
      hessianMatrix(2,3) = 846.7252769646717
      hessianMatrix(3,1) = 835.9831203358196
      hessianMatrix(3,2) = 846.7252769646717
      hessianMatrix(3,3) = 4260.374867575962
      CALL DSYEVJ3(hessianMatrix,eigvecs,eigvals)
      PRINT *, "DSYEVJ3 Produces"
      PRINT *, "eigvec 1"
      PRINT  *, eigvecs(1,:)
      PRINT *, "eigvec 2"
      PRINT  *, eigvecs(2,:)
      PRINT *, "eigvec 3"
      PRINT  *, eigvecs(3,:)
      PRINT *, "eigvals are"
      PRINT *, eigvals
      hessianMatrix(1,1) = 4570.188002120697
      hessianMatrix(1,2) = 852.1787784831047
      hessianMatrix(1,3) = 835.9831203358196
      hessianMatrix(2,1) = 852.1787784831047
      hessianMatrix(2,2) = 7811.922523438933
      hessianMatrix(2,3) = 846.7252769646717
      hessianMatrix(3,1) = 835.9831203358196
      hessianMatrix(3,2) = 846.7252769646717
      hessianMatrix(3,3) = 4260.374867575962
      CALL jacobi_eigenvalue(3,hessianMatrix,9999,eigvecs,eigvals,&
        it_num, rot_num )
      PRINT *, "jacobi Produces"
      PRINT *, "eigvec 1"
      PRINT  *, eigvecs(1,:)
      PRINT *, "eigvec 2"
      PRINT  *, eigvecs(2,:)
      PRINT *, "eigvec 3"
      PRINT  *, eigvecs(3,:)
      PRINT *, "eigvals are"
      PRINT *, eigvals
      PRINT *, "EigvalCharPoly produces"
      hessianMatrix(1,1) = 4570.188002120697
      hessianMatrix(1,2) = 852.1787784831047
      hessianMatrix(1,3) = 835.9831203358196
      hessianMatrix(2,1) = 852.1787784831047
      hessianMatrix(2,2) = 7811.922523438933
      hessianMatrix(2,3) = 846.7252769646717
      hessianMatrix(3,1) = 835.9831203358196
      hessianMatrix(3,2) = 846.7252769646717
      hessianMatrix(3,3) = 4260.374867575962
      CALL EigvalCharPoly(hessianMatrix,eigvals,eigvecs)
      IF (eigvecs(1,1) /= -0.5694273336689270 .OR. &
          eigvecs(1,2) /= -0.7140883142791939 .OR. & 
          eigvecs(1,3) /= -0.4072227781946826 .OR. &
          eigvecs(2,1) /= -0.6717960446248690 .OR. &
          eigvecs(2,2) /= -0.1747987305973604 .OR. &
          eigvecs(2,3) /= 0.7198162808716766 .OR. & 
          eigvecs(3,1) /= -0.7198162808716766 .OR. &
          eigvecs(3,2) /= 0.05039746667582087 .OR. &
          eigvecs(3,3) /= 0.5695024751675655) THEN
        PRINT *, "Eigenvalues produced are not consistent with Mathematica!"
        PRINT *, "eigvec 1"
        PRINT  *, eigvecs(1,:)
        PRINT *, "eigvec 2"
        PRINT  *, eigvecs(2,:)
        PRINT *, "eigvec 3"
        PRINT  *, eigvecs(3,:)
        PRINT *, "eigvals are"
        PRINT *, eigvals
      END IF
    END SUBROUTINE DebugCoreFunctionsCheck

  END MODULE


  
