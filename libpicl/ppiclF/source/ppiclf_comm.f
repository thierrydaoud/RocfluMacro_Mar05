!-----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_comm_InitMPI(comm,id,np)
     > bind(C, name="ppiclc_comm_InitMPI")
#else
      SUBROUTINE ppiclf_comm_InitMPI(comm,id,np)
#endif
!
!     This subroutine is called from rocflu/RFLU_InitFlowSolver.F90
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input: 
!
      INTEGER*4 comm
      INTEGER*4 id
      INTEGER*4 np
!
! Code:
!
      ! Ensures a later subroutine init wasn't called out of order
      IF (PPICLF_LINIT .OR. PPICLF_LFILT .OR. PPICLF_OVERLAP)
     >   CALL ppiclf_exittr('InitMPI must be called first$',0.0d0,0)

      ! set ppiclf_processor information
      ppiclf_comm = comm
      ppiclf_nid  = id
      ppiclf_np   = np

      ! GSlib call
      CALL ppiclf_prints('   *Begin InitCrystal$')
         CALL ppiclf_comm_InitCrystal
      CALL ppiclf_prints('    End InitCrystal$')

      ! check to make sure subroutine is called in correct order later
      ! on in the code sequence
      PPICLF_LCOMM = .true.

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_InitCrystal
!
!     This subroutine is called form ppiclf_comm_InitMPI
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input: 
!

!
! Code:
!
      ! GSlib call
      CALL pfgslib_crystal_setup(ppiclf_cr_hndl,ppiclf_comm,ppiclf_np)

      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_CreateBin
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      INTEGER*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >                            nfacegp, nedgegp, ncornergp
      INTEGER*4 ix, iy, iz, iperiodicx, iperiodicy, iperiodicz, 
     >          npt_total, j, i, idum, jdum, kdum, total_bin, 
     >          sum_value, count, targetTotBin, idealBin(3), iBin(3),
     >          iBinTot, temp,nBinMax,nBinMed,nBinMin, m, l, k,
     >          LBMax,LBMin,
     >          ivx, ivy, ivz
      REAL*8 xmin, ymin, zmin, xmax, ymax, zmax, rduml, rdumr, rthresh,
     >       rmiddle, rdiff, binb_length(3),temp1,temp2,
     >       distchk
      INTEGER*4 ppiclf_iglsum
      EXTERNAL ppiclf_iglsum
      REAL*8 ppiclf_glmin,ppiclf_glmax,ppiclf_glsum
      EXTERNAL ppiclf_glmin,ppiclf_glmax,ppiclf_glsum
!

! face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/
     >                 -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                 -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      IF(ppiclf_ndim .GT. 2) THEN
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      END IF

      ix = PPICLF_JX
      iy = PPICLF_JY
      iz = 1
      IF(ppiclf_ndim .EQ. 3)
     >iz = PPICLF_JZ
      
      ivx = PPICLF_JVX
      ivy = PPICLF_JVY
      ivz = PPICLF_JVZ

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)
        
      ! iglsum is integer addition across MPI ranks.
      npt_total = ppiclf_iglsum(ppiclf_npart,1)
      ! compute bin boundaries
      xmin = 1E10
      ymin = 1E10
      zmin = 1E10
      xmax = -1E10
      ymax = -1E10
      zmax = -1E10
      ! Looping through particles on this processor
      ! Find the bin boundaries based on the REAL particles
      DO i=1,ppiclf_npart
         ! Finding min/max REAL particle extremes.
         ! Need to consider filter/neighborwidths
         ! to ensure ppiclf_bins_dx > ppiclf_d2chk(1)
         rduml = ppiclf_y(ix,i) - ppiclf_d2chk(1)
         rdumr = ppiclf_y(ix,i) + ppiclf_d2chk(1)
         IF(rduml .LT. xmin) xmin = rduml
         IF(rdumr .GT. xmax) xmax = rdumr

         rduml = ppiclf_y(iy,i) - ppiclf_d2chk(1)
         rdumr = ppiclf_y(iy,i) + ppiclf_d2chk(1)
         IF(rduml .LT. ymin) ymin = rduml
         IF(rdumr .GT. ymax) ymax = rdumr

         IF(ppiclf_ndim .EQ. 3) THEN
            rduml = ppiclf_y(iz,i) - ppiclf_d2chk(1)
            rdumr = ppiclf_y(iz,i) + ppiclf_d2chk(1)
            IF(rduml .LT. zmin) zmin = rduml
            IF(rdumr .GT. zmax) zmax = rdumr
         END IF
      END DO

      ! Finds global max/mins across MPI ranks
      ppiclf_binb(1) = ppiclf_glmin(xmin,1)
      ppiclf_binb(2) = ppiclf_glmax(xmax,1)
      ppiclf_binb(3) = ppiclf_glmin(ymin,1)
      ppiclf_binb(4) = ppiclf_glmax(ymax,1)
      ppiclf_binb(5) = 0.0d0
      ppiclf_binb(6) = 0.0d0
      IF(ppiclf_ndim .GT. 2) ppiclf_binb(5) = ppiclf_glmin(zmin,1)
      IF(ppiclf_ndim .GT. 2) ppiclf_binb(6) = ppiclf_glmax(zmax,1)




      ! For linear periodicity:
      !     1) bin boundary greater than periodic domain      -> set it equal to periodic domain
      !     2) bin boundary within distchk to periodic domain -> set it equal to periodic domain
      !     3) ghost creation is now taken care of by the periodic planes

      ! when ppiclf_linear_binbx,y,z is .true. -> possibility of linear particle colliding on the wrapped side
      !                                           with another particle
      !                                        -> evauate linear ghost algorithm subroutine

      distchk = ppiclf_d2chk(1)                                    ! max(filter_width, neighbor_width)

      if((ppiclf_binb(1).lt.ppiclf_xdrange(1,1)) .or.            ! min binb .lt. periodic domain
     >  (abs(ppiclf_binb(1)-ppiclf_xdrange(1,1)).le.distchk)) then ! min binb within distchk to min periodic domain
          ppiclf_binb(1) = ppiclf_xdrange(1,1)
        if( iperiodicx .eq. 0) ppiclf_linear_bxmin = .true.
      endif
      if((ppiclf_binb(2) .gt. ppiclf_xdrange(2,1)) .or.          ! max binb .gt. periodic domain
     >  (abs(ppiclf_binb(2)-ppiclf_xdrange(2,1)).le.distchk)) then ! max binb within distchk to max periodic domain
          ppiclf_binb(2) = ppiclf_xdrange(2,1)
        if( iperiodicx .eq. 0) ppiclf_linear_bxmax = .true.
      endif

        if(ppiclf_linear_bxmin.and.ppiclf_linear_bxmax)          ! bool check to evaluate linear periodic ghost 
     >   ppiclf_linear_bx = .true.                                 ! only when this is true


      if((ppiclf_binb(3) .lt. ppiclf_xdrange(1,2)) .or.
     > (abs(ppiclf_binb(3)-ppiclf_xdrange(1,2)).le.distchk)) then
          ppiclf_binb(3) = ppiclf_xdrange(1,2)
        if(iperiodicy .eq. 0) ppiclf_linear_bymin = .true.
      endif

      if((ppiclf_binb(4) .gt. ppiclf_xdrange(2,2)) .or. 
     > (abs(ppiclf_binb(4)-ppiclf_xdrange(2,2)).le.distchk)) then
          ppiclf_binb(4) = ppiclf_xdrange(2,2)
        if(iperiodicy .eq. 0) ppiclf_linear_bymax = .true.
      endif

        if(ppiclf_linear_bymin.and.ppiclf_linear_bymax)          ! bool check to evaluate linear periodic ghost 
     >   ppiclf_linear_by = .true.                                 ! only when this is true

      if(ppiclf_ndim .gt. 2) then                                  ! 3D runs
        if((ppiclf_binb(5) .lt. ppiclf_xdrange(1,3)) .or.
     >   (abs(ppiclf_binb(5)-ppiclf_xdrange(1,3)).le.distchk)) then
            ppiclf_binb(5) = ppiclf_xdrange(1,3)
          if(iperiodicz .eq. 0) ppiclf_linear_bzmin = .true.
        endif
        if((ppiclf_binb(6) .gt. ppiclf_xdrange(2,3)) .or. 
     >   (abs(ppiclf_binb(6)-ppiclf_xdrange(2,3)).le.distchk)) then
            ppiclf_binb(6) = ppiclf_xdrange(2,3)
          if(iperiodicz .eq. 0) ppiclf_linear_bzmax = .true.
        endif
        if(ppiclf_linear_bzmin.and.ppiclf_linear_bzmax)          ! bool check to evaluate linear periodic ghost 
     >   ppiclf_linear_bz = .true.                                 ! only when this is true
      endif ! ppiclf_ndim


!      if (npt_total .gt. 0) then
!      do i=1,ppiclf_ndim
!         if (ppiclf_bins_balance(i) .eq. 1) then
!            rmiddle = 0.0
!            do j=1,ppiclf_npart
!               rmiddle = rmiddle + ppiclf_y(i,j)
!            enddo
!            rmiddle = ppiclf_glsum(rmiddle,1)
!            rmiddle = rmiddle/npt_total
!
!            rdiff =  max(abs(rmiddle-ppiclf_binb(2*(i-1)+1)),
!     >                   abs(ppiclf_binb(2*(i-1)+2)-rmiddle))
!            ppiclf_binb(2*(i-1)+1) = rmiddle - rdiff
!            ppiclf_binb(2*(i-1)+2) = rmiddle + rdiff
!         endif
!      enddo
!      endif

      ! David's Algorithm sets the bin boundaries as big as 
      ! the periodic direction when periodicity is invoked
      ! we comment it out now as this is taken care of by invoking
      ! the periodic plane for both linear and angular periodicity
!      if (ppiclf_xdrange(2,1) .lt. ppiclf_binb(2) .or.
!     >    ppiclf_xdrange(1,1) .gt. ppiclf_binb(1) .or. 
!     >    iperiodicx .eq. 0) then
!         ppiclf_binb(1) = ppiclf_xdrange(1,1)
!         ppiclf_binb(2) = ppiclf_xdrange(2,1)
!      endif
!
!      if (ppiclf_xdrange(2,2) .lt. ppiclf_binb(4) .or.
!     >    ppiclf_xdrange(1,2) .gt. ppiclf_binb(3) .or.
!     >    iperiodicy .eq. 0) then
!         ppiclf_binb(3) = ppiclf_xdrange(1,2)
!         ppiclf_binb(4) = ppiclf_xdrange(2,2)
!      endif
!      
!      if (ppiclf_ndim .gt. 2) then
!      if (ppiclf_xdrange(2,3) .lt. ppiclf_binb(6) .or.
!     >    ppiclf_xdrange(1,3) .gt. ppiclf_binb(5) .or. 
!     >    iperiodicz .eq. 0) then
!         ppiclf_binb(5) = ppiclf_xdrange(1,3)
!         ppiclf_binb(6) = ppiclf_xdrange(2,3)
!      endif ! ndim
!      endif ! xdrange

      ! End subroutine if no particles present      
      IF(npt_total .LT. 1) RETURN
      LBMax = 0
      LBMin = 0
      temp1 = 1.0D-10
      temp2 = 1.0D10
      ! Find ppiclf bin domain lengths
      ! and Max,Med,Min dimensions
      DO l = 1,3
        binb_length(l) = ppiclf_binb(2*l) -
     >                         ppiclf_binb(2*l-1)
        IF(binb_length(l).GT.temp1) THEN
          temp1 = binb_length(l)
          LBMax = l
        END IF
        IF(binb_length(l).LT.temp2) THEN
          temp2 = binb_length(l)
          LBMin = l
        END IF
      END DO

      IF(ppiclf_ndim .LT. 3)
     >   CALL ppiclf_exittr('CreateBins only supports 3D Grids',0.0D0,0)
      
      ! Update with targetTotBin = ActiveBinNum
      targetTotBin = ppiclf_np

      ! Number of bins calculated based on bin surface
      ! area minimization and bin aspect ratio close to 1
      ppiclf_n_bins(1) = INT((targetTotBin**(1.0D0/3.0D0))*
     >                   (binb_length(1)**(2.0D0/3.0D0))/ 
     >                   ((binb_length(2)**(1.0D0/3.0D0))*
     >                   (binb_length(3))**(1.0D0/3.0D0)))
      
      ppiclf_n_bins(2) = INT((targetTotBin**(1.0D0/3.0D0))*
     >                   (binb_length(2)**(2.0D0/3.0D0))/ 
     >                   ((binb_length(1)**(1.0D0/3.0D0))*
     >                   (binb_length(3))**(1.0D0/3.0D0)))
     
      ppiclf_n_bins(3) = INT((targetTotBin**(1.0D0/3.0D0))*
     >                   (binb_length(3)**(2.0D0/3.0D0))/ 
     >                   ((binb_length(2)**(1.0D0/3.0D0))*
     >                   (binb_length(1))**(1.0D0/3.0D0)))


      ! Since INT trucates, make sure n_bins at least 1 
      DO l = 1,3
        IF(ppiclf_n_bins(l) .LT. 1) ppiclf_n_bins(l) = 1
      END DO

      iBinTot = 0

!      print*, "**************************************"
!      print*, "In CreateBin, loc 0"
!      print*, "ppiclf_nid", ppiclf_nid
!      print*, "ppiclf_binb =", ppiclf_binb
!      print*, 'ppiclf_xdrange(1,:) =', ppiclf_xdrange(1,:)
!      print*, 'ppiclf_xdrange(2,:) =', ppiclf_xdrange(2,:)
!      print*, 'ppiclf_bins_dx =', ppiclf_bins_dx
!      print*, 'ppiclf_n_bins =', ppiclf_n_bins
!      print*, 'binb_length(1:3) =', binb_length
!      print*, 'targetTotBin =', targetTotBin
!      print*, 'ppiclf_d2chk(1) =', ppiclf_d2chk(1)
!      print*, "**************************************"

      ! Filterwidth criteria check.  ppiclf_d2chk(2) automatically
      ! set to be at least 2 fluid cell widths in
      ! PICL_TEMP_InitSolver.F90

      DO l = 1,3
          ! Ensure ppiclf_bin_dx(l) > ppiclf_d2chk(1) 
          IF((binb_length(l)/ppiclf_n_bins(l)) .LT. ppiclf_d2chk(1)) 
     >      ppiclf_n_bins(l) = FLOOR(binb_length(l)/ppiclf_d2chk(1))
!     >      ppiclf_n_bins(l) = INT(ppiclf_n_bins(l)/ppiclf_d2chk(1))
          IF(ppiclf_n_bins(l) .LT. 1)  
     >  CALL ppiclf_exittr('ppiclf_d2chk(1) criteria violated.',0.0D0,0)
        idealBin(l) = ppiclf_n_bins(l)
      END DO

!      print*, "**************************************"
!      print*, "In CreateBin, loc 1"
!      print*, "ppiclf_nid", ppiclf_nid
!      print*, "ppiclf_binb =", ppiclf_binb
!      print*, 'ppiclf_xdrange(1,:) =', ppiclf_xdrange(1,:)
!      print*, 'ppiclf_xdrange(2,:) =', ppiclf_xdrange(2,:)
!      print*, 'ppiclf_bins_dx =', ppiclf_bins_dx
!      print*, 'ppiclf_n_bins =', ppiclf_n_bins
!      print*, 'binb_length(1:3) =', binb_length
!      print*, 'targetTotBin =', targetTotBin
!      print*, 'ppiclf_d2chk(1) =', ppiclf_d2chk(1)
!      print*, "**************************************"

      ! Since bin must be an integer, check -1, +0, +1 number of bins for each bin dimension
      ! ideal number of bins will be max value while less than number of total target of bins.
      ! Will not check total bin value (cycle do loop) if
      ! ppiclf_d2chk(1) criteria is violated or ppiclf_n_bins < 1

      total_bin = 0 
      DO ix = 1,3
        iBin(1) = ppiclf_n_bins(1) + (ix-2)
        ppiclf_bins_dx(1) = binb_length(1)/iBin(1)
        IF(ppiclf_bins_dx(1) .LT. ppiclf_d2chk(1) .OR.
     >                           iBin(1) .LT. 1) CYCLE
        DO iy = 1,3
          iBin(2) = ppiclf_n_bins(2) + (iy-2)
          ppiclf_bins_dx(2) = binb_length(2)/iBin(2)
          IF(ppiclf_bins_dx(2) .LT. ppiclf_d2chk(1) .OR.
     >                             iBin(2) .LT. 1) CYCLE
          DO iz = 1,3
            iBin(3) = ppiclf_n_bins(3) + (iz-2)
            ppiclf_bins_dx(3) = binb_length(3)/iBin(3)
            IF(ppiclf_bins_dx(3) .LT. ppiclf_d2chk(1) .OR.
     >                               iBin(3) .LT. 1) CYCLE
            iBinTot = iBin(1)*iBin(2)*iBin(3)
            IF(iBinTot .GT. total_bin .AND.
     >                     iBinTot .LE. targetTotBin) THEN
              total_bin = 1
              DO l = 1,3
                idealBin(l) = iBin(l)
                total_bin = total_bin*idealBin(l)
              END DO
              ! These loops are to make sure the dimension with the longest
              ! ppiclf_binb length gets more bins in the case where two or
              ! more dimensions are within 1 bin division of each other.
              temp = 0
              nBinMax = MAX(idealBin(1),idealBin(2),idealBin(3))
              nBinMin = MIN(idealBin(1),idealBin(2),idealBin(3))
              nBinMed = -99
              DO l = 1,3
                IF(idealBin(l).LT.nBinMax .AND. idealBin(l).GT.nBinMin)
     >             nBinMed = idealBin(l)
              END DO
              IF(nBinMed.EQ. -99) THEN !two number of bins are equal
                DO l = 1,3
                  IF(idealBin(l).EQ.nBinMax) temp = temp + 1
                  IF(idealBin(l).EQ.nBinMin) temp = temp + 10
                END DO
                IF(temp .EQ. 2) THEN
                  nBinMed = nBinMax
                ELSE ! Either two nBinMin or all 3 equal
                  nBinMed = nBinMin
                END IF
              END IF
              DO l = 1,3
                IF(l.EQ.LBMax) THEN
                  idealBin(l)=nBinMax
                ELSE IF(l.EQ.LBMin) THEN
                  idealBin(l)=nBinMin
                ELSE
                  idealBin(l)=nBinMed 
                END IF
              END DO 
            END IF
          END DO !iz
        END DO !iy
      END DO !ix

!      print*, "**************************************"
!      print*, "In CreateBin, loc 2"
!      print*, "ppiclf_nid", ppiclf_nid
!      print*, "ppiclf_binb =", ppiclf_binb
!      print*, 'ppiclf_xdrange(1,:) =', ppiclf_xdrange(1,:)
!      print*, 'ppiclf_xdrange(2,:) =', ppiclf_xdrange(2,:)
!      print*, 'ppiclf_bins_dx =', ppiclf_bins_dx
!      print*, 'ppiclf_n_bins =', ppiclf_n_bins
!      print*, 'binb_length(1:3) =', binb_length
!      print*, 'targetTotBin =', targetTotBin
!      print*, 'ppiclf_d2chk(1) =', ppiclf_d2chk(1)
!      print*, "**************************************"

      ! Set common ppiclf arrays based on above calculation
      DO l = 1,3
        ppiclf_n_bins(l) = idealBin(l)
        ppiclf_bins_dx(l) = binb_length(l)/ppiclf_n_bins(l)
      END DO


      ! Loop to see if we can add one to dimension with largest number of bins
      ! Choose this dimension because it is smallest incremental increase to total bins 
      DO
        IF((total_bin/ppiclf_n_bins(LBMax))*
     >      (ppiclf_n_bins(LBMax)+1) .LT. targetTotBin) THEN
          ! Add a bin and set new bin dx length
          ppiclf_n_bins(LBMax) = ppiclf_n_bins(LBMax)+1
          ppiclf_bins_dx(LBMax) = binb_length(LBMax)/
     >                              ppiclf_n_bins(LBMax)
          IF(ppiclf_bins_dx(LBMax) .LT. ppiclf_d2chk(1)) THEN
            ! If ppiclf_d2chk criteria violated, return to previous bin configuration
            ppiclf_n_bins(LBMax) = ppiclf_n_bins(LBMax)-1
            ppiclf_bins_dx(LBMax) = binb_length(LBMax)/
     >                                ppiclf_n_bins(LBMax)
            EXIT
          END IF
          total_bin = 1
          DO l = 1,3
            total_bin = total_bin*ppiclf_n_bins(l)
          END DO
        ELSE
          EXIT
        END IF
      END DO
      IF(total_bin .GT. ppiclf_np) THEN
        PRINT*, 'ERROR: Num Bins > NumProcessors',total_bin,ppiclf_np
        CALL ppiclf_exittr('Error in Createbins',0.0,0)
      END IF
          

!      print*, "**************************************"
!      print*, "In CreateBin, loc 3"
!      print*, "ppiclf_nid", ppiclf_nid
!      print*, "ppiclf_binb =", ppiclf_binb
!      print*, 'ppiclf_xdrange(1,:) =', ppiclf_xdrange(1,:)
!      print*, 'ppiclf_xdrange(2,:) =', ppiclf_xdrange(2,:)
!      print*, 'ppiclf_bins_dx =', ppiclf_bins_dx
!      print*, 'ppiclf_n_bins =', ppiclf_n_bins
!      print*, 'binb_length(1:3) =', binb_length
!      print*, 'targetTotBin =', targetTotBin
!      print*, 'ppiclf_d2chk(1) =', ppiclf_d2chk(1)
!      print*, "**************************************"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! David's old binning method left below for now
!      finished(1) = 0
!      finished(2) = 0
!      finished(3) = 0
!      total_bin = 1 
!
!      do i=1,ppiclf_ndim
!         finished(i) = 0
!         exit_1_array(i) = ppiclf_bins_set(i)
!         exit_2_array(i) = 0
!         if (ppiclf_bins_set(i) .ne. 1) ppiclf_n_bins(i) = 1
!         ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                        ppiclf_binb(2*(i-1)+1)  ) / 
!     >                       ppiclf_n_bins(i)
!         ! Make sure exit_2 is not violated by user input
!         if (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1)) then
!            do while (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1))
!               ppiclf_n_bins(i) = max(1, ppiclf_n_bins(i)-1)
!               ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                              ppiclf_binb(2*(i-1)+1)  ) / 
!     >                             ppiclf_n_bins(i)
!         WRITE(*,*) "Inf. loop in CreateBin", i, 
!     >              ppiclf_bins_dx(i), ppiclf_d2chk(1)
!         call ppiclf_exittr('Inf. loop in CreateBin$',0.0,0)
!            enddo
!         endif
!         total_bin = total_bin*ppiclf_n_bins(i)
!      enddo
!
!      ! Make sure exit_1 is not violated by user input
!      count = 0
!      do while (total_bin > ppiclf_np)
!          count = count + 1;
!          i = modulo((ppiclf_ndim-1)+count,ppiclf_ndim)+1
!          ppiclf_n_bins(i) = max(ppiclf_n_bins(i)-1,1)
!          ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                         ppiclf_binb(2*(i-1)+1)  ) / 
!     >                        ppiclf_n_bins(i)
!          total_bin = 1
!          do j=1,ppiclf_ndim
!             total_bin = total_bin*ppiclf_n_bins(j)
!          enddo
!          if (total_bin .le. ppiclf_np) exit
!       enddo
!
!       exit_1 = .false.
!       exit_2 = .false.
!
!       do while (.not. exit_1 .and. .not. exit_2)
!          do i=1,ppiclf_ndim
!             if (exit_1_array(i) .eq. 0) then
!                ppiclf_n_bins(i) = ppiclf_n_bins(i) + 1
!                ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                               ppiclf_binb(2*(i-1)+1)  ) / 
!     >                              ppiclf_n_bins(i)
!
!                ! Check conditions
!                ! exit_1
!                total_bin = 1
!                do j=1,ppiclf_ndim
!                   total_bin = total_bin*ppiclf_n_bins(j)
!                enddo
!                if (total_bin .gt. ppiclf_np) then
!                   ! two exit arrays aren't necessary for now, but
!                   ! to make sure exit_2 doesn't slip through, we
!                   ! set both for now
!                   exit_1_array(i) = 1
!                   exit_2_array(i) = 1
!                   ppiclf_n_bins(i) = ppiclf_n_bins(i) - 1
!                   ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                                  ppiclf_binb(2*(i-1)+1)  ) / 
!     >                                  ppiclf_n_bins(i)
!                   exit
!                endif
!                
!                ! exit_2
!                if (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1)) then
!                   ! two exit arrays aren't necessary for now, but
!                   ! to make sure exit_2 doesn't slip through, we
!                   ! set both for now
!                   exit_1_array(i) = 1
!                   exit_2_array(i) = 1
!                   ppiclf_n_bins(i) = ppiclf_n_bins(i) - 1
!                   ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                                  ppiclf_binb(2*(i-1)+1)  ) / 
!     >                                  ppiclf_n_bins(i)
!                   exit
!                endif
!             endif
!          enddo
!
!          ! full exit_1
!          sum_value = 0
!          do i=1,ppiclf_ndim
!             sum_value = sum_value + exit_1_array(i)
!          enddo
!          if (sum_value .eq. ppiclf_ndim) then
!             exit_1 = .true.
!          endif
!
!          ! full exit_2
!          sum_value = 0
!          do i=1,ppiclf_ndim
!             sum_value = sum_value + exit_2_array(i)
!          enddo
!          if (sum_value .eq. ppiclf_ndim) then
!             exit_2 = .true.
!          endif
!       enddo
!      ! Check for too small bins 
!      rthresh = 1E-12
!      total_bin = 1
!      do i=1,ppiclf_ndim
!         total_bin = total_bin*ppiclf_n_bins(i)
!         if (ppiclf_bins_dx(i) .lt. rthresh) ppiclf_bins_dx(i) = 1.0
!      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! -------------------------------------------------------
! SETUP 3D BACKGROUND GRID PARAMETERS FOR GHOST PARTICLES
! -------------------------------------------------------

!     current box coordinates
      IF(ppiclf_nid .LE. total_bin-1) THEN
         idum = modulo(ppiclf_nid,ppiclf_n_bins(1))
         jdum = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
         kdum = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
         IF(ppiclf_ndim .LT. 3) kdum = 0
         ppiclf_binx(1,1) = ppiclf_binb(1) + idum    *ppiclf_bins_dx(1)
         ppiclf_binx(2,1) = ppiclf_binb(1) + (idum+1)*ppiclf_bins_dx(1)
         ppiclf_biny(1,1) = ppiclf_binb(3) + jdum    *ppiclf_bins_dx(2)
         ppiclf_biny(2,1) = ppiclf_binb(3) + (jdum+1)*ppiclf_bins_dx(2)
         ppiclf_binz(1,1) = 0.0d0
         ppiclf_binz(2,1) = 0.0d0
         IF(ppiclf_ndim .GT. 2) THEN
            ppiclf_binz(1,1) = ppiclf_binb(5)+kdum    *ppiclf_bins_dx(3)
            ppiclf_binz(2,1) = ppiclf_binb(5)+(kdum+1)*ppiclf_bins_dx(3)
         END IF
      END IF

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateSubBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 nbin, idum, jdum, kdum, ndumx, ndumy, itmp, jtmp, ktmp,
     >          i, j, k
!

      nbin = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)

c     current box coordinates
      if (ppiclf_nid .le. nbin-1) then
         idum = modulo(ppiclf_nid,ppiclf_n_bins(1))
         jdum = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
         kdum = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
         if (ppiclf_ndim .lt. 3) kdum = 0
         ! interior grid of each bin
         ! +1 for making mesh smaller and +1 since these are vertice counts
         ppiclf_bx = floor(ppiclf_bins_dx(1)/ppiclf_filter) + 1 + 1
         ppiclf_by = floor(ppiclf_bins_dx(2)/ppiclf_filter) + 1 + 1
         ppiclf_bz = 1
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_bz = floor(ppiclf_bins_dx(3)/ppiclf_filter) + 1 + 1

         ppiclf_bx = ppiclf_bx*ppiclf_ngrids
         ppiclf_by = ppiclf_by*ppiclf_ngrids
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_bz = ppiclf_bz*ppiclf_ngrids

         if (ppiclf_bx .gt. PPICLF_BX1)
     >      call ppiclf_exittr('Increase PPICLF_BX1$',0.,ppiclf_bx)
         if (ppiclf_by .gt. PPICLF_BY1)
     >      call ppiclf_exittr('Increase PPICLF_BY1$',0.,ppiclf_by)
         if (ppiclf_bz .gt. PPICLF_BZ1)
     >      call ppiclf_exittr('Increase PPICLF_BZ1$',0.,ppiclf_bz)

         ppiclf_rdx = ppiclf_bins_dx(1)/(ppiclf_bx-1)
         ppiclf_rdy = ppiclf_bins_dx(2)/(ppiclf_by-1)
         ppiclf_rdz = 0
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_rdz = ppiclf_bins_dx(3)/(ppiclf_bz-1)

         ndumx = ppiclf_n_bins(1)*(ppiclf_bx-1) + 1
         ndumy = ppiclf_n_bins(2)*(ppiclf_by-1) + 1
    
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            ppiclf_grid_x(i,j,k) = sngl(ppiclf_binx(1,1) +
     >                                  (i-1)*ppiclf_rdx)
            ppiclf_grid_y(i,j,k) = sngl(ppiclf_biny(1,1) +
     >                                  (j-1)*ppiclf_rdy)
            ppiclf_grid_z(i,j,k) = sngl(ppiclf_binz(1,1) +
     >                                  (k-1)*ppiclf_rdz)

            itmp = idum*(ppiclf_bx-1) + (i-1)
            jtmp = jdum*(ppiclf_by-1) + (j-1)
            ktmp = kdum*(ppiclf_bz-1) + (k-1)
    
            ppiclf_grid_i(i,j,k)  = itmp + ndumx*jtmp + ndumx*ndumy*ktmp

         enddo
         enddo
         enddo

      endif

      return
      end
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MapOverlapMesh
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE 'mpif.h'
!
! Internal:
!
      INTEGER*4 icalld
      SAVE      icalld
      DATA      icalld /0/
      INTEGER*4 nkey(2), i, j, k,l, ie, iee, ii, jj, kk, ndum, nrank,
     >          nl, nii, njj, nrr, ilow, jlow, klow, nxyz, il,
     >          ihigh, jhigh, khigh, ierr
      INTEGER*4 ix, iy, iz, ixLow, ixHigh, iyLow,
     >          iyHigh, izLow, izHigh 
      REAL*8    rxval, ryval, rzval, EleSizei(3), MaxPoint(3),
     >          MinPoint(3), ppiclf_vlmin, ppiclf_vlmax,
     >          centeri(3)
      LOGICAL   partl, ErrorFound
      EXTERNAL  ppiclf_vlmin, ppiclf_vlmax
!
! Code Start:
!

      nxyz = PPICLF_LEZ*PPICLF_LEY*PPICLF_LEX !Num of vertices per cell
      ppiclf_neltb = 0 !counts number of Rocflu elements on this processor
                       !that are within the ppiclf bounds domain
      DO ie=1,ppiclf_nee ! Number of fluid cells in Rocflu Grid domain
        ! Find fluid cell max x,y,z lengths and centroid
        DO l=1,3
          centeri(l)  = 0.0D0
          MaxPoint(l) = -1000000.0D0
          MinPoint(l) =  1000000.0D0 
          EleSizei(l) =  0.0D0
          DO k=1,PPICLF_LEZ
            DO j=1,PPICLF_LEY
              DO i=1,PPICLF_LEX
                centeri(l) = centeri(l) + ppiclf_xm1bs(i,j,k,l,ie)
                IF (ppiclf_xm1bs(i,j,k,l,ie) .GT. MaxPoint(l)) 
     >              MaxPoint(l) = ppiclf_xm1bs(i,j,k,l,ie)
                IF (ppiclf_xm1bs(i,j,k,l,ie) .LT. MinPoint(l)) 
     >              MinPoint(l) = ppiclf_xm1bs(i,j,k,l,ie)
              END DO !i
            END DO !j
          END DO !k
          centeri(l) = centeri(l)/nxyz

          ! 1.2 times the cell length to ensure that one layers of cells
          ! outside of the ppiclf bin are mapped for interpolation.
          ! This could be changed based on the frequency of ppiclf bin
          ! creation and mapping.
          EleSizei(l) = 1.2*(MaxPoint(l) - MinPoint(l))
          IF(EleSizei(l) .GT. ppiclf_bins_dx(l)) THEN
            PRINT*, "ppiclf_nid", ppiclf_nid
            PRINT*,'EleSizei', EleSizei(l), '>', ppiclf_bins_dx(l),
     >         'ppiclf_bins_dx in MapOverlapMesh'
            PRINT*, 'Dimension:',l
            PRINT*, 'ppiclf_bins_dx(1:3)', ppiclf_bins_dx(1:3)
            PRINT*, 'ppiclf_n_bins(1:3)', ppiclf_n_bins(1:3)
            PRINT*, " " 
            CALL ppiclf_exittr('Error: EleSizei(OverlapMesh)',0.0D0,l)
          END IF
        END DO !l 

        ! Fluid Cell vertex position without additional length
        rxval = centeri(1)
        ryval = centeri(2)
        rzval = 0.0D0
        IF(ppiclf_ndim .GT. 2) rzval = centeri(3)
      
        ! Exits if fluid cell vertex is outside of all bin 
        ! boundaries + Exchange Ghost Fluid Cell Buffer (EleSizei)
        IF (rxval .GT. (ppiclf_binb(2))) CYCLE
        IF (rxval .LT. (ppiclf_binb(1))) CYCLE
        IF (ryval .GT. (ppiclf_binb(4))) CYCLE
        IF (ryval .LT. (ppiclf_binb(3))) CYCLE
        IF (ppiclf_ndim .GT. 2 .AND. rzval .GT. 
     >      (ppiclf_binb(6))) CYCLE
        IF (ppiclf_ndim.GT.2 .AND. rzval .LT.
     >      (ppiclf_binb(5))) CYCLE
 
        ! Determines what bin the fluid cell is nominally mapped to
        ii    = FLOOR((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
        jj    = FLOOR((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
        kk    = FLOOR((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3))

        ! Default is Do loop with ix=iy=iz=2 for fluid cells not near
        ! bin boundary

        ixLow =2
        ixHigh=2
        iyLow =2
        iyHigh=2
        izLow =2
        izHigh=2

        ! These series of if statements check if bin mapping changes
        ! when adding/subtracting multiple of fluid cell length defined
        ! by EleSizei(l). 
        ! This is used to map fluid cells slightly outside of the ppiclf
        ! bin boundary.  If any .NE. 2, then fluid cell is mapped to
        ! multiple ppiclf bins. 
        
        IF (FLOOR((rxval + EleSizei(1) - ppiclf_binb(1))
     >       /ppiclf_bins_dx(1)) .NE. ii)  ixHigh = 3

        IF (FLOOR((rxval - EleSizei(1) - ppiclf_binb(1))
     >       /ppiclf_bins_dx(1)) .NE. ii)  ixLow = 1
        IF (FLOOR((ryval + EleSizei(2) - ppiclf_binb(3))
     >       /ppiclf_bins_dx(2)) .NE. jj)  iyHigh = 3

        IF (FLOOR((ryval - EleSizei(2) - ppiclf_binb(3))
     >       /ppiclf_bins_dx(2)) .NE. jj)  iyLow = 1

        IF (ppiclf_ndim .GT. 2 .AND. FLOOR((rzval + EleSizei(3)
     >    - ppiclf_binb(5))/ppiclf_bins_dx(3)) .NE. kk)  izHigh = 3

        IF (ppiclf_ndim .GT. 2 .AND. FLOOR((rzval - EleSizei(3)
     >    - ppiclf_binb(5))/ppiclf_bins_dx(3)) .NE. kk)  izLow = 1

        DO ix=ixLow,ixHigh
          DO iy=iyLow,iyHigh
            DO iz=izLow,izHigh
              ! Change cell position by EleSizei if ix,iy,or iz NE 2
              rxval = centeri(1) + (ix-2)*EleSizei(1)
              ryval = centeri(2) + (iy-2)*EleSizei(2)
              rzval = 0.0D0
              IF(ppiclf_ndim.GT.2) rzval = centeri(3)
     >                + (iz-2)*EleSizei(3)
              ! Find bin for adjusted rval
              ii    = FLOOR((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
              jj    = FLOOR((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
              kk    = FLOOR((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
              IF (ppiclf_ndim.LT.3) kk = 0
              

              ! This covers ghost exchanged cells for linear periodicity
              ! Maps cells greater than ppiclf bin domain to first bin
              ! Maps cells less than ppiclf bin domain to last bin
              IF (x_per_flag .EQ. 1) THEN
                IF (ii .EQ. ppiclf_n_bins(1)) ii = 0
                IF (ii .EQ. -1) ii = ppiclf_n_bins(1) - 1
              END IF
              IF (y_per_flag .EQ. 1) THEN
                IF (jj .EQ. ppiclf_n_bins(2)) jj = 0
                IF (jj .EQ. -1) jj = ppiclf_n_bins(2) - 1
              END IF
              IF (z_per_flag .EQ. 1) THEN
                IF (kk .EQ. ppiclf_n_bins(3)) kk = 0
                IF (kk .EQ. -1) kk = ppiclf_n_bins(3) - 1
              END IF
              
              ! Ensures duplicate cells don't get sent to same processor
              IF (ii .LT. 0 .OR. ii .GT. ppiclf_n_bins(1)-1) CYCLE
              IF (jj .LT. 0 .OR. jj .GT. ppiclf_n_bins(2)-1) CYCLE
              IF (kk .LT. 0 .OR. kk .GT. ppiclf_n_bins(3)-1) CYCLE


              ! Calculates processor rank
              ndum  = ii + ppiclf_n_bins(1)*jj + 
     >                     ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
              nrank = ndum

                            ppiclf_neltb = ppiclf_neltb + 1
              IF(ppiclf_neltb .GT. PPICLF_LEE) THEN
                PRINT*, '***ERROR*** PPICLF_LEE',PPICLF_LEE, 'in', 
     >           'MapOverlapMesh must be greater than', ppiclf_neltb 
                CALL ppiclf_exittr('Increase PPICLF_LEE$ (MapOverlap)',0.0D0
     >               ,ppiclf_neltb)
              END IF
              ! Stores element to rank mapping.
              ppiclf_er_map(1,ppiclf_neltb) = ie
              ppiclf_er_map(2,ppiclf_neltb) = ppiclf_nid
              ppiclf_er_map(3,ppiclf_neltb) = ndum
              ppiclf_er_map(4,ppiclf_neltb) = nrank
              ppiclf_er_map(5,ppiclf_neltb) = nrank
              ppiclf_er_map(6,ppiclf_neltb) = nrank

!              The loop makes this subroutine 10x slower
!              Replaced with tempCheck below, since cell would 
!              be duplicated in sequential order. It shouldn't happen,
!              so implemented as error vs standard fix in loop.

!              IF (ppiclf_neltb .GT. 1) THEN
!              DO il=1,ppiclf_neltb-1
!                 IF (ppiclf_er_map(1,il) .EQ. ie) THEN
!                 IF (ppiclf_er_map(4,il) .EQ. nrank) THEN
!                    PRINT*, 'AVERY - NELTB Loop remover still used!'
!                    ppiclf_neltb = ppiclf_neltb - 1
!                     CYCLE
!                 END IF
!                 END IF
!              END DO
!              END IF

            END DO !iz
          END DO !iy
        END DO !ix
      END DO !ie
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      DO ie=1,ppiclf_neltb !Number of fluid cells in ppiclf bin domain
       ! These copy all nxyz vertecies since Fortran is column-major
       iee = ppiclf_er_map(1,ie)
       CALL ppiclf_copy(ppiclf_xm1b(1,1,1,1,ie)
     >                 ,ppiclf_xm1bs(1,1,1,1,iee),nxyz)
       CALL ppiclf_copy(ppiclf_xm1b(1,1,1,2,ie)
     >                 ,ppiclf_xm1bs(1,1,1,2,iee),nxyz)
       CALL ppiclf_copy(ppiclf_xm1b(1,1,1,3,ie)
     >                 ,ppiclf_xm1bs(1,1,1,3,iee),nxyz)
      END DO

      ppiclf_neltbb = ppiclf_neltb
      DO ie=1,ppiclf_neltbb
         ! Copies element to rank mapping (integer copy)
         CALL ppiclf_icopy(ppiclf_er_maps(1,ie),ppiclf_er_map(1,ie)
     >             ,PPICLF_LRMAX)
      END DO

      ! GSLIB required info
      ! neltb - number of columns to transfer
      ! PPICLF_LEE - number of columns declared
      ! nl - partl row size (dummy logical variable)
      nl   = 0
      ! nii - ppiclf_er_maps row size declared
      nii  = PPICLF_LRMAX
      ! njj - Row index of ppiclf_er_maps with processor/rank number
      njj  = 6
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      ! nrr - ppiclf_xm1b row size declared
      nrr  = nxyz*3
      ! Defines sorting order
      nkey(1) = 2
      nkey(2) = 1

      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl,ppiclf_neltb
     >     ,PPICLF_LEE,ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,njj)
      CALL pfgslib_crystal_tuple_sort    (ppiclf_cr_hndl,ppiclf_neltb
     >     ,ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,nkey,2)

!*************
! This is only needed for multi-element projection.
! We are currently doing single element projection (hardcoded in rocpicl)
! 
!      DO ie=1,ppiclf_neltb
!      DO k=1,PPICLF_LEZ
!      DO j=1,PPICLF_LEY
!      DO i=1,PPICLF_LEX
!         rxval = ppiclf_xm1b(i,j,k,1,ie)
!         ryval = ppiclf_xm1b(i,j,k,2,ie)
!         rzval = 0.0D0
!         IF(ppiclf_ndim.GT.2) rzval = ppiclf_xm1b(i,j,k,3,ie)
!         
!         ii    = FLOOR((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
!         jj    = FLOOR((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
!         kk    = FLOOR((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
!         IF (ppiclf_ndim.EQ.2) kk = 0
!          IF (ii .EQ. ppiclf_n_bins(1)) ii = ppiclf_n_bins(1) - 1
!          IF (jj .EQ. ppiclf_n_bins(2)) jj = ppiclf_n_bins(2) - 1
!          IF (kk .EQ. ppiclf_n_bins(3)) kk = ppiclf_n_bins(3) - 1
!          IF (ii .EQ. -1) ii = 0
!          IF (jj .EQ. -1) jj = 0
!          IF (kk .EQ. -1) kk = 0
!          ndum  = ii + ppiclf_n_bins(1)*jj + 
!     >                 ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
!
!         ppiclf_modgp(i,j,k,ie,1) = ii
!         ppiclf_modgp(i,j,k,ie,2) = jj
!         ppiclf_modgp(i,j,k,ie,3) = kk
!         ppiclf_modgp(i,j,k,ie,4) = ndum
!   
!      END DO
!      END DO
!      END DO
!      END DO
!**************

      DO ie=1,ppiclf_neltb
         ! Finds minimum and maximum vertex in x,y,z of cell
         ppiclf_xerange(1,1,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,1,ie),nxyz)
         ppiclf_xerange(2,1,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,1,ie),nxyz)
         ppiclf_xerange(1,2,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,2,ie),nxyz)
         ppiclf_xerange(2,2,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,2,ie),nxyz)
         ppiclf_xerange(1,3,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,3,ie),nxyz)
         ppiclf_xerange(2,3,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,3,ie),nxyz)
         
         ! Finds the ppiclf bin that the max/min cell vertex resides in
         ilow  = 
     >     FLOOR((ppiclf_xerange(1,1,ie) - ppiclf_binb(1))/
     >                                             ppiclf_bins_dx(1))
         ihigh = 
     >     FLOOR((ppiclf_xerange(2,1,ie) - ppiclf_binb(1))/
     >                                             ppiclf_bins_dx(1))
         jlow  = 
     >     FLOOR((ppiclf_xerange(1,2,ie) - ppiclf_binb(3))/
     >                                             ppiclf_bins_dx(2))
         jhigh = 
     >     FLOOR((ppiclf_xerange(2,2,ie) - ppiclf_binb(3))/
     >                                             ppiclf_bins_dx(2))
         klow  = 
     >     FLOOR((ppiclf_xerange(1,3,ie) - ppiclf_binb(5))/
     >                                             ppiclf_bins_dx(3))
         khigh = 
     >     FLOOR((ppiclf_xerange(2,3,ie) - ppiclf_binb(5))/
     >                                             ppiclf_bins_dx(3))
         IF (ppiclf_ndim.LT.3) THEN
            klow = 0
            khigh = 0
         END IF

         ! Maps the cell to bin rank range (1,2) and min/max bins in
         ! x,y,z (3-8).  If ppiclf_el_map(1:8,ie) are the same, then 
         ! fluid cell is only in 1 bin.
         ppiclf_el_map(1,ie) = ilow  + ppiclf_n_bins(1)*jlow  
     >                         + ppiclf_n_bins(1)*ppiclf_n_bins(2)*klow
         ppiclf_el_map(2,ie) = ihigh + ppiclf_n_bins(1)*jhigh 
     >                         + ppiclf_n_bins(1)*ppiclf_n_bins(2)*khigh
         ppiclf_el_map(3,ie) = ilow
         ppiclf_el_map(4,ie) = ihigh
         ppiclf_el_map(5,ie) = jlow
         ppiclf_el_map(6,ie) = jhigh
         ppiclf_el_map(7,ie) = klow
         ppiclf_el_map(8,ie) = khigh
      END DO

      IF (icalld .EQ. 0) THEN 
         icalld = icalld + 1
         CALL ppiclf_prints('   *Begin mpi_comm_split$')
            CALL mpi_comm_split(ppiclf_comm
     >                         ,ppiclf_nid
     >                         ,0
     >                         ,ppiclf_comm_nid
     >                         ,ierr)
         CALL ppiclf_prints('    End mpi_comm_split$')
         CALL ppiclf_io_OutputDiagGrid
      END IF

      RETURN
      END 
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_comm_InitOverlapMesh(ncell,lx1,ly1,lz1,
     >                                       xgrid,ygrid,zgrid)
     > bind(C, name="ppiclc_comm_InitOverlapMesh")
#else
      subroutine ppiclf_comm_InitOverlapMesh(ncell,lx1,ly1,lz1,
     >                                       xgrid,ygrid,zgrid)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 ncell
      integer*4 lx1
      integer*4 ly1
      integer*4 lz1
      real*8    xgrid(*)
      real*8    ygrid(*)
      real*8    zgrid(*)
!
! External:
!
      integer*4 nxyz, i, j, ie
      integer*4 k, jj, icont
!
      ppiclf_overlap = .true.

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitOverlap$',0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitOverlap$'
     >                  ,0.0d0,0)

      if (ncell .gt. PPICLF_LEE .or. ncell .lt. 0) then
        PRINT*, '***ERROR*** PPICLF_LEE', PPICLF_LEE, 'in', 
     > 'InitMapOverlapMesh must be greater than', ncell 
        call ppiclf_exittr('Increase LEE in InitOverlap$',0.0d0,ncell)
      endif
      if (lx1 .ne. PPICLF_LEX) 
     >   call ppiclf_exittr('LX1 != LEX in InitOverlap$',0.0d0,ncell)
      if (ly1 .ne. PPICLF_LEY)
     >   call ppiclf_exittr('LY1 != LEY in InitOverlap$',0.0d0,ncell)
      if (lz1 .ne. PPICLF_LEZ)
     >   call ppiclf_exittr('LZ1 != LEZ in InitOverlap$',0.0d0,ncell)

      ppiclf_nee = ncell
      nxyz       = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ

      do ie=1,ppiclf_nee
         ! TLJ changing loop structure
         !do i=1,nxyz
         !   j = (ie-1)*nxyz + i
         !   ppiclf_xm1bs(i,1,1,1,ie) = xgrid(j)
         !   ppiclf_xm1bs(i,1,1,2,ie) = ygrid(j)
         !   ppiclf_xm1bs(i,1,1,3,ie) = zgrid(j)
         !enddo
         icont = 0
         do k=1,PPICLF_LEZ
         do j=1,PPICLF_LEY
         do i=1,PPICLF_LEX
            icont = icont + 1
            jj = (ie-1)*nxyz + icont
            ppiclf_xm1bs(i,j,k,1,ie) = xgrid(jj)
            ppiclf_xm1bs(i,j,k,2,ie) = ygrid(jj)
            ppiclf_xm1bs(i,j,k,3,ie) = zgrid(jj)
         enddo
         enddo
         enddo
      enddo
      
      call ppiclf_solve_InitSolve

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_FindParticle
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 ix, iy, iz, i, ii, jj, kk, ndum, nrank
!
      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim.eq.3)
     >iz = 3

      do i=1,ppiclf_npart
         ! check if particles are greater or less than binb bounds....
         ii  = floor((ppiclf_y(ix,i)-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
         jj  = floor((ppiclf_y(iy,i)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
         kk  = floor((ppiclf_y(iz,i)-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
         if (ppiclf_ndim .lt. 3) kk = 0
         ndum  = ii + ppiclf_n_bins(1)*jj + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
         nrank = ndum
         
         print*, "FindParticle, nid, ppiclf_binb, bins_dx", 
     >    ppiclf_nid, ppiclf_binb, ppiclf_bins_dx
         print*, "FindParticle, nid, ii, jj, kk, nrank",
     >           ppiclf_nid, ii, jj, kk, nrank        

         ppiclf_iprop(8,i)  = ii
         ppiclf_iprop(9,i)  = jj
         ppiclf_iprop(10,i) = kk
         ppiclf_iprop(11,i) = ndum

         ppiclf_iprop(3,i)  = nrank ! where particle is actually moved
         ppiclf_iprop(4,i)  = nrank ! where particle is actually moved
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_MoveParticle
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      logical partl    
      integer*4 lrf
      parameter(lrf = PPICLF_LRS*4 + PPICLF_LRP + PPICLF_LRP2
     >       + PPICLF_LRP3)
      real*8 rwork(lrf,PPICLF_LPART)
      integer*4 i, ic, j0
!

      print*, "Move Particle - copying data"

      do i=1,ppiclf_npart
         ic = 1
         call ppiclf_copy(rwork(ic,i),ppiclf_y(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_y1((i-1)*PPICLF_LRS+1)
     >                   ,PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_ydot(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_ydotc(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop(1,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop2(1,i),PPICLF_LRP2)
         ic = ic + PPICLF_LRP2
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop3(1,i),PPICLF_LRP3)
      enddo

      print*, "Move Particle - calling pfgslib"

      j0 = 4

      print*, ppiclf_nid, ppiclf_cr_hndl, ppiclf_npart,PPICLF_LPART,
     >        PPICLF_LIP, partl, lrf, j0
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >                                  ,ppiclf_npart,PPICLF_LPART
     >                                  ,ppiclf_iprop,PPICLF_LIP
     >                                  ,partl,0
     >                                  ,rwork,lrf
     >                                  ,j0)

      print*, "Move Particle - Done calling gslib"
      if (ppiclf_npart .gt. PPICLF_LPART .or. ppiclf_npart .lt. 0) then
        print*, "***ERROR In MoveParticle, PPICLF_LPART", 
     >   PPICLF_LPART, "smaller than ", ppiclf_npart
        call ppiclf_exittr('Increase LPART$',0.0d0,ppiclf_npart)
      endif

      do i=1,ppiclf_npart
         ic = 1
         call ppiclf_copy(ppiclf_y(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_y1((i-1)*PPICLF_LRS+1),rwork(ic,i)
     >                   ,PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_ydot(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_ydotc(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_rprop(1,i),rwork(ic,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call ppiclf_copy(ppiclf_rprop2(1,i),rwork(ic,i),PPICLF_LRP2)
         ic = ic + PPICLF_LRP2
         call ppiclf_copy(ppiclf_rprop3(1,i),rwork(ic,i),PPICLF_LRP3)
      enddo
        
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 rxval, ryval, rzval, rxl, ryl, rzl, rxr, ryr, rzr, 
     >       distchk, dist, map(PPICLF_LRP_PRO)
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >           gpsave(27), nfacegp, nedgegp, ncornergp, 
     >           ip, idum, iip, jjp, kkp, ii1, jj1, kk1, iig, jjg, kkg,
     >           isave, ndumn, nrank,  i, j, k, ifc, ist
!

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      ppiclf_npart_gp = 0

      do ip=1,ppiclf_npart

         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))

         idum = 0
         do j=1,PPICLF_LRS
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip)
         enddo
         idum = PPICLF_LRS
         do j=1,PPICLF_LRP
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip)
         enddo
         idum = PPICLF_LRS+PPICLF_LRP
         do j=1,PPICLF_LRP_PRO
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = map(j)
         enddo

         rxval = ppiclf_cp_map(PPICLF_JX,ip) ! ppiclf_y(PPICLF_JX,ip)
         ryval = ppiclf_cp_map(PPICLF_JY,ip) ! ppiclf_y(PPICLF_JY,ip)
         rzval = 0.0d0
         if (ppiclf_ndim .gt. 2) rzval = ppiclf_cp_map(PPICLF_JZ,ip) ! ppiclf_y(PPICLF_JZ,ip)

         iip    = ppiclf_iprop(8,ip) ! i-index of bin where ip particle is located
         jjp    = ppiclf_iprop(9,ip) ! j-index of bin where ip particle is located
         kkp    = ppiclf_iprop(10,ip)! k-index of bin where ip particle is located

         rxl = ppiclf_binb(1) + ppiclf_bins_dx(1)*iip ! min-x of bin where ip is located
         rxr = rxl + ppiclf_bins_dx(1)                ! max-x of bin where ip is located
         ryl = ppiclf_binb(3) + ppiclf_bins_dx(2)*jjp ! min-y of bin where ip is located
         ryr = ryl + ppiclf_bins_dx(2)                ! max-y of bin where ip is located
         rzl = 0.0d0
         rzr = 0.0d0
         if (ppiclf_ndim .gt. 2) then
            rzl = ppiclf_binb(5) + ppiclf_bins_dx(3)*kkp ! min-z of bin where ip is located
            rzr = rzl + ppiclf_bins_dx(3)                ! max-z of bin where ip is located
         endif

         isave = 0

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) ! i-index of bin face  
            jj1 = jjp + el_face_num(ist+2) ! j-index of bin face
            kk1 = kkp + el_face_num(ist+3) ! k-index of bin face

            iig = ii1
            jjg = jj1
            kkg = kk1

            ! Step 1 - Check if particle is within user-input max(neighbor-width, neighbor-width) from bin boundary
            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2     ! calculate x-distance check based on user-input neighbor-width
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2 ! calculate x-distance of ip particle from x-bin boundary
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2 ! calculate y-distance of ip particle from y-bin boundary
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2 ! calculate z-distance of ip particle from z-bin boundary
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
              endif
            endif
            distchk = sqrt(distchk) ! calculate distance magnitude check based on neighbor-width
            dist = sqrt(dist)    ! calculate distance magnitude of ip particle from bin boundary

            ! if ip particle is farther away from bin boundary than the user-input neighbor-width -> dont create ghost
            if (dist .gt. distchk) cycle

            ! Step 2 - Check whether particle is about to leave the periodic domain

            ! David previously used to fix the bin boundaries as big as the periodic domain 
            ! and evaluate the ghost particle using the 3 if statements below. 

            ! This is now commented out as we check for periodicity in a separate subroutine 
            ! by initializing linear periodic plane and evaluating the distance of the real particle
            ! relative to the periodic planes. 

            ! David's ghost algorithm is just used for mapping the real particle (creating ghosts for it)
            ! when needed in the nearby neighboring bins. 

            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
                cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
                cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
                cycle
            endif

            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            ! if ghost rank is the same processor of the real particle -> do not create ghost
            if (nrank .eq. ppiclf_nid) cycle

            ! gpsave keeps track of which nrank destinations we've already created ghosts for
            ! if a ghost for destination nrank was already created earlier, and this is not a boundary crossing -> don't duplicate
            ! it.
            do i=1,isave
               if (gpsave(i) .eq. nrank) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            do k=1,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
!            write(1920+ppiclf_nid,*) "DavidFaces", ppiclf_time, ! 0-1   
!     >              ppiclf_nid, nrank, ppiclf_npart_gp,       ! 2-4
!     >              ppiclf_rprop_gp(1:6, ppiclf_npart_gp),    ! 5-7, 8,10
!     >              ppiclf_y(1:3,ip), ppiclf_y(4:6,ip)        ! 11-13, 14-16
  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
                cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
                cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
                cycle
            endif

            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            do k=1,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
!            write(1920+ppiclf_nid,*) "DavidEdges", ppiclf_time, ! 0-1   
!     >              ppiclf_nid, nrank, ppiclf_npart_gp,       ! 2-4
!     >              ppiclf_rprop_gp(1:6, ppiclf_npart_gp),    ! 5-7, 8,10
!     >              ppiclf_y(1:3,ip), ppiclf_y(4:6,ip)        ! 11-13, 14-16
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
                cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
                cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
                cycle
            endif

            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank) goto 333
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            do k=1,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
!            write(1920+ppiclf_nid,*) "DavidCorners", ppiclf_time, ! 0-1   
!     >              ppiclf_nid, nrank, ppiclf_npart_gp,       ! 2-4
!     >              ppiclf_rprop_gp(1:6, ppiclf_npart_gp),    ! 5-7, 8,10
!     >              ppiclf_y(1:3,ip), ppiclf_y(4:6,ip)        ! 11-13, 14-16
  333 continue
         enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateLinearGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 rxnew(3), rxval, ryval, rzval, rxl, ryl, rzl, rxr, ryr, 
     >       rzr, distchk, dist, distmin(3), distmax(3), 
     >       map(PPICLF_LRP_PRO)
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >           gpsave(27), nfacegp, nedgegp, ncornergp, 
     >           iperiodicx, iperiodicy, iperiodicz, 
     >           ip, idum, iip, jjp, kkp, ii1, jj1, kk1, iig, jjg, kkg,
     >           isave, ndumn, nrank, i, j, k, ifc, ist

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------

      ! Initialize dist as crazy values for periodicity check
      distmin(1) = 1E20
      distmin(2) = 1E20
      distmin(3) = 1E20
      distmax(1) = 1E20
      distmax(2) = 1E20
      distmax(3) = 1E20

      do ip=1,ppiclf_npart

        ! Step 1 - Get the mapped particle properties
         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))

         idum = 0
         do j=1,PPICLF_LRS
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip)
         enddo
         idum = PPICLF_LRS
         do j=1,PPICLF_LRP
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip)
         enddo
         idum = PPICLF_LRS+PPICLF_LRP
         do j=1,PPICLF_LRP_PRO
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = map(j)
         enddo

         ! Step 2 - Evalute distance of particle relative to linear periodic planes
         !          distmin -> distance between particle and min periodic planes             
         !          distmax -> distance between particle and max periodic planes              

         rxval = ppiclf_cp_map(PPICLF_JX,ip) ! ppiclf_y(PPICLF_JX,ip)
         ryval = ppiclf_cp_map(PPICLF_JY,ip) ! ppiclf_y(PPICLF_JY,ip)
         rzval = 0.0d0
         if (ppiclf_ndim .gt. 2) rzval = ppiclf_cp_map(PPICLF_JZ,ip) ! ppiclf_y(PPICLF_JZ,ip)
         rxnew(1) = rxval
         rxnew(2) = ryval
         rxnew(3) = rzval

         if(iperiodicx .eq. 0) then
          distmin(1) = abs(A_xmin*rxval + B_xmin*ryval 
     >                   + C_xmin*rzval + D_xmin)
          distmin(1) = distmin(1) / 
     >                 sqrt(A_xmin**2 + B_xmin**2 + C_xmin**2)

          distmax(1) = abs(A_xmax*rxval + B_xmax*ryval
     >                   + C_xmax*rzval + D_xmax)
          distmax(1) = distmax(1) / 
     >                 sqrt(A_xmax**2 + B_xmax**2 + C_xmax**2)
         endif

         if(iperiodicy .eq. 0) then
          distmin(2) = abs(A_ymin*rxval + B_ymin*ryval 
     >                   + C_ymin*rzval + D_ymin)
          distmin(2) = distmin(2) / 
     >                 sqrt(A_ymin**2 + B_ymin**2 + C_ymin**2)

          distmax(2) = abs(A_ymax*rxval + B_ymax*ryval 
     >                   + C_ymax*rzval +D_ymax)
          distmax(2) = distmax(2) / 
     >                 sqrt(A_ymax**2 + B_ymax**2 + C_ymax**2)
         endif

         if(iperiodicz .eq. 0) then
          distmin(3) = abs(A_zmin*rxval + B_zmin*ryval 
     >                   + C_zmin*rzval + D_zmin)
          distmin(3) = distmin(3) /
     >                 sqrt(A_zmin**2 + B_zmin**2 + C_zmin**2)

          distmax(3) = abs(A_zmax*rxval + B_zmax*ryval 
     >                   + C_zmax*rzval + D_zmax)
          distmax(3) = distmax(3) / 
     >                 sqrt(A_zmax**2 + B_zmax**2 + C_zmax**2)
         endif

         distchk = ppiclf_d2chk(1)
        
         ! particle far away from all periodic planes 
         if((distmin(1).gt.distchk).and.(distmax(1).gt.distchk).and.
     >      (distmin(2).gt.distchk).and.(distmax(2).gt.distchk).and.
     >      (distmin(3).gt.distchk).and.(distmax(3).gt.distchk))then
         cycle
         ! particle close to at least one periodic wall
         elseif((distmin(1).le.distchk).or.(distmax(1).le.distchk).or.
     >          (distmin(2).le.distchk).or.(distmax(2).le.distchk).or.
     >          (distmin(3).le.distchk).or.(distmax(3).le.distchk))then
         call ppiclf_comm_CheckLinearPeriodic(distmin,distmax,rxnew)

         else
           print*, "***ERROR Ghost Periodic Linear Plane!"
           print*, "distmin =", distmin
           print*, "distmax =", distmax
           print*, "distchk =", distchk
           call ppiclf_exittr('Error Ghost Periodic 
     >                                 Linear Plane!')
         endif

         ! Step 3 - map this linearly periodic mirror particle to nearby faces, edges, and corners.
         !        - DO NOT ADD MIRROR PARTICLE TO GHOST LIST FOR LINEAR PERIODICITY

         iip    = floor((rxnew(1)-ppiclf_binb(1))/ppiclf_bins_dx(1)) ! i-index of bin where mirror particle is located
         jjp    = floor((rxnew(2)-ppiclf_binb(3))/ppiclf_bins_dx(2)) ! j-index of bin where mirror particle is located
         kkp    = floor((rxnew(3)-ppiclf_binb(5))/ppiclf_bins_dx(3)) ! k-index of bin where mirror particle is located

         rxval = rxnew(1)
         ryval = rxnew(2)
         rzval = rxnew(3)

         rxl = ppiclf_binb(1) + ppiclf_bins_dx(1)*iip ! min-x of bin where ip is located
         rxr = rxl + ppiclf_bins_dx(1)                ! max-x of bin where ip is located
         ryl = ppiclf_binb(3) + ppiclf_bins_dx(2)*jjp ! min-y of bin where ip is located
         ryr = ryl + ppiclf_bins_dx(2)                ! max-y of bin where ip is located
         rzl = 0.0d0
         rzr = 0.0d0
         if (ppiclf_ndim .gt. 2) then
            rzl = ppiclf_binb(5) + ppiclf_bins_dx(3)*kkp ! min-z of bin where ip is located
            rzr = rzl + ppiclf_bins_dx(3)                ! max-z of bin where ip is located
         endif

         isave = 0

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) ! i-index of bin face  
            jj1 = jjp + el_face_num(ist+2) ! j-index of bin face
            kk1 = kkp + el_face_num(ist+3) ! k-index of bin face

            iig = ii1
            jjg = jj1
            kkg = kk1

            ! Step 1 - Check if particle is within user-input neighbor-width from bin boundary
            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2     ! calculate x-distance check based on user-input neighbor-width
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2 ! calculate x-distance of ip particle from x-bin boundary
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2 ! calculate y-distance of ip particle from y-bin boundary
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2 ! calculate z-distance of ip particle from z-bin boundary
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk) ! calculate distance magnitude check based on neighbor-width
            dist = sqrt(dist)    ! calculate distance magnitude of ip particle from bin boundary

            ! if ip particle is farther away from bin boundary than the user-input neighbor-width -> dont create ghost
            if (dist .gt. distchk) cycle

            ! Step 2 - Check whether particle is outside the bi

            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
                cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
                cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
                cycle
            endif

            ! if particle is crossing the periodic domain in each direction -> iflgsum = 1 crossing in 1-direction
            !                                                     otherwise -> iflgsum = 0
            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            ! if ghost rank is the same processor of the real particle -> do not create a ghost
            if (nrank .eq. ppiclf_nid) cycle

            ! gpsave keeps track of which nrank destinations we've already created ghosts for
            ! if a ghost for destination nrank was already created earlier -> don't duplicate
            ! it.
            do i=1,isave
               if (gpsave(i) .eq. nrank) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ! Step 3 - Modify ghost properties accordingly
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(PPICLF_JX,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(PPICLF_JY,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(PPICLF_JZ,ppiclf_npart_gp) = rxnew(3)

            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
!            write(1920+ppiclf_nid,*) "LinearFaces", ppiclf_time,                ! 0-1   
!     >              ip, ifc, ist, nrank,                           ! 2-5 
!     >              ii1, jj1, kk1, dist, distchk,                  ! 6-10
!     >              ppiclf_npart_gp,                               ! 11
!     >              iig, jjg, kkg, ndumn,                          ! 12-15
!     >              rxnew(1:3), ppiclf_y(4:6,ip), ppiclf_cp_map(4:6,ip),  ! 16-24
!     >              ppiclf_rprop_gp(4:6, ppiclf_npart_gp),         ! 25-27
!     >              ppiclf_y(1:3,ip), ppiclf_y(4:6,ip),            ! 28-30, 31-33
!     >              ppiclf_n_bins, ppiclf_bins_dx,                 ! 34-36, 37-39
!     >              ppiclf_binb,                                   ! 40-45
!     >              ppiclf_nid                                     ! 46
  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
              cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
              cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
              cycle
            endif

            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(PPICLF_JX,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(PPICLF_JY,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(PPICLF_JZ,ppiclf_npart_gp) = rxnew(3)

            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
!            write(1920+ppiclf_nid,*) "LinearEdges", ppiclf_time,                ! 0-1   
!     >              ip, ifc, ist, nrank,                           ! 2-5 
!     >              ii1, jj1, kk1, dist, distchk,                  ! 6-10
!     >              ppiclf_npart_gp,                               ! 11
!     >              iig, jjg, kkg, ndumn,                          ! 12-15
!     >              rxnew(1:3), ppiclf_y(4:6,ip), ppiclf_cp_map(4:6,ip),  ! 16-24
!     >              ppiclf_rprop_gp(4:6, ppiclf_npart_gp),         ! 25-27
!     >              ppiclf_y(1:3,ip), ppiclf_y(4:6,ip)             ! 28-30, 31-33
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
              cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
              cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
              cycle
            endif

            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank) goto 333
            enddo
            gpsave(isave) = nrank
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(PPICLF_JX,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(PPICLF_JY,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(PPICLF_JZ,ppiclf_npart_gp) = rxnew(3)

            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
!            write(1920+ppiclf_nid,*) "LinearCorners", ppiclf_time, ! 0-1   
!     >              ppiclf_nid, nrank, ppiclf_npart_gp,       ! 2-4
!     >              ppiclf_rprop_gp(1:6, ppiclf_npart_gp),    ! 5-7, 8,10
!     >              ppiclf_y(1:3,ip), ppiclf_y(4:6,ip)        ! 11-13, 14-16
  333 continue
         enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateAngularGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 rval(3), rxval, ryval, rzval, 
     >       rxl, ryl, rzl, rxr, ryr, rzr, distchk, dist,
     >       dist1, dist2, map(PPICLF_LRP_PRO) 
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >           nfacegp, nedgegp, ncornergp,
     >           ip, idum, iip, jjp, kkp, ii1, jj1, kk1, isave, 
     >           iig, jjg, kkg, ndumn, nrank, i, ifc, ist, j, k,
     >           gpsave(27)
      logical CW
!

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      CW = .false. ! logical statement to use CW or CCW rotation matrix
! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      do ip=1,ppiclf_npart

        ! Step 1 - Get the mapped particle properties
         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))

         idum = 0
         do j=1,PPICLF_LRS
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip)
         enddo
         idum = PPICLF_LRS
         do j=1,PPICLF_LRP
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip)
         enddo
         idum = PPICLF_LRS+PPICLF_LRP
         do j=1,PPICLF_LRP_PRO
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = map(j)
         enddo
         
         rval(1) = ppiclf_cp_map(PPICLF_JX,ip)
         rval(2) = ppiclf_cp_map(PPICLF_JY,ip)
         rval(3) = ppiclf_cp_map(PPICLF_JZ,ip)

         rxval = rval(1)
         ryval = rval(2)
         rzval = rval(3)

         ! Step 2 - Evalute distance of particle relative to angular periodic planes
         !          dist1 -> distance to lower plane
         !          dist2 -> distance to upper plane

          dist1 = abs(Ap1*rxval + Bp1*ryval + Cp1*rzval)
          dist1 = dist1/sqrt(Ap1**2 + Bp1**2 + Cp1**2)

          dist2 = abs(Ap2*rxval + Bp2*ryval + Cp2*rzval)
          dist2 = dist2/sqrt(Ap2**2 + Bp2**2 + Cp2**2)


         distchk = ppiclf_d2chk(1)
         ! particle farther from both angular planes
         if((dist1.gt.distchk) .and. (dist2.gt.distchk)) then
!           print*, "dist1 dist2 cycling"
!           print*, "dist1 =", dist1
!           print*, "dist2 =", dist2
!           print*, "distchk =", distchk
           cycle
         ! particle closer to lower angular periodic plane -> rotate CCW
         elseif((dist1.gt.distchk) .and. (dist2.lt.distchk)) then
           rval = MATMUL(rotCCW, rval)
!           print*, "Particle Closer to Lower Angular Plane"

         ! particle closer to upper angular periodic plane -> rotate CW
         elseif((dist1.lt.distchk) .and. (dist2.gt.distchk)) then
           rval = MATMUL(rotCW, rval)
           CW = .true.
!           print*, "Particle Closer to Upper Angular Plane"
         else
           print*, "***ERROR Ghost Periodic Angular Plane!"
           print*, "dist1 =", dist1
           print*, "dist2 =", dist2
           print*, "distchk =", distchk
           call ppiclf_exittr('Error Ghost Periodic 
     >                                 Angular Plane!')
         endif

         ! These are the angularly rotated coordinates of real particle
         !                           i.e. coordinates of the mirror particle

         rxval = rval(1)
         ryval = rval(2)
         rzval = rval(3)

         ! Step 3 - check if mirror particle is within bin bounds

         iip = floor((rval(1)-ppiclf_binb(1))/ppiclf_bins_dx(1)) ! i-index of bin where mirror particle is located
         jjp = floor((rval(2)-ppiclf_binb(3))/ppiclf_bins_dx(2)) ! j-index of bin where mirror particle is located
         kkp = floor((rval(3)-ppiclf_binb(5))/ppiclf_bins_dx(3)) ! k-index of bin where mirror particle is located

         iig = iip
         jjg = jjp
         kkg = kkp

         ! If mirror particle outside of bin bounds -> cycle and do not add to ghost list
         if (iig .lt. -1 .or. iig .gt. ppiclf_n_bins(1)) then
             cycle
         endif
         if (jjg .lt. -1 .or. jjg .gt. ppiclf_n_bins(2)) then
             cycle
         endif
         if (kkg .lt. -1 .or. kkg .gt. ppiclf_n_bins(3)) then
             cycle
         endif
         
         ! Step 4 - add the mirror particle as a ghost on the ghost list since within bin bounds

         ndumn = iig + ppiclf_n_bins(1)*jjg
     >               + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg

         nrank = ndumn

         ppiclf_npart_gp = ppiclf_npart_gp + 1
         ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
         ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
         ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
         ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
         ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

         ppiclf_rprop_gp(PPICLF_JX,ppiclf_npart_gp) = rval(1) ! angularly rotated coordinates of mirror particle
         ppiclf_rprop_gp(PPICLF_JY,ppiclf_npart_gp) = rval(2)
         ppiclf_rprop_gp(PPICLF_JZ,ppiclf_npart_gp) = rval(3)

         ! Store mirror particle properties and Angularly Rotate properties that need to be rotated
         call ppiclf_comm_RotateAngularGhostProperties(ppiclf_npart_gp
     >                                                            ,CW)
         ! Step 5 - check the neighborhood of the mirror ghost particle by looping over the faces, edges, and corners.

         rxl = ppiclf_binb(1) + ppiclf_bins_dx(1)*iip ! min-x of bin where mirror is located
         rxr = rxl + ppiclf_bins_dx(1)                ! max-x of bin where mirror is located
         ryl = ppiclf_binb(3) + ppiclf_bins_dx(2)*jjp ! min-y of bin where mirror is located
         ryr = ryl + ppiclf_bins_dx(2)                ! max-y of bin where mirror is located
         rzl = 0.0d0
         rzr = 0.0d0
         if (ppiclf_ndim .gt. 2) then
            rzl = ppiclf_binb(5) + ppiclf_bins_dx(3)*kkp ! min-z of bin where ip is located
            rzr = rzl + ppiclf_bins_dx(3)                ! max-z of bin where ip is located
         endif

         isave = 0

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) ! i-index of bin face  
            jj1 = jjp + el_face_num(ist+2) ! j-index of bin face
            kk1 = kkp + el_face_num(ist+3) ! k-index of bin face

            iig = ii1
            jjg = jj1
            kkg = kk1

            ! Check if particle is within user-input neighbor-width from bin boundary
            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
              ! calculate x-distance check based on user-input max(filter-width, neighbor-width)
               distchk = distchk + ppiclf_d2chk(1)**2     
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2 ! calculate x-distance of ip particle from x-bin boundary
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2 ! calculate y-distance of ip particle from y-bin boundary
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2 ! calculate z-distance of ip particle from z-bin boundary
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk) ! calculate distance magnitude check based on neighbor-width
            dist = sqrt(dist)    ! calculate distance magnitude of ip particle from bin boundary

            ! if ip particle is farther away from bin boundary than the user-input neighbor-width -> dont create ghost
            if (dist .gt. distchk) cycle

            ! if ghost particle not it bin bounds -> don't create it
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
              cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
              cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
              cycle
            endif

            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn
          
            ! if particle is still in the same processorn -> do not create a ghost
            if (nrank .eq. ppiclf_nid) cycle

            ! gpsave keeps track of which nrank destinations we've already created ghosts for                                     
            ! if a ghost for destination nrank was already created earlier, and this is not a boundary crossing -> don't duplicate
            ! it.                                                                                                                 
            do i=1,isave
               if (gpsave(i) .eq. nrank) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ! Step 6 - Modify ghost properties accordingly
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ! copy the properties of the mirror particle or the previous created ghost particle
            do k=1,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = 
     >                           ppiclf_rprop_gp(k,ppiclf_npart_gp-1)
            enddo
!            write(1920+ppiclf_nid,*) "NewFaces", ppiclf_time, ! 0-1   
!     >              ppiclf_nid, nrank, ppiclf_npart_gp,       ! 2-4
!     >              ppiclf_rprop_gp(1:6, ppiclf_npart_gp),    ! 5-7, 8,10
!     >              ppiclf_y(1:3,ip), ppiclf_y(4:6,ip)        ! 11-13, 14-16


  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
              cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
              cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
              cycle
            endif

            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ! Step 6 - Modify ghost properties accordingly
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ! copy the properties of the mirror particle or the previous created ghost particle
            do k=1,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = 
     >                           ppiclf_rprop_gp(k,ppiclf_npart_gp-1)
            enddo

!            write(1920+ppiclf_nid,*) "NewEdges", ppiclf_time, ! 0-1   
!     >              ppiclf_nid, nrank, ppiclf_npart_gp,       ! 2-4
!     >              ppiclf_rprop_gp(1:6, ppiclf_npart_gp),    ! 5-7, 8,10
!     >              ppiclf_y(1:3,ip), ppiclf_y(4:6,ip)        ! 11-13, 14-16
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + ppiclf_d2chk(1)**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
              cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
              cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
              cycle
            endif

            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank) goto 333
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ! Step 6 - Modify ghost properties accordingly
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            do k=1,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = 
     >                           ppiclf_rprop_gp(k,ppiclf_npart_gp-1)
            enddo

!            write(1920+ppiclf_nid,*) "NewCorners", ppiclf_time, ! 0-1   
!     >              ppiclf_nid, nrank, ppiclf_npart_gp,       ! 2-4
!     >              ppiclf_rprop_gp(1:6, ppiclf_npart_gp),    ! 5-7, 8,10
!     >              ppiclf_y(1:3,ip), ppiclf_y(4:6,ip)        ! 11-13, 14-16
  333 continue
         enddo

      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_CheckLinearPeriodic(distmin,distmax,rxnew)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      real*8 distmin(3), distmax(3)
! 
! Local:
!  
      real*8 distchk, rxdrng(3)
! Input/Output:
!
      real*8 rxnew(3)
!

      distchk = ppiclf_d2chk(1)
      rxdrng(1) = ppiclf_xdrange(2,1) - ppiclf_xdrange(1,1)
      rxdrng(2) = ppiclf_xdrange(2,2) - ppiclf_xdrange(1,2)
      rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)

      ! particle leaving from min x periodic face
      if (distmin(1) .le. distchk) then
         rxnew(1) = rxnew(1) + rxdrng(1)
         goto 123
      endif
      ! particle leaving from max x periodic face
      if (distmax(1) .le. distchk) then
         rxnew(1) = rxnew(1) - rxdrng(1)
         goto 123
      endif

  123 continue    
      
      ! particle leaving from min y periodic face
      if (distmin(2) .le. distchk) then
         rxnew(2) = rxnew(2) + rxdrng(2)
         goto 124
      endif
      ! particle leaving from max y periodic face
      if (distmax(2) .le. distchk) then
         rxnew(2) = rxnew(2) - rxdrng(2)
         goto 124
      endif

  124 continue

      if (ppiclf_ndim .gt. 2) then
        ! particle leaving from min z periodic face
        if (distmin(3) .le. distchk) then
           rxnew(3) = rxnew(3) + rxdrng(3)
           goto 125
        endif
        ! particle leaving from max z periodic face
        if (distmax(3) .le. distchk) then
           rxnew(3) = rxnew(3) - rxdrng(3)
           goto 125
        endif
      endif

  125 continue

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_RotateAngularGhostProperties(i,CW)
!
      implicit none
!
      include "PPICLF"
!
! Input:
      integer*4 i
      logical CW
!
! Internal:
      integer*4 nvec
      parameter nvec = 11
      integer*4 jvec(nvec), k, j
      real*8 rotmat(3,3), rprop(3)
      ! Define all the base indices (x-component) of the 3D vector properties
      ! jvec contains the starting indices (x-components) of 3D vectors to be rotated, in ascending order
      ! To rotate new vectors, add their x-component index here and increase nvec accordingly
      parameter jvec = (/ PPICLF_JVX,
     >           PPICLF_JOX, 
     >           PPICLF_LRS + PPICLF_R_JUX, 
     >           PPICLF_LRS + PPICLF_R_WDOTX,
     >           PPICLF_LRS + PPICLF_R_FQSX,
     >           PPICLF_LRS + PPICLF_R_FAMX,
     >           PPICLF_LRS + PPICLF_R_FAMBX, 
     >           PPICLF_LRS + PPICLF_R_FCX,
     >           PPICLF_LRS + PPICLF_R_FVUX,
     >           PPICLF_LRS + PPICLF_R_FPGX,
     >           PPICLF_LRS + PPICLF_LRP + PPICLF_P_JFX /)
!      
      ! Choose rotation matrix based on where real particle is leaving the domain
      if(CW) then
        rotmat = rotCW
      elseif(.not. CW) then
        rotmat = rotCCW
      else
        call ppiclf_exittr('Error in RotateAngularGhostProperties 
     >    Rotation Matrix Choice', 0.0d0, 0)
      endif

  
      ! Properties currently being rotated:
      !
      ! ppiclf_y: 
      !        velocity, JVX JVY JVZ
      !        angular velocity, JOX JOY JOZ
      ! ppiclf_rprop:
      !        interpolated fluid velocity, JUX JUY JUZ
      !        density weighted particle acceleration, WDOTX
      !        quasi-steady force, FQSX
      !        added-mass force, FAMX
      !        binary added-mass force, FAMBX
      !        collision  force, FCX
      !        viscous unsteady force, FVUX
      !        pressure gradient, FPGCX
      ! ppiclf_cp_map: 
      !        JFX, JFY, JFZ



      j = 1
      ! Loop over all ppiclf_rprop_gp properties and rotate the ones specified in jvec
      do k=4,PPICLF_LRP_GP
        ! rotate 3D vector property if specified in jvec
        if(k == jvec(j)) then
          rprop(1) = ppiclf_cp_map(k  , i)
          rprop(2) = ppiclf_cp_map(k+1, i)
          rprop(3) = ppiclf_cp_map(k+2, i)

          rprop = MATMUL(rotmat, rprop)

          ppiclf_rprop_gp(k  , i) = rprop(1)
          ppiclf_rprop_gp(k+1, i) = rprop(2)
          ppiclf_rprop_gp(k+2, i) = rprop(3)
          j = j +1

        ! copy scalar or vector-component  property from ppiclf_cp_map 
        ! to ppiclf_rprop_gp as-is if not in jvec
        else
          ppiclf_rprop_gp(k,i) = ppiclf_cp_map(k,i)
        endif
      end do

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_MoveGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      logical partl         
!
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >                                  ,ppiclf_npart_gp,PPICLF_LPART_GP
     >                                  ,ppiclf_iprop_gp,PPICLF_LIP_GP
     >                                  ,partl,0
     >                                  ,ppiclf_rprop_gp,PPICLF_LRP_GP
     >                                  ,1)

      return
      end
c----------------------------------------------------------------------
!      subroutine ppiclf_comm_DavidCreateGhost
!!
!      implicit none
!!
!      include "PPICLF"
!!
!! Internal:
!!
!      real*8 xdlen,ydlen,zdlen,rxdrng(3),rxnew(3), rxval, ryval,
!     >       rzval, rxl, ryl, rzl, rxr, ryr, rzr, distchk, dist
!      integer*4 iadd(3),gpsave(27)
!      real*8 map(PPICLF_LRP_PRO)
!      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
!     >           nfacegp, nedgegp, ncornergp, iperiodicx, iperiodicy,
!     >           iperiodicz, jx, jy, jz, ip, idum, iip, jjp, kkp, ii1,
!     >           jj1, kk1, iig, jjg, kkg, iflgx, iflgy, iflgz,
!     >           isave, iflgsum, ndumn, nrank, ibctype, i, ifc, ist, j,
!     >           k
!!
!
!c     face, edge, and corner number, x,y,z are all inline, so stride=3
!      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
!      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
!     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
!     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
!      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
!     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)
!
!      nfacegp   = 4  ! number of faces
!      nedgegp   = 4  ! number of edges
!      ncornergp = 0  ! number of corners
!
!      if (ppiclf_ndim .gt. 2) then
!         nfacegp   = 6  ! number of faces
!         nedgegp   = 12 ! number of edges
!         ncornergp = 8  ! number of corners
!      endif
!
!      iperiodicx = ppiclf_iperiodic(1)
!      iperiodicy = ppiclf_iperiodic(2)
!      iperiodicz = ppiclf_iperiodic(3)
!
!! ------------------------
!c CREATING GHOST PARTICLES
!! ------------------------
!      jx    = 1
!      jy    = 2
!      jz    = 3
!
!      xdlen = ppiclf_binb(2) - ppiclf_binb(1)
!      ydlen = ppiclf_binb(4) - ppiclf_binb(3)
!      
!      zdlen = -1.
!      if (ppiclf_ndim .gt. 2) 
!     >   zdlen = ppiclf_binb(6) - ppiclf_binb(5)
!      if (iperiodicx .ne. 0) xdlen = -1
!      if (iperiodicy .ne. 0) ydlen = -1
!      if (iperiodicz .ne. 0) zdlen = -1
!
!      rxdrng(1) = xdlen
!      rxdrng(2) = ydlen
!      rxdrng(3) = zdlen
!      
!      ppiclf_npart_gp = 0
!
!      do ip=1,ppiclf_npart
!
!         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
!     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))
!
!c        idum = 1
!c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
!c        idum = 2
!c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
!c        idum = 3
!c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
!
!         idum = 0
!         do j=1,PPICLF_LRS
!            idum = idum + 1
!            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip)
!         enddo
!         idum = PPICLF_LRS
!         do j=1,PPICLF_LRP
!            idum = idum + 1
!            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip)
!         enddo
!         idum = PPICLF_LRS+PPICLF_LRP
!         do j=1,PPICLF_LRP_PRO
!            idum = idum + 1
!            ppiclf_cp_map(idum,ip) = map(j)
!         enddo
!
!         rxval = ppiclf_cp_map(1,ip) ! ppiclf_y(PPICLF_JX,ip)
!         ryval = ppiclf_cp_map(2,ip) ! ppiclf_y(PPICLF_JY,ip)
!         rzval = 0.0d0
!         if (ppiclf_ndim .gt. 2) rzval = ppiclf_cp_map(3,ip) ! ppiclf_y(PPICLF_JZ,ip)
!
!         iip    = ppiclf_iprop(8,ip) ! i-index of bin where ip particle is located
!         jjp    = ppiclf_iprop(9,ip) ! j-index of bin where ip particle is located
!         kkp    = ppiclf_iprop(10,ip)! k-index of bin where ip particle is located
!
!         rxl = ppiclf_binb(1) + ppiclf_bins_dx(1)*iip ! min-x of bin where ip is located
!         rxr = rxl + ppiclf_bins_dx(1)                ! max-x of bin where ip is located
!         ryl = ppiclf_binb(3) + ppiclf_bins_dx(2)*jjp ! min-y of bin where ip is located
!         ryr = ryl + ppiclf_bins_dx(2)                ! max-y of bin where ip is located
!         rzl = 0.0d0
!         rzr = 0.0d0
!         if (ppiclf_ndim .gt. 2) then
!            rzl = ppiclf_binb(5) + ppiclf_bins_dx(3)*kkp ! min-z of bin where ip is located
!            rzr = rzl + ppiclf_bins_dx(3)                ! max-z of bin where ip is located
!         endif
!
!         isave = 0
!
!         ! faces
!         do ifc=1,nfacegp
!            ist = (ifc-1)*3
!            ii1 = iip + el_face_num(ist+1) ! i-index of bin face  
!            jj1 = jjp + el_face_num(ist+2) ! j-index of bin face
!            kk1 = kkp + el_face_num(ist+3) ! k-index of bin face
!
!            iig = ii1
!            jjg = jj1
!            kkg = kk1
!
!            ! Step 1 - Check if particle is within user-input neighbor-width from bin boundary
!            distchk = 0.0d0
!            dist = 0.0d0
!            if (ii1-iip .ne. 0) then
!               distchk = distchk + ppiclf_d2chk(1)**2     ! calculate x-distance check based on user-input neighbor-width
!               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2 ! calculate x-distance of ip particle from x-bin boundary
!               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
!            endif
!            if (jj1-jjp .ne. 0) then
!               distchk = distchk + ppiclf_d2chk(1)**2
!               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2 ! calculate y-distance of ip particle from y-bin boundary
!               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
!            endif
!            if (ppiclf_ndim .gt. 2) then
!            if (kk1-kkp .ne. 0) then
!               distchk = distchk + ppiclf_d2chk(1)**2
!               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2 ! calculate z-distance of ip particle from z-bin boundary
!               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
!              endif
!            endif
!            distchk = sqrt(distchk) ! calculate distance magnitude check based on neighbor-width
!            dist = sqrt(dist)    ! calculate distance magnitude of ip particle from bin boundary
!
!            ! if ip particle is farther away from bin boundary than the user-input neighbor-width -> dont create ghost
!            if (dist .gt. distchk) cycle
!
!            iflgx = 0
!            iflgy = 0
!            iflgz = 0
!            
!            ! Step 2 - Check whether particle is about to leave the periodic domain
!
!            ! David previously used to fix the bin boundaries as big as the periodic domain 
!            ! and evaluate the ghost particle using the 3 if statements below. 
!
!            ! This is now commented out as we check for periodicity in a separate subroutine 
!            ! by initializing linear periodic plane and evaluating the distance of the real particle
!            ! relative to the periodic planes. 
!
!            ! David's ghost algorithm is just used for mapping the real particle (creating ghosts for it)
!            ! when needed in the nearby neighboring bins. 
!
!            ! periodic if out of domain - add some ifsss
!            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
!               iflgx = 1
!               iig =modulo(iig,ppiclf_n_bins(1))
!               if (iperiodicx .ne. 0) cycle
!            endif
!            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
!               iflgy = 1
!               jjg =modulo(jjg,ppiclf_n_bins(2))
!               if (iperiodicy .ne. 0) cycle
!            endif
!            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
!               iflgz = 1  
!               kkg =modulo(kkg,ppiclf_n_bins(3))
!               if (iperiodicz .ne. 0) cycle
!            endif
!
!            ! if particle is crossing the periodic domain in each direction -> iflgsum = 1 crossing in 1-direction
!            !                                                     otherwise -> iflgsum = 0
!            iflgsum = iflgx + iflgy + iflgz
!            ndumn = iig + ppiclf_n_bins(1)*jjg 
!     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
!            nrank = ndumn
!
!            ! if ghost rank is the same processor of the real particle .and. not crossing the periodic domain -> do not create a ghost
!            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle
!
!            ! gpsave keeps track of which nrank destinations we've already created ghosts for
!            ! if a ghost for destination nrank was already created earlier, and this is not a boundary crossing -> don't duplicate
!            ! it.
!            do i=1,isave
!               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
!            enddo
!            isave = isave + 1
!            gpsave(isave) = nrank
!
!            ibctype = iflgx+iflgy+iflgz
!                 
!            ! Step 3 - Modify ghost properties accordingly
!            
!            rxnew(1) = rxval
!            rxnew(2) = ryval
!            rxnew(3) = rzval
!       
!            iadd(1) = ii1
!            iadd(2) = jj1
!            iadd(3) = kk1
!
!            ! Linear periodic properties are now checked and evaluated in a separate subroutine
!
!            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
!                 
!            ppiclf_npart_gp = ppiclf_npart_gp + 1
!            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
!            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
!            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
!            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
!            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn
!
!            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
!            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
!            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)
!
!            do k=4,PPICLF_LRP_GP
!               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
!            enddo
!  111 continue
!         enddo
!
!         ! edges
!         do ifc=1,nedgegp
!            ist = (ifc-1)*3
!            ii1 = iip + el_edge_num(ist+1) 
!            jj1 = jjp + el_edge_num(ist+2)
!            kk1 = kkp + el_edge_num(ist+3)
!
!            iig = ii1
!            jjg = jj1
!            kkg = kk1
!
!            distchk = 0.0d0
!            dist = 0.0d0
!            if (ii1-iip .ne. 0) then
!               distchk = distchk + ppiclf_d2chk(1)**2
!               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
!               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
!            endif
!            if (jj1-jjp .ne. 0) then
!               distchk = distchk + ppiclf_d2chk(1)**2
!               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
!               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
!            endif
!            if (ppiclf_ndim .gt. 2) then
!            if (kk1-kkp .ne. 0) then
!               distchk = distchk + ppiclf_d2chk(1)**2
!               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
!               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
!            endif
!            endif
!            distchk = sqrt(distchk)
!            dist = sqrt(dist)
!            if (dist .gt. distchk) cycle
!
!            iflgx = 0
!            iflgy = 0
!            iflgz = 0
!            ! periodic if out of domain - add some ifsss
!            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
!               iflgx = 1
!               iig =modulo(iig,ppiclf_n_bins(1))
!               if (iperiodicx .ne. 0) cycle
!            endif
!            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
!               iflgy = 1
!               jjg =modulo(jjg,ppiclf_n_bins(2))
!               if (iperiodicy .ne. 0) cycle
!            endif
!            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
!               iflgz = 1  
!               kkg =modulo(kkg,ppiclf_n_bins(3))
!               if (iperiodicz .ne. 0) cycle
!            endif
!
!            iflgsum = iflgx + iflgy + iflgz
!            ndumn = iig + ppiclf_n_bins(1)*jjg 
!     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
!            nrank = ndumn
!
!            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle
!
!            do i=1,isave
!               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
!            enddo
!            isave = isave + 1
!            gpsave(isave) = nrank
!
!            ibctype = iflgx+iflgy+iflgz
!                 
!            rxnew(1) = rxval
!            rxnew(2) = ryval
!            rxnew(3) = rzval
!       
!            iadd(1) = ii1
!            iadd(2) = jj1
!            iadd(3) = kk1
!
!            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
!                 
!            ppiclf_npart_gp = ppiclf_npart_gp + 1
!            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
!            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
!            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
!            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
!            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn
!
!            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
!            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
!            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)
!
!            do k=4,PPICLF_LRP_GP
!               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
!            enddo
!  222 continue
!         enddo
!
!         ! corners
!         do ifc=1,ncornergp
!            ist = (ifc-1)*3
!            ii1 = iip + el_corner_num(ist+1) 
!            jj1 = jjp + el_corner_num(ist+2)
!            kk1 = kkp + el_corner_num(ist+3)
!
!            iig = ii1
!            jjg = jj1
!            kkg = kk1
!
!            distchk = 0.0d0
!            dist = 0.0d0
!            if (ii1-iip .ne. 0) then
!               distchk = distchk + ppiclf_d2chk(1)**2
!               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
!               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
!            endif
!            if (jj1-jjp .ne. 0) then
!               distchk = distchk + ppiclf_d2chk(1)**2
!               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
!               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
!            endif
!            if (ppiclf_ndim .gt. 2) then
!            if (kk1-kkp .ne. 0) then
!               distchk = distchk + ppiclf_d2chk(1)**2
!               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
!               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
!            endif
!            endif
!            distchk = sqrt(distchk)
!            dist = sqrt(dist)
!            if (dist .gt. distchk) cycle
!
!            iflgx = 0
!            iflgy = 0
!            iflgz = 0
!            ! periodic if out of domain - add some ifsss
!            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
!               iflgx = 1
!               iig =modulo(iig,ppiclf_n_bins(1))
!               if (iperiodicx .ne. 0) cycle
!            endif
!            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
!               iflgy = 1
!               jjg =modulo(jjg,ppiclf_n_bins(2))
!               if (iperiodicy .ne. 0) cycle
!            endif
!            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
!               iflgz = 1  
!               kkg =modulo(kkg,ppiclf_n_bins(3))
!               if (iperiodicz .ne. 0) cycle
!            endif
!
!            iflgsum = iflgx + iflgy + iflgz
!            ndumn = iig + ppiclf_n_bins(1)*jjg 
!     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
!            nrank = ndumn
!
!            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle
!
!            do i=1,isave
!               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
!            enddo
!            isave = isave + 1
!            gpsave(isave) = nrank
!
!            ibctype = iflgx+iflgy+iflgz
!                 
!            rxnew(1) = rxval
!            rxnew(2) = ryval
!            rxnew(3) = rzval
!       
!            iadd(1) = ii1
!            iadd(2) = jj1
!            iadd(3) = kk1
!
!            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
!                 
!            ppiclf_npart_gp = ppiclf_npart_gp + 1
!            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
!            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
!            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
!            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
!            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn
!
!            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
!            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
!            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)
!
!            do k=4,PPICLF_LRP_GP
!               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
!            enddo
!  333 continue
!         enddo
!
!      enddo
!
!      return
!      end
c-----------------------------------------------------------------------
!      subroutine ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
!!
!      implicit none
!!
!      include "PPICLF"
!!
!! Input:
!!
!      real*8 rxdrng(3)
!      integer*4 iadd(3)
!!
!! Input/Output:
!!
!      real*8 rxnew(3)
!!
!      ! rxdrng(1) = ppiclf_binb(2) - ppiclf_binb(1)
!      ! rxdrng(1) = -1.0  if not periodic in X
!      if (rxdrng(1) .gt. 0 ) then
!      ! particle leaving from max x periodic face
!      if (iadd(1) .ge. ppiclf_n_bins(1)) then
!         rxnew(1) = rxnew(1) - rxdrng(1)
!         goto 123
!      endif
!      endif
!      ! particle leaving from min x periodic face
!      if (rxdrng(1) .gt. 0 ) then
!      if (iadd(1) .lt. 0) then
!         rxnew(1) = rxnew(1) + rxdrng(1)
!         goto 123
!      endif
!      endif
!
!  123 continue    
!      ! rxdrng(2) = ppiclf_xdrange(2,2) - ppiclf_xdrange(1,2)
!      ! rxdrng(2) = -1.0  if not periodic in Y
!      if (rxdrng(2) .gt. 0 ) then
!      ! particle leaving from max y periodic face
!      if (iadd(2) .ge. ppiclf_n_bins(2)) then
!         rxnew(2) = rxnew(2) - rxdrng(2)
!         goto 124
!      endif
!      endif
!      if (rxdrng(2) .gt. 0 ) then
!      ! particle leaving from min y periodic face
!      if (iadd(2) .lt. 0) then
!         rxnew(2) = rxnew(2) + rxdrng(2)
!         goto 124
!      endif
!      endif
!  124 continue
!
!      if (ppiclf_ndim .gt. 2) then
!        ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
!        ! rxdrng(3) = -1.0  if not periodic in Z
!         if (rxdrng(3) .gt. 0 ) then
!      ! particle leaving from max z periodic face
!         if (iadd(3) .ge. ppiclf_n_bins(3)) then
!            rxnew(3) = rxnew(3) - rxdrng(3)
!            goto 125
!         endif
!         endif
!         if (rxdrng(3) .gt. 0 ) then
!      ! particle leaving from min z periodic face
!         if (iadd(3) .lt. 0) then
!            rxnew(3) = rxnew(3) + rxdrng(3)
!            goto 125
!         endif
!         endif
!      endif
!  125 continue
!
!      return
!      end
c----------------------------------------------------------------------
