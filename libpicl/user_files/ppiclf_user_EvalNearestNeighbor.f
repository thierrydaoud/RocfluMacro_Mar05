!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine to find nearest neighbor for particle-particle and
!     particle-wall interactions. Subroutine also includes
!     the new added-mass binary terms, developed by Sam Briney.
!
! if collisional_flag = 1  F = Fn
!                     = 2  F = Fn + Ft + Tt
!                     = 3  F = Fn + Ft + Tt + Th + Tr
! where Tt = collisional torque
!       Th = hydrodynamic torque
!       Tr = rolling torque
!
! Note that the collisional and rolling torques are due to
!     particle-particle interactions and thus are evaulated
!     in ppiclf_user_EvalNearestNeighbor. Therefore, only the
!     hydrodynamic torque is left to be calculated.
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_EvalNearestNeighbor
     >    (i,j,yi,rpropi,yj,rpropj)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData, ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, write_forces
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData, ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, write_forces

      integer*4 i
      integer*4 j
      real*8 yi    (PPICLF_LRS)    
      real*8 rpropi(PPICLF_LRP)
      real*8 yj    (PPICLF_LRS)    
      real*8 rpropj(PPICLF_LRP)
!
! Internal:
!
      real*8 pi2, rthresh, rxdiff, rydiff, rzdiff, rdiff, rm1, rm2,
     >       rmult, eta_n, rbot, rn_12x, rn_12y, rn_12z, rdelta12,
     >       rv12_mag, rv12_mage, rksp_max, rnmag, rksp_wall, rextra,
     >       JDP_i,JDP_j
      real*8 eps

      ! Sam - from TLJ for box filter
      integer*4 ifilt
      real*8 adptfilter, dpl, phip, dist2, xdist2, ydist2, zdist2
      real*8 dist, rsig
      real*8 sig2, gkern, pi

      ! 06/06/2024 - Thierry - Added Mass code
      integer*4 k, l, kk, ll
      real*8 alpha_local, rad
      real*8 rxdiff1, rydiff1, rzdiff1

      ! 07/16/2024 - TLJ added tangential component
      ! - does not take into account angular velocity
      real*8 unx, uny, unz, utx, uty, utz, ut_mag, rn_mag
      real*8 rt_12x, rt_12y, rt_12z
      real*8 Fn_mag, Ftmin
      real*8 eta_t, mu_c
      real*8 A12x, A12y, A12z 
      real*8 rad1, rad2
      real*8 u12x, u12y, u12z
      real*8 tcx, tcy, tcz
      real*8 trx, try, trz
      real*8 thetar, dp1, dp2, r12
      real*8 omgrx, omgry, omgrz, omgr_mag
      real*8 Ftx, Fty, Ftz

!
! Code:
!
      pi   = acos(-1.0d0)
      pi2  = pi*pi

      ! other particles
      if (j .ne. 0) then
         !Added spload and radius factor

         ! Compute mean particle diameter between i and j
         rthresh  = 0.5d0*(rpropi(PPICLF_R_JDP) + rpropj(PPICLF_R_JDP))

         ! Compute distance between centers of particles i and j
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = rxdiff**2 + rydiff**2 + rzdiff**2
         rdiff = sqrt(rdiff)

!-----------------------------------------------------------------------
!
         ! For binary added-mass for Briney model

         ! 06/06/2024 - Thierry - Added Mass code continues here
         ! 07/09/2024 - TLJ - Updated
         ! 07/14/2024 - Thierry - Updated overlapping particles if statement

         ! Filter widths are input values from *.inp
         ! ppiclf_d2chk(1 & 2) are 1/2 filter width - user defined
         !   = FILTERWIDTH / 2
         ! ppiclf_d2chk(3) is neighbor width - user defined
         !   = max(NEIGHBORWIDTH,4*Dp)
        
         if (am_flag == 2 .and. rdiff <= ppiclf_d2chk(3)) then
            ! Do not overwrite rxdiff, rydiff, rzdiff
            rxdiff1 = rxdiff
            rydiff1 = rydiff
            rzdiff1 = rzdiff

            ! Check if particles are overlapping, replace value
            !   if yes, since we get crazy resistance values
            if (rdiff < rthresh) then
               rxdiff1 = rthresh
               rydiff1 = rthresh
               rzdiff1 = rthresh 
            endif

            ! Model only valid for local volume fraction
            ! less than 0.4, so we limit it here without
            ! over riding rphip
            ! limit alpha to mitigate misuse
            alpha_local = min(0.4, rphip) 
            
            ! Compute the resistance matrix
            ! Only valid for monodispersed particles
            rad = 0.5d0*rpropi(PPICLF_R_JDP)

            call resistance_pair(rxdiff1, rydiff1, rzdiff1, 
     >           alpha_local, rad, R_pair)
            
            ! accumulate number of neighbors
            nneighbors = nneighbors + 1
            
            do k=1,3
               do l=1,3
                 ll = PPICLF_R_WDOTX + (l-1)
                 ! added mass
                 Fam(k) = Fam(k) + R_pair(k,l)   * rpropi(ll)
                 ! induced added mass
                 Fam(k) = Fam(k) + R_pair(k,l+3) * rpropj(ll)
               end do ! l-loop

               ! accumulate neighbor acceleration
               kk = PPICLF_R_WDOTX + (k-1)
               Wdot_neighbor_mean(k) = Wdot_neighbor_mean(k)
     >                               + rpropj(kk)
            end do ! k-loop
         end if ! am_flag==2 .and. rdiff <= ppiclf_d2chk(3)

!-----------------------------------------------------------------------
!
         ! For particle-particle collision

         ! Cycle if rdiff > rthresh + repi
         ! For eps, see Capecelatro etal, JCP, 2013
         eps = 0.075*min(rpropi(PPICLF_R_JDP),rpropj(PPICLF_R_JDP))
         !eps = 0.0d0

         if (rdiff .lt. rthresh+eps) then

            rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
            rm2 = rpropj(PPICLF_R_JRHOP)*rpropj(PPICLF_R_JVOLP)
         
            rmult = 1.0d0/(1.0d0/rm1+1.0d0/rm2)
            eta_n = -2.0d0*sqrt(ksp)*log(erest)/sqrt(log(erest)**2+pi2)
     >              *sqrt(rmult)

            ! Compute unit normal vector along line of contact 
            !   pointing from particle i to particle j
            rbot = 1.0d0/rdiff
            rn_12x = rxdiff*rbot
            rn_12y = rydiff*rbot
            rn_12z = rzdiff*rbot
            rn_mag = rdiff
         
            ! Relative velocity in normal direction
            u12x = yi(PPICLF_JVX)-yj(PPICLF_JVX)
            u12y = yi(PPICLF_JVY)-yj(PPICLF_JVY)
            u12z = yi(PPICLF_JVZ)-yj(PPICLF_JVZ)

            if (collisional_flag>=2) then
               ! Add contribution from angular velocity
               rad1 = 0.5d0*rpropi(PPICLF_R_JDP)
               rad2 = 0.5d0*rpropj(PPICLF_R_JDP)
               A12x = rad1*yi(PPICLF_JOX) + rad2*yj(PPICLF_JOX)
               A12y = rad1*yi(PPICLF_JOY) + rad2*yj(PPICLF_JOY)
               A12z = rad1*yi(PPICLF_JOZ) + rad2*yj(PPICLF_JOZ)

               u12x = u12x + (A12y*rn_12z - A12z*rn_12y)
               u12y = u12y + (A12z*rn_12x - A12x*rn_12z)
               u12z = u12z + (A12x*rn_12y - A12y*rn_12x)
            endif

            ! Compute (u_ij \cdot n_ij)
            rv12_mag = u12x*rn_12x
     >               + u12y*rn_12y
     >               + u12z*rn_12z
         
            ! Compute delta_12 and normal parameters
            rdelta12 = rthresh - rdiff
            rksp_max  = ksp*rdelta12
            rv12_mage = rv12_mag*eta_n
            rnmag     = -rksp_max - rv12_mage

            Fn_mag = abs(rnmag)

            ! Compute tangential unit vector
            unx = rv12_mag*rn_12x
            uny = rv12_mag*rn_12y
            unz = rv12_mag*rn_12z
            utx = u12x - unx
            uty = u12y - uny
            utz = u12z - unz
            ut_mag = sqrt(utx*utx + uty*uty + utz*utz)
            ut_mag = max(ut_mag,1.0d-8)
            rt_12x = utx/ut_mag
            rt_12y = uty/ut_mag
            rt_12z = utz/ut_mag

            ! Compute tangential collision force
            if (ut_mag > 0) then
               mu_c  = 0.4d0  ! Dimensionless; Coulomb
               eta_t = eta_n  ! Set to normal; damping
               Ftmin  = -min(mu_c*Fn_mag,eta_t*ut_mag)  
            endif
            if (collisional_flag==1) then ! Normal component only
               Ftmin = 0.0d0
            endif

            ! Compute contributions to angular velocities
            tcx = 0.0d0; tcy = 0.0d0; tcz = 0.0d0;
            trx = 0.0d0; try = 0.0d0; trz = 0.0d0;
            if (collisional_flag>=2) then

               ! Collision torque contribution
               Ftx = Ftmin*rt_12x
               Fty = Ftmin*rt_12y
               Ftz = Ftmin*rt_12z
               rad1 = 0.5d0*rpropi(PPICLF_R_JDP)
               tcx = rad1*(rn_12y*Ftz - rn_12z*Fty)
               tcy = rad1*(rn_12z*Ftx - rn_12x*Ftz)
               tcz = rad1*(rn_12x*Fty - rn_12y*Ftx)

               if (collisional_flag==3) then
                  ! Rolling torque contribution
                  thetar = 0.06  ! Needs to be calibrated
                  dp1 = rpropi(PPICLF_R_JDP)
                  dp2 = rpropj(PPICLF_R_JDP)
                  r12 = 0.5d0*(dp1*dp2)/(dp1+dp2)
                  omgrx = yi(PPICLF_JOX) - yj(PPICLF_JOX)
                  omgry = yi(PPICLF_JOY) - yj(PPICLF_JOY)
                  omgrz = yi(PPICLF_JOZ) - yj(PPICLF_JOZ)
                  omgr_mag = sqrt(omgrx*omgrx+omgry*omgry+omgrz*omgrz)
                  omgr_mag = max(omgr_mag,1.d-8)
                  trx = -thetar*Fn_mag*r12*omgrx/omgr_mag
                  try = -thetar*Fn_mag*r12*omgry/omgr_mag
                  trz = -thetar*Fn_mag*r12*omgrz/omgr_mag
               endif
            endif


            ! Now update that part of the RHS of equations 
            !   that involve nearest neighbors

            ! Particle velocities
            ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
     >                                 + rnmag*rn_12x
     >                                 + Ftmin*rt_12x
            ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
     >                                 + rnmag*rn_12y
     >                                 + Ftmin*rt_12y
            ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
     >                                 + rnmag*rn_12z
     >                                 + Ftmin*rt_12z

            ! Particle angular velocities
            ppiclf_ydotc(PPICLF_JOX,i) = ppiclf_ydotc(PPICLF_JOX,i)
     >                                 + tcx + trx
            ppiclf_ydotc(PPICLF_JOY,i) = ppiclf_ydotc(PPICLF_JOY,i)
     >                                 + tcy + try
            ppiclf_ydotc(PPICLF_JOZ,i) = ppiclf_ydotc(PPICLF_JOZ,i)
     >                                 + tcz + trz

         end if ! rdiff lt rthresh + eps

!-----------------------------------------------------------------------
!
         ! Feedback fluctuation mean

         dist2 = ppiclf_d2chk(2)

         ! Box filter half-width dist2
         if (qs_fluct_filter_adapt_flag==0) then
            dist2 = ppiclf_d2chk(2)
         else if (qs_fluct_filter_adapt_flag>=1) then
            ! Adaptive filter defined wrt particle i
            ! Used for adaptive box or gaussian
            dpl = rpropi(PPICLF_R_JDP)
            phip = rpropi(PPICLF_R_JPHIP)
            adptfilter = ( 10.*(dpl**3)/max(1.e-4,phip) )**(1./3.)
            adptfilter = adptfilter/2.0
            dist2 = max(ppiclf_d2chk(2),adptfilter)
         endif

         ! Check if particle lies inside box or gaussian filter
         xdist2 = abs(yi(PPICLF_JX)-yj(PPICLF_JX))
         if (xdist2 .gt. dist2) return

         ydist2 = abs(yi(PPICLF_JY)-yj(PPICLF_JY))
         if (ydist2 .gt. dist2) return

         if (ppiclf_ndim .eq. 3) then
           zdist2 = abs(yi(PPICLF_JZ)-yj(PPICLF_JZ))
           if (zdist2 .gt. dist2) return
         endif

         !
         ! The mean is calcuated according to Lattanzi etal,
         !   Physical Review Fluids, 2022.
         !
         if (j.ne.0) then
         if (qs_fluct_filter_flag==0) then
           upmean   = upmean + yj(PPICLF_JVX)
           vpmean   = vpmean + yj(PPICLF_JVY)
           wpmean   = wpmean + yj(PPICLF_JVZ)
           u2pmean  = u2pmean + yj(PPICLF_JVX)**2
           v2pmean  = v2pmean + yj(PPICLF_JVY)**2
           w2pmean  = w2pmean + yj(PPICLF_JVZ)**2
           icpmean  = icpmean + 1
         else if (qs_fluct_filter_flag==1) then
           ! See https://dpzwick.github.io/ppiclF-doc/algorithms/overlap_mesh.html
           dist = sqrt(xdist2**2 + ydist2**2 + zdist2**2)
           gkern = sqrt(pi*ppiclf_filter**2/
     >              (4.0d0*log(2.0d0)))**(-ppiclf_ndim) * 
     >              exp(-dist**2/(ppiclf_filter**2/(4.0d0*log(2.0d0))))

           phipmean = phipmean + gkern*rpropj(PPICLF_R_JVOLP)
           upmean   = upmean +
     >                gkern*yj(PPICLF_JVX)*rpropj(PPICLF_R_JVOLP)
           vpmean   = vpmean +
     >                gkern*yj(PPICLF_JVY)*rpropj(PPICLF_R_JVOLP)
           wpmean   = wpmean +
     >                gkern*yj(PPICLF_JVZ)*rpropj(PPICLF_R_JVOLP)
           u2pmean  = u2pmean +
     >               gkern*(yj(PPICLF_JVX)**2)*rpropj(PPICLF_R_JVOLP)
           v2pmean  = v2pmean +
     >               gkern*(yj(PPICLF_JVY)**2)*rpropj(PPICLF_R_JVOLP)
           w2pmean  = w2pmean +
     >               gkern*(yj(PPICLF_JVZ)**2)*rpropj(PPICLF_R_JVOLP)
           icpmean = icpmean + 1
         end if
         end if


!-----------------------------------------------------------------------
!
      ! boundaries
      elseif (j .eq. 0) then

         rksp_wall = ksp
         !rksp_wall = 1000

         ! give a bit larger collision threshold for walls
         rextra   = 0.05d0 !
         ! add sploading and radius factor 
         rthresh  = (0.5d0+rextra)*rpropi(PPICLF_R_JDP)
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = rxdiff**2 + rydiff**2 + rzdiff**2
         rdiff = sqrt(rdiff)
         
         if (rdiff .gt. rthresh) return

         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         
         rmult = sqrt(rm1)
         eta_n = 2.0d0*sqrt(rksp_wall)*log(erest)
     >           /sqrt(log(erest)**2+pi2)*rmult
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = -yi(PPICLF_JVX)*rn_12x
     >              -yi(PPICLF_JVY)*rn_12y
     >              -yi(PPICLF_JVZ)*rn_12z

         rv12_mage = rv12_mag*eta_n
         rksp_max  = rksp_wall*rdelta12
         rnmag     = -rksp_max - rv12_mage

         
         ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
     >                              + rnmag*rn_12y
         ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
     >                              + rnmag*rn_12z
        
       !write(*,*) "Wall NEAR",i,ppiclf_ydotc(PPICLF_JVY,i)  
      endif


      return
      end
