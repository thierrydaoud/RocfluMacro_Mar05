!-----------------------------------------------------------------------
!
! Created Spet. 14, 2024
!
! Subroutine for unit test problems
!      
! stationary = -1: azimuthal velocity only
!            = -2: radial velocity only
!            = -3: radial + azimuthal velocity
!            = -4: x-velocity only
!            = -5: y-velocity only
!
! NOTE: Time step ppiclf_dt is hardcoded in ppiclf_solve_IntegrateRK3s_Rocflu
!       subroutine to be ppiclf_dt = 5.0000000000000004E-008      
!-----------------------------------------------------------------------
      subroutine ppiclf_user_unit_tests(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, flow_model
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, flow_model

      integer*4 i, iStage
      real*8 famx, famy, famz
!
! Code:
!
      if (stationary==-1) then
         call ppiclf_user_unit01(i,iStage,famx,famy,famz)
      elseif (stationary==-2) then
         call ppiclf_user_unit02(i,iStage,famx,famy,famz)
      elseif (stationary==-3) then
         call ppiclf_user_unit03(i,iStage,famx,famy,famz)
      elseif (stationary==-4) then
         call ppiclf_user_unit04(i,iStage,famx,famy,famz)
      elseif (stationary==-5) then
         call ppiclf_user_unit05(i,iStage,famx,famy,famz)
      endif

      return
      end
!      
!-----------------------------------------------------------------------
!
! Azimuthal(theta) velocity only
!
      subroutine ppiclf_user_unit01(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
      real*8 pi, omg00, omg02, tau

!
! Code:
!
      pi = acos(-1.0d0)
      tau = 0.1d-3 ! time for 2pi revolution. Can be max simulation time
      omg00 = (2.0d0*pi)/tau
      omg02 = omg00*omg00

      if(ppiclf_time .eq. 0.0) then
        ppiclf_y(PPICLF_JVX,i) = -omg00*ppiclf_y(PPICLF_JY,i)
        ppiclf_y(PPICLF_JVY,i) =  omg00*ppiclf_y(PPICLF_JX,i)
        ppiclf_y(PPICLF_JVZ,i) =  0.0d0
      endif
      
      famx = -omg02*ppiclf_y(PPICLF_JX,i)
      famy = -omg02*ppiclf_y(PPICLF_JY,i)
      famz = 0.0d0

      return
      end
!
!-----------------------------------------------------------------------
!
! Radial velocity only
!
      subroutine ppiclf_user_unit02(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
      real*8 drdt00, drdt02, theta
!
! Code:
!
      ! dx/dt = dr/dt * cos theta
      ! dy/dt = dr/dt * sin theta
      drdt00 = 79.928163327808903
      drdt02 = drdt00*drdt00
      
      ! angle between particle and +ve x-axis
      theta = atan2(ppiclf_y(PPICLF_JY,i), ppiclf_y(PPICLF_JX,i))

      if(ppiclf_time .eq. 0.0) then
        ppiclf_y(PPICLF_JVX,i) =  -drdt00*cos(theta)
        ppiclf_y(PPICLF_JVY,i) =  -drdt00*sin(theta)
        ppiclf_y(PPICLF_JVZ,i) =  0.0d0
      endif
      
      famx = -drdt02*cos(theta)
      famy = -drdt02*sin(theta)
      famz = 0.0d0


      return
      end
!
!-----------------------------------------------------------------------
!
! Radial + Azimuthal velocity
!
      subroutine ppiclf_user_unit03(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
      real*8 pi, tau, omg00, omg02, drdt00, theta
!
! Code:
!
      ! dx/dt = dr/dt * cos theta
      ! dy/dt = dr/dt * sin theta
      
      pi = acos(-1.0d0)

      tau = 0.1d-3
      omg00 = (2.0d0*pi)/tau
      omg02 = omg00*omg00
      drdt00 = 79.928163327808903
      
      ! angle between particle and +ve x-axis
      theta = atan2(ppiclf_y(PPICLF_JY,i), ppiclf_y(PPICLF_JX,i))

      if(ppiclf_time .eq. 0.0) then
        ppiclf_y(PPICLF_JVX,i) = -omg00*ppiclf_y(PPICLF_JY,i)
     >                           -drdt00*cos(theta)
        ppiclf_y(PPICLF_JVY,i) =  omg00*ppiclf_y(PPICLF_JX,i) 
     >                           -drdt00*sin(theta)
        ppiclf_y(PPICLF_JVZ,i) =  0.0d0
      endif
      
      famx = 2.0d0*drdt00*omg00*sin(theta) 
     >       -omg02*ppiclf_y(PPICLF_JX,i)
      famy = -2.0d0*drdt00*omg00*cos(theta) 
     >       -omg02*ppiclf_y(PPICLF_JY,i)
      famz = 0.0d0

      return
      end
!
!-----------------------------------------------------------------------
!
! X-velocity only
!
      subroutine ppiclf_user_unit04(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
      real*8 drdt00, drdt02, theta
!
! Code:
!
      ! dx/dt = dr/dt
      ! dy/dt = 0.0
      drdt00 = 799.28163327808903
      drdt02 = drdt00*drdt00
      
      if(ppiclf_time .eq. 0.0) then
        ppiclf_y(PPICLF_JVX,i) =  -drdt00
        ppiclf_y(PPICLF_JVY,i) =  0.0d0
        ppiclf_y(PPICLF_JVZ,i) =  0.0d0
      endif
      
      famx = -drdt02
      famy = 0.0d0
      famz = 0.0d0


      return
      end
!
!-----------------------------------------------------------------------
!
! Y-velocity only
!
      subroutine ppiclf_user_unit05(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
      real*8 drdt00, drdt02, theta
!
! Code:
!
      ! dx/dt = 0.0d0
      ! dy/dt = dr/dt
      drdt00 = 2.0d0*79.928163327808903
      drdt02 = drdt00*drdt00
      
      if(ppiclf_time .eq. 0.0) then
        ppiclf_y(PPICLF_JVX,i) =  0.0d0
        ppiclf_y(PPICLF_JVY,i) =  -drdt00
        ppiclf_y(PPICLF_JVZ,i) =  0.0d0
      endif
      
      famx = 0.0d0
      famy = -drdt02
      famz = 0.0d0


      return
      end
!
!-----------------------------------------------------------------------
