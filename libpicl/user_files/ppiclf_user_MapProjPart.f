!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for projection map
!
!
! This subroutine maps certain ppiclf terms onto the fluid mesh
!
!    This routine is needed for
!       (i)  feedback_flag = 1 (rocpicl/PICL_TEMP_Runge.F90)
!       (ii) # PROBE is in the *.inp file (libfloflu/WriteProbe.F90)

!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)
!
      implicit none
      include "PPICLF"
!
! Input:
!
      real*8 y    (PPICLF_LRS)
      real*8 ydot (PPICLF_LRS)
      real*8 ydotc(PPICLF_LRS)
      real*8 rprop(PPICLF_LRP)
!
! Output:
!
      real*8 map  (PPICLF_LRP_PRO)
!
! Code:
!
!
      !Add here sploading factor         
      map(PPICLF_P_JPHIP) = rprop(PPICLF_R_JVOLP)*rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JFX)   = ydotc(PPICLF_JVX) * rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JFY)   = ydotc(PPICLF_JVY) * rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JFZ)   = ydotc(PPICLF_JVZ) * rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JE)    = ydotc(PPICLF_JT)  * rprop(PPICLF_R_JSPL)

      ! TLJ - modified 12/21/2024
      map(PPICLF_P_JPHIPD) = rprop(PPICLF_R_JVOLP)*rprop(PPICLF_R_JRHOP)
      map(PPICLF_P_JPHIPU) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JVX)
      map(PPICLF_P_JPHIPV) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JVY)
      map(PPICLF_P_JPHIPW) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JVZ)
      map(PPICLF_P_JPHIPT) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JT)


      return
      end
