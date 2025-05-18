!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
!******************************************************************************
!
! Purpose: 
!
! Description: none.
!
! Input: 
!
! Output:
!
! Notes: 
!
!******************************************************************************
!
! $Id: PICL_F90,v 1.0 2022/05/08 bdurant Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PICL_TEMP_Runge( pRegion)

!  USE 

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_level,t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt

  USE RFLU_ModDifferentiationCells
  USE RFLU_ModLimiters, ONLY: RFLU_CreateLimiter, &
                              RFLU_ComputeLimiterBarthJesp, &
                              RFLU_ComputeLimiterVenkat, &
                              RFLU_LimitGradCells, &
                              RFLU_LimitGradCellsSimple, &
                              RFLU_DestroyLimiter
  USE RFLU_ModWENO, ONLY: RFLU_WENOGradCellsWrapper, &
                          RFLU_WENOGradCellsXYZWrapper
#ifdef PICL
USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                             RFLU_ConvertCvPrim2Cons

 USE ModInterfaces, ONLY: RFLU_DecideWrite !BRAD added for picl
 
#endif



#ifdef PICL
!DEC$ NOFREEFORM
#include "../libpicl/ppiclF/source/PPICLF_USER.h"
#include "../libpicl/ppiclF/source/PPICLF_STD.h"
!DEC$ FREEFORM
#endif


  IMPLICIT NONE


! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString


TYPE(t_global), POINTER :: global
TYPE(t_level), POINTER :: levels(:)
TYPE(t_region), POINTER :: pRegion
TYPE(t_grid), POINTER :: pGrid
!INTEGER :: errorFlag

#ifdef PICL
  LOGICAL :: doWrite      
  INTEGER(KIND=4) :: i,piclIO,nCells,lx,ly,lz
  INTEGER :: errorFlag,icg      
  REAL(KIND=8) :: piclDtMin,piclCurrentTime, &
          temp_dudtMixt,temp_dvdtMixt,temp_dwdtMixt,energydotg
  REAL(KIND=8) :: dudx,dudy,dudz
  REAL(KIND=8) :: dvdx,dvdy,dvdz
  REAL(KIND=8) :: dwdx,dwdy,dwdz
  REAL(KIND=8) :: vFrac

  REAL(KIND=8), DIMENSION(3) :: ug      
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: rhoF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: uxF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: uyF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: uzF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: csF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: tpF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: ppF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: vfP
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: dpxF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: dpyF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: dpzF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: SDRX
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: SDRY
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: SDRZ
  REAL(KIND=8), DIMENSION(:,:,:), POINTER :: pGc 
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: rhsR        
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: pGcX 
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: pGcY
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: pGcZ
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: JFX
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFXCell
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: JFY
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFYCell
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: JFZ
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFZCell
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: JFE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFECell
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: PhiP
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: YTEMP
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: domgdx
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: domgdy
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: domgdz
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: drhodx
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: drhody
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: drhodz
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: dpvxF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: dpvyF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: dpvzF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: SDOX
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: SDOY
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: SDOZ

  ! TLJ - added for Feedback term - 04/01/2025
  INTEGER, DIMENSION(:), ALLOCATABLE :: varInfoPicl
  INTEGER, DIMENSION(:), POINTER :: piclcvInfo
  REAL(KIND=8) :: dodx, dody, dodz,     &
                  omgx, omgy, omgz,     &
                  divu,                 &
                  dprdx, dprdy, dprdz,  &
                  dpdx, dpdy, dpdz,     &
                  phirho, ir, ir2 ,     &
                  dfxdx, dfxdy, dfxdz,  &
                  dfydx, dfydy, dfydz,  &
                  dfzdx, dfzdy, dfzdz   

  !REAL(KIND=8) :: ppiclf

#endif


   
!******************************************************************************

  RCSIdentString = '$RCSfile: PICL_TEMP_Runge.F90,v $ $Revision: 1.0 $'
 
  global => pRegion%global
  
  CALL RegisterFunction( global, 'PICL_TEMP_Runge',__FILE__ )



! Set pointers ----------------------------------------------------------------

    !pRegion => regions!pLevel%regions(iReg)
    pGrid   => pRegion%grid

!PPICLF Integration
#ifdef PICL

     piclIO = 100000000
     piclDtMin = REAL(global%dtMin,8)
     piclCurrentTime = REAL(global%currentTime,8)

     ! TLJ - 11/23/2024
     !     - This has now been removed
     doWrite = RFLU_DecideWrite(global)
     !Figure out piclIO call, might need to look into timestepping
     IF ( (doWrite .EQV. .TRUE.)) piclIO = 1


!PARTICLE stuff possbile needed
!    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)


!allocate arrays to send to picl
    nCells = pRegion%grid%nCells
    ALLOCATE(rhoF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uxF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uyF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uzF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(csF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(tpF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error    

    ALLOCATE(ppF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error    

    ALLOCATE(vfP(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpxF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
    
    ALLOCATE(dpyF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpzF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRX(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRY(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRZ(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(rhsR(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcX(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcY(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcZ(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFX(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFXCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFY(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFYCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFZ(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFZCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFECell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFE(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(PhiP(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    IF (pRegion%mixtInput%axiFlag) THEN
      ALLOCATE(YTEMP(2,2,2,nCells),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
      END IF ! global%error
    ENDIF

    ALLOCATE(domgdx(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
    
    ALLOCATE(domgdy(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(domgdz(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(drhodx(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(drhody(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(drhodz(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpvxF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpvyF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpvzF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDOX(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDOY(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDOZ(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error



!Might need to update prim like plag does
pGc => pRegion%mixt%gradCell

    ! 04/01/2025 - TLJ - we need feedback terms and their gradients to
    !       calculate the undisturbed torque component
    ! Internal definitions; some redundancy but just ignore
    ! We do not need energy, but might in the future
    DO i = 1,pRegion%grid%nCells
       JFXCell(i) = 0.0_RFREAL
       JFYCell(i) = 0.0_RFREAL
       JFZCell(i) = 0.0_RFREAL
       do lz=1,2
       do ly=1,2
       do lx=1,2 
          call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JFX,JFX(lx,ly,lz,i))  
          call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JFY,JFY(lx,ly,lz,i))
          call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JFZ,JFZ(lx,ly,lz,i))
          JFXCell(i) = JFXCell(i) + JFX(lx,ly,lz,i)
          JFYCell(i) = JFYCell(i) + JFY(lx,ly,lz,i) 
          JFZCell(i) = JFZCell(i) + JFZ(lx,ly,lz,i) 
       end do
       end do
       end do 
       !! Do not multiply by cell volume like what is done for
       !! the right hand side of the Euler/NS equations
       JFXCell(i) = JFXCell(i) * 0.125 * pregion%grid%vol(i)
       JFYCell(i) = JFYCell(i) * 0.125 * pregion%grid%vol(i)
       JFZCell(i) = JFZCell(i) * 0.125 * pregion%grid%vol(i)
       pregion%mixt%piclFeedback(1,i) = JFXCell(i)
       pregion%mixt%piclFeedback(2,i) = JFYCell(i)
       pregion%mixt%piclFeedback(3,i) = JFZCell(i)
    ENDDO
    ! Now calculate the gradient of the feedback force
    ALLOCATE(varInfoPicl(3),STAT=errorFlag)
    ALLOCATE(piclcvInfo(3),STAT=errorFlag)
    varInfoPicl(1) = 1
    varInfoPicl(2) = 2
    varInfoPicl(3) = 3
    piclcvInfo = varInfoPicl
    CALL RFLU_ComputeGradCellsWrapper(pRegion,1,3,1,3,varInfoPicl, &
                                      pRegion%mixt%piclFeedback,&
                                      pRegion%mixt%piclgradFeedback)
    CALL RFLU_WENOGradCellsXYZWrapper(pRegion,1,3, &
                                      pRegion%mixt%piclgradFeedback)
    CALL RFLU_LimitGradCellsSimple(pRegion,1,3,1,3, &
                                   pRegion%mixt%piclFeedback,&
                                   piclcvInfo,&
                                   pRegion%mixt%piclgradFeedback)
    DEALLOCATE(varInfoPicl,STAT=errorFlag)
    DEALLOCATE(piclcvInfo,STAT=errorFlag)
    ! END - TLJ calculating gradient of feedback force

!Fill arrays for interp field
    DO i = 1,pRegion%grid%nCells
!Zero out phip
       PhiP(i) = 0.0_RFREAL

       ug(XCOORD) = pRegion%mixt%cv(CV_MIXT_XMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(YCOORD) = pRegion%mixt%cv(CV_MIXT_YMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(ZCOORD) = pRegion%mixt%cv(CV_MIXT_ZMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ! 03/11/2025 - Thierry - Du/Dt, Dv/Dt, Dw/Dt (not weighted by phi^g or rho^g)

       temp_dudtMixt  = (-pRegion%mixt%rhs(CV_MIXT_XMOM,i)/pRegion%grid%vol(i)& 
                         +ug(XCOORD)*pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i))&
                         /pRegion%mixt%cv(CV_MIXT_DENS,i)&
                         +DOT_PRODUCT(ug,pGc(:,2,i))

       temp_dvdtMixt  = (-pRegion%mixt%rhs(CV_MIXT_YMOM,i)/pRegion%grid%vol(i)& 
                         +ug(YCOORD)*pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i))&
                         /pRegion%mixt%cv(CV_MIXT_DENS,i)&
                         +DOT_PRODUCT(ug,pGc(:,3,i))
                         
       temp_dwdtMixt  = (-pRegion%mixt%rhs(CV_MIXT_ZMOM,i)/pRegion%grid%vol(i)& 
                         +ug(ZCOORD)*pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i))&
                         /pRegion%mixt%cv(CV_MIXT_DENS,i)&
                         +DOT_PRODUCT(ug,pGc(:,4,i))

       do lz=1,2
       do ly=1,2
       do lx=1,2 
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JPHIP,vfP(lx,ly,lz,i))
       PhiP(i) = PhiP(i) +  (0.125*vfP(lx,ly,lz,i))*pRegion%grid%vol(i)

       ! TLJ - 02/07/2025 scaled conserved density by gas-phase volume fraction
       vFrac = 1.0_RFREAL - pRegion%mixt%piclVF(i)
       rhoF(lx,ly,lz,i) = pRegion%mixt%cv(CV_MIXT_DENS,i) / vFrac
       uxF(lx,ly,lz,i) = pRegion%mixt%cv(CV_MIXT_XMOM,i) &
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)
       uyF(lx,ly,lz,i) = pRegion%mixt%cv(CV_MIXT_YMOM,i) &
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)
       uzF(lx,ly,lz,i) = pRegion%mixt%cv(CV_MIXT_ZMOM,i) &
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       csF(lx,ly,lz,i) = pRegion%mixt%dv(DV_MIXT_SOUN,i)
       tpF(lx,ly,lz,i) = pRegion%mixt%dv(DV_MIXT_TEMP,i) 
       ! Davin - added pressure to interpolation values 02/22/2025
       ppF(lx,ly,lz,i) = pRegion%mixt%dv(DV_MIXT_PRES,i) 

       dpxF(lx,ly,lz,i) = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_PRES,i) ! dp/dx
       dpyF(lx,ly,lz,i) = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_PRES,i) ! dp/dy
       dpzF(lx,ly,lz,i) = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_PRES,i) ! dp/dz

       dudx = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_XVEL,i)
       dudy = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_XVEL,i)
       dudz = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_XVEL,i)

       dvdx = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_YVEL,i)
       dvdy = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_YVEL,i)
       dvdz = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_YVEL,i)

       dwdx = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_ZVEL,i)
       dwdy = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_ZVEL,i)
       dwdz = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_ZVEL,i)

       domgdx(lx,ly,lz,i) = dwdy - dvdz
       domgdy(lx,ly,lz,i) = dudz - dwdx
       domgdz(lx,ly,lz,i) = dvdx - dudy

       ! 04/01/2025 - TLJ - Calculate the substantial derivative of vorticity
       ! Internal definitions; some redundancy but just ignore
       dodx   = 0.0_RFREAL ! D(Omega_x)/DT
       dody   = 0.0_RFREAL ! D(Omega_y)/DT
       dodz   = 0.0_RFREAL ! D(Omega_z)/DT
       omgx   = dwdy - dvdz ! Omega_x
       omgy   = dudz - dwdx ! Omega_y
       omgz   = dvdx - dudy ! Omega_z
       divu   = dudx + dvdy + dwdz ! u_x+v_y+w_z; divergence of velocity
       dprdx  = pGc(XCOORD,1,i) ! d(rho phi)/dx
       dprdy  = pGc(YCOORD,1,i) ! d(rho phi)/dy
       dprdz  = pGc(ZCOORD,1,i) ! d(rho phi)/dz
       dpdx   = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_PRES,i) ! dp/dx
       dpdy   = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_PRES,i) ! dp/dy
       dpdz   = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_PRES,i) ! dp/dz
       dfxdx  = pRegion%mixt%piclgradFeedback(XCOORD,1,i) ! dFx/dx
       dfxdy  = pRegion%mixt%piclgradFeedback(YCOORD,1,i) ! dFx/dy
       dfxdz  = pRegion%mixt%piclgradFeedback(ZCOORD,1,i) ! dFx/dz
       dfydx  = pRegion%mixt%piclgradFeedback(XCOORD,2,i) ! dFy/dx
       dfydy  = pRegion%mixt%piclgradFeedback(YCOORD,2,i) ! dFy/dy
       dfydz  = pRegion%mixt%piclgradFeedback(ZCOORD,2,i) ! dFy/dz
       dfzdx  = pRegion%mixt%piclgradFeedback(XCOORD,3,i) ! dFz/dx
       dfzdy  = pRegion%mixt%piclgradFeedback(YCOORD,3,i) ! dFz/dy
       dfzdz  = pRegion%mixt%piclgradFeedback(ZCOORD,3,i) ! dFz/dz
       phirho = pRegion%mixt%cv(CV_MIXT_DENS,i) ! phi_g*rho_g
       ir     = 1.0_RFREAL / phirho
       ir2    = ir*ir
       ! 1. Vortex stretching
       dodx = omgx*dudx + omgy*dudy + omgz*dudz
       dody = omgx*dvdx + omgy*dvdy + omgz*dvdz
       dodz = omgx*dwdx + omgy*dwdy + omgz*dwdz
       ! 2. Vortex dilatation
       dodx = dodx - omgx*divu
       dody = dody - omgy*divu
       dodz = dodz - omgz*divu
       ! 3. Baroclinic
       dodx = dodx + (dprdy*dpdz - dprdz*dpdy)*ir2
       dody = dody + (dprdz*dpdx - dprdx*dpdz)*ir2
       dodz = dodz + (dprdx*dpdy - dprdy*dpdx)*ir2
       ! 4. Torque due to feedback force
       dodx = dodx + (dfzdy - dfydz)*ir
       dody = dody + (dfxdz - dfzdx)*ir
       dodz = dodz + (dfydx - dfxdy)*ir
       ! 5. Misalignment of phi*rho and feedback force
       dodx = dodx + (dprdy*JFZCell(i) - dprdz*JFYCell(i))*ir2
       dody = dody + (dprdz*JFXCell(i) - dprdx*JFZCell(i))*ir2
       dodz = dodz + (dprdx*JFYCell(i) - dprdy*JFXCell(i))*ir2
       ! 6. Add terms and store
       SDOX(lx,ly,lz,i) = dodx
       SDOY(lx,ly,lz,i) = dody
       SDOZ(lx,ly,lz,i) = dodz
       ! End - TLJ - Calculate the substantial derivative of vorticity

       ! Substantial derivative of gas-phase velocity
       SDRX(lx,ly,lz,i) = temp_dudtMixt ! Du/Dt
       SDRY(lx,ly,lz,i) = temp_dvdtMixt ! Dv/Dt
       SDRZ(lx,ly,lz,i) = temp_dwdtMixt ! Dw/Dt

       rhsR(lx,ly,lz,i) = -pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i) ! \p(rho*phi)/\p(t)

       pGcX(lx,ly,lz,i) = pGc(XCOORD,1,i) ! d(rho phi)/dx
       pGcY(lx,ly,lz,i) = pGc(YCOORD,1,i) ! d(rho phi)/dy
       pGcz(lx,ly,lz,i) = pGc(ZCOORD,1,i) ! d(rho phi)/dz

       ! Gradient of rho^g of mixture (not weighted by phi^g!)
       ! Using grad(rhog) directly
       drhodx(lx,ly,lz,i) = pRegion%mixt%piclgradRhog(1,1,i) ! d(rho)/dx
       drhody(lx,ly,lz,i) = pRegion%mixt%piclgradRhog(2,1,i) ! d(rho)/dy
       drhodz(lx,ly,lz,i) = pRegion%mixt%piclgradRhog(3,1,i) ! d(rho)/dz

       ! Viscous term of pressure gradient (divergence of tau)
       dpvxF(lx,ly,lz,i) = pRegion%mixt%diss(CV_MIXT_XMOM,i)/pRegion%grid%vol(i)
       dpvyF(lx,ly,lz,i) = pRegion%mixt%diss(CV_MIXT_YMOM,i)/pRegion%grid%vol(i)
       dpvzF(lx,ly,lz,i) = pRegion%mixt%diss(CV_MIXT_ZMOM,i)/pRegion%grid%vol(i)

       end do
       end do
       end do 
       
       !Dump back VolFrac
       !VOL Frac cap
       PhiP(i) = PhiP(i) / (pRegion%grid%vol(i))
       if (PhiP(i) .gt. 0.62) PhiP(i) = 0.62
       do lz=1,2
       do ly=1,2
       do lx=1,2 
          vfp(lx,ly,lz,i) = PhiP(i)      
       end do
       end do
       end do   

    END DO

! Interp field calls
! TLJ - interpolates various fluid quantities onto the 
!       the ppiclf particle locations
! TLJ PPICLF_LRP_INT in PPICLF_USER.h must match the number
!     of calls to ppiclf_solve_InterpFieldUser
! Davin - added pressure 02/22/2025
      IF (PPICLF_LRP_INT .NE. 30) THEN
         write(*,*) "Error: PPICLF_LRP_INT must be set to 30"
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:LRP_INT')
      endif

      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOF,rhoF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUX,uxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUY,uyF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUZ,uzF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDX,dpxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDY,dpyF)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDZ,dpzF)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JCS,csF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JT,tpF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPHIP,vfP)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRX,SDRX)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRY,SDRY)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRZ,SDRZ)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHSR,rhsR)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCX,pGcX) 
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCY,pGcY) 
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCZ,pGcZ) 
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JXVOR,domgdx)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JYVOR,domgdy)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JZVOR,domgdz)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JP,ppF)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOGX,drhodx)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOGY,drhody)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOGZ,drhodz)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPVDX,dpvxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPVDY,dpvyF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPVDZ,dpvzF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDOX,SDOX)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDOY,SDOY)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDOZ,SDOZ)  


!FEED BACK TERM
!Fill arrays for interp field
IF (global%piclFeedbackFlag == 1) THEN
    DO i = 1,pRegion%grid%nCells
       ug(XCOORD) = pRegion%mixt%cv(CV_MIXT_XMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(YCOORD) = pRegion%mixt%cv(CV_MIXT_YMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(ZCOORD) = pRegion%mixt%cv(CV_MIXT_ZMOM,i)&
                          /pRegion%mixt%cv(CV_MIXT_DENS,i)

       JFXCell(i) = 0.0_RFREAL
       JFYCell(i) = 0.0_RFREAL
       JFZCell(i) = 0.0_RFREAL
       JFECell(i) = 0.0_RFREAL
       do lz=1,2
       do ly=1,2
       do lx=1,2 
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JFX,JFX(lx,ly,lz,i))  
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JFY,JFY(lx,ly,lz,i))
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JFZ,JFZ(lx,ly,lz,i))
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JE,JFE(lx,ly,lz,i))      
       !call get energy  
       JFXCell(i) = JFXCell(i) + JFX(lx,ly,lz,i) ! / pRegion%grid%vol(i)    
       JFYCell(i) = JFYCell(i) + JFY(lx,ly,lz,i) 
       JFZCell(i) = JFZCell(i) + JFZ(lx,ly,lz,i) 
       JFECell(i) = JFECell(i) + JFE(lx,ly,lz,i)  
       !Jenergy = +...
       end do
       end do
       end do 
       JFXCell(i) = JFXCell(i) * 0.125 * pregion%grid%vol(i)
       JFYCell(i) = JFYCell(i) * 0.125 * pregion%grid%vol(i)
       JFZCell(i) = JFZCell(i) * 0.125 * pregion%grid%vol(i)
       !JE correction

       JFECell(i) = JFECell(i) * 0.125 * pregion%grid%vol(i)

       !energydotg = JFXCell(i) * ug(1) + JFYCell(i) * ug(2) + JFECell(i)

       energydotg = JFECell(i) ! includs KE feedback already

IF (IsNan(JFXCell(i)) .EQV. .TRUE.) THEN
        write(*,*) "BROKEN-PX",i,JFXCell(i),ug(1),ug(2),ug(3)
        write(*,*) "JFY",i,JFYCell(i)
        write(*,*) "JFZ",i,JFZCell(i)
        write(*,*) "pregionvol", pregion%grid%vol(i)
        CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PX')
endif
IF (IsNan(JFYCell(i)) .EQV. .TRUE.) THEN
        write(*,*) "BROKEN-PY",i,JFYCell(i),ug(1),ug(2),ug(3)
        write(*,*) "pregionvol", pregion%grid%vol(i)
        CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PY')
endif
IF (IsNan(JFZCell(i)) .EQV. .TRUE.) THEN
        write(*,*) "BROKEN-PZ",i,JFZCell(i),ug(1),ug(2),ug(3)
        write(*,*) "pregionvol", pregion%grid%vol(i)
        CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PY')
endif
IF (IsNan(energydotg) .EQV. .TRUE.) THEN
        write(*,*) "BROKEN-PE",energydotg,i,JFXCell(i),ug(1),JFYCell(i),ug(2),pregion%grid%vol(i),pRegion%mixt%piclGeom
        write(*,*) "pregionvol", pregion%grid%vol(i)
        CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PE')
endif

        pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                         + JFXCell(i)
        
        pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                         + JFYCell(i)

        pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                         + JFZCell(i)

        pRegion%mixt%rhs(CV_MIXT_ENER,i) &
                         = pRegion%mixt%rhs(CV_MIXT_ENER,i) &
                         + energydotg
    END DO

END IF ! global%piclFeedbackFlag

!SOLVE
     CALL ppiclf_solve_IntegrateParticle(1,piclIO,piclDtMin,piclCurrentTime)

!

!Due to moving particle integration stuff stoping this for now
DO i = 1,pRegion%grid%nCells
!zero out PhiP
       PhiP(i) = 0 
       do lz=1,2
       do ly=1,2
       do lx=1,2 
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JPHIP,vfP(lx,ly,lz,i))
       PhiP(i) = PhiP(i) +  (0.125*vfP(lx,ly,lz,i))*(pRegion%grid%vol(i))
       end do
       end do
       end do 
       !Particles have moved  
       !Dump back VolFrac
       PhiP(i) = PhiP(i) / (pRegion%grid%vol(i))
!VOL Frac Cap
       if (PhiP(i) .gt. 0.62) PhiP(i) = 0.62
       pRegion%mixt%piclVF(i) = PhiP(i) 
end DO


!Deallocate arrays

    IF (pRegion%mixtInput%axiFlag) THEN
      DEALLOCATE(YTEMP,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
      END IF ! global%error
    ENDIF

    DEALLOCATE(rhoF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uxF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uyF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uzF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(csF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(tpF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(ppF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(vfP,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpxF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpyF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpzF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
        
    DEALLOCATE(SDRX,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDRY,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDRZ,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(rhsR,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(pGcX,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(pGcY,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF !global%error    

    DEALLOCATE(pGcZ,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF !global%error    

    DEALLOCATE(JFX,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFXCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFY,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFYCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFZ,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFZCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFE,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFECell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error    

    DEALLOCATE(PhiP,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(domgdx,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(domgdy,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(domgdz,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(drhodx,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(drhody,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(drhodz,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpvxF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpvyF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpvzF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDOX,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDOY,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDOZ,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error



#endif
!PPICLF Integration END

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PICL_TEMP_Runge

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PICL_.F90,v $
!
!
!******************************************************************************

