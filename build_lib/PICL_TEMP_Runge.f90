










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
USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                             RFLU_ConvertCvPrim2Cons

 USE ModInterfaces, ONLY: RFLU_DecideWrite !BRAD added for picl
 



!DEC$ NOFREEFORM

! y, y1, ydot, ydotc: 12

! rprop: 36

! map: 10





















!DEC$ FREEFORM


  IMPLICIT NONE


! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString


TYPE(t_global), POINTER :: global
TYPE(t_level), POINTER :: levels(:)
TYPE(t_region), POINTER :: pRegion
TYPE(t_grid), POINTER :: pGrid
!INTEGER :: errorFlag

  LOGICAL :: doWrite      
  INTEGER(KIND=4) :: i,piclIO,nCells,lx,ly,lz
  INTEGER :: errorFlag,icg      
  REAL(KIND=8) :: piclDtMin,piclCurrentTime, &
          temp_drudtMixt,temp_drvdtMixt,temp_drwdtMixt,energydotg
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

  REAL(KIND=8) :: ppiclf



   
!******************************************************************************

  RCSIdentString = '$RCSfile: PICL_TEMP_Runge.F90,v $ $Revision: 1.0 $'
 
  global => pRegion%global
  
  CALL RegisterFunction( global, 'PICL_TEMP_Runge',"../rocpicl/PICL_TEMP_Runge.F90" )



! Set pointers ----------------------------------------------------------------

    !pRegion => regions!pLevel%regions(iReg)
    pGrid   => pRegion%grid

!PPICLF Integration

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
      CALL ErrorStop(global,ERR_ALLOCATE,204,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uxF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,210,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uyF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,216,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uzF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,222,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(csF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,228,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(tpF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,234,'PPICLF:xGrid')
    END IF ! global%error    

    ALLOCATE(ppF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,240,'PPICLF:xGrid')
    END IF ! global%error    

    ALLOCATE(vfP(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,246,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpxF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,252,'PPICLF:xGrid')
    END IF ! global%error
    
    ALLOCATE(dpyF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,258,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpzF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,264,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRX(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,270,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRY(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,276,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRZ(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,282,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(rhsR(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,288,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcX(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,294,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcY(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,300,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcZ(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,306,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFX(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,312,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFXCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,318,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFY(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,324,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFYCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,330,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFZ(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,336,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFZCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,342,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFECell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,348,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFE(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,354,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(PhiP(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,360,'PPICLF:xGrid')
    END IF ! global%error

    IF (pRegion%mixtInput%axiFlag) THEN
      ALLOCATE(YTEMP(2,2,2,nCells),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,367,'PPICLF:xGrid')
      END IF ! global%error
    ENDIF

    ALLOCATE(domgdx(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,374,'PPICLF:xGrid')
    END IF ! global%error
    
    ALLOCATE(domgdy(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,380,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(domgdz(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,386,'PPICLF:xGrid')
    END IF ! global%error


!Might need to update prim like plag does
pGc => pRegion%mixt%gradCell
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

        temp_drudtMixt = -pRegion%mixt%rhs(CV_MIXT_XMOM,i)&
                  +pRegion%mixt%cv(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,2,i)) &
                  +ug(XCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))

        temp_drvdtMixt = -pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                  +pRegion%mixt%cv(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,3,i))&
                  +ug(YCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))

        temp_drwdtMixt = -pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                +pRegion%mixt%cv(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,4,i))&
                +ug(ZCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))



       do lz=1,2
       do ly=1,2
       do lx=1,2 
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,1,vfP(lx,ly,lz,i))
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

       dpxF(lx,ly,lz,i) = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_PRES,i)
       dpyF(lx,ly,lz,i) = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_PRES,i) 
       dpzF(lx,ly,lz,i) = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_PRES,i) 

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

 
       SDRX(lx,ly,lz,i) = temp_drudtMixt 
       SDRY(lx,ly,lz,i) = temp_drvdtMixt 
       SDRZ(lx,ly,lz,i) = temp_drwdtMixt 

       rhsR(lx,ly,lz,i) = -pRegion%mixt%rhs(CV_MIXT_DENS,i)

       pGcX(lx,ly,lz,i) = pGc(XCOORD,1,i) ! dr/dx
       pGcY(lx,ly,lz,i) = pGc(YCOORD,1,i) ! dr/dy
       pGcz(lx,ly,lz,i) = pGc(ZCOORD,1,i) ! dr/dz

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
! TLJ 21 in PPICLF_USER.h must match the number
!     of calls to ppiclf_solve_InterpFieldUser
! Davin - added pressure 02/22/2025
      IF (21 .NE. 21) THEN
         write(*,*) "Error: PPICLF_LRP_INT must be set to 21"
         CALL ErrorStop(global,ERR_INVALID_VALUE ,497,'PPICLF:LRP_INT')
      endif

      CALL ppiclf_solve_InterpFieldUser(2,rhoF)
      CALL ppiclf_solve_InterpFieldUser(6,uxF)
      CALL ppiclf_solve_InterpFieldUser(7,uyF)
      CALL ppiclf_solve_InterpFieldUser(8,uzF)
      CALL ppiclf_solve_InterpFieldUser(10,dpxF)
      CALL ppiclf_solve_InterpFieldUser(11,dpyF)  
      CALL ppiclf_solve_InterpFieldUser(12,dpzF)  
      CALL ppiclf_solve_InterpFieldUser(9,csF)
      CALL ppiclf_solve_InterpFieldUser(24,tpF)
      CALL ppiclf_solve_InterpFieldUser(5,vfP)  
      CALL ppiclf_solve_InterpFieldUser(13,SDRX)
      CALL ppiclf_solve_InterpFieldUser(14,SDRY)  
      CALL ppiclf_solve_InterpFieldUser(15,SDRZ)  
      CALL ppiclf_solve_InterpFieldUser(16,rhsR)  
      CALL ppiclf_solve_InterpFieldUser(17,pGcX) 
      CALL ppiclf_solve_InterpFieldUser(18,pGcY) 
      CALL ppiclf_solve_InterpFieldUser(19,pGcZ) 
      CALL ppiclf_solve_InterpFieldUser(31,domgdx)
      CALL ppiclf_solve_InterpFieldUser(32,domgdy)  
      CALL ppiclf_solve_InterpFieldUser(33,domgdz)  
      CALL ppiclf_solve_InterpFieldUser(36,ppF)  


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
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,2,JFX(lx,ly,lz,i))  
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,3,JFY(lx,ly,lz,i))
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,4,JFZ(lx,ly,lz,i))
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,5,JFE(lx,ly,lz,i))      
       !call get energy  
       JFXCell(i) = JFXCell(i) + JFX(lx,ly,lz,i) ! / pRegion%grid%vol(i)    
       JFYCell(i) = JFYCell(i) + JFY(lx,ly,lz,i) 
       JFZCell(i) = JFZCell(i) + JFZ(lx,ly,lz,i) 
       JFECell(i) = JFECell(i) + JFE(lx,ly,lz,i)  
        !Jenergy = +...
       end do
       end do
       end do 
        JFXCell(i) = JFXCell(i) *0.125 * pregion%grid%vol(i)
        JFYCell(i) = JFYCell(i) *0.125 * pregion%grid%vol(i)
        JFZCell(i) = JFZCell(i) *0.125 * pregion%grid%vol(i)
        !JE correction

        JFECell(i) = JFECell(i) * 0.125 * pregion%grid%vol(i)

        !energydotg = JFXCell(i) * ug(1) + JFYCell(i) * ug(2) + JFECell(i)

        energydotg = JFECell(i) ! includs KE feedback already

IF (IsNan(JFXCell(i)) .EQV. .TRUE.) THEN
        write(*,*) "BROKEN-PX",i,JFXCell(i),ug(1),ug(2),ug(3)
        write(*,*) "JFY",i,JFYCell(i)
        write(*,*) "JFZ",i,JFZCell(i)
        write(*,*) "pregionvol", pregion%grid%vol(i)
        CALL ErrorStop(global,ERR_INVALID_VALUE ,571,'PPICLF:Broken PX')
endif
IF (IsNan(JFYCell(i)) .EQV. .TRUE.) THEN
        write(*,*) "BROKEN-PY",i,JFYCell(i),ug(1),ug(2),ug(3)
        write(*,*) "pregionvol", pregion%grid%vol(i)
        CALL ErrorStop(global,ERR_INVALID_VALUE ,576,'PPICLF:Broken PY')
endif
IF (IsNan(JFZCell(i)) .EQV. .TRUE.) THEN
        write(*,*) "BROKEN-PZ",i,JFZCell(i),ug(1),ug(2),ug(3)
        write(*,*) "pregionvol", pregion%grid%vol(i)
        CALL ErrorStop(global,ERR_INVALID_VALUE ,581,'PPICLF:Broken PY')
endif
IF (IsNan(energydotg) .EQV. .TRUE.) THEN
        write(*,*) "BROKEN-PE",energydotg,i,JFXCell(i),ug(1),JFYCell(i),ug(2),pregion%grid%vol(i),pRegion%mixt%piclGeom
        write(*,*) "pregionvol", pregion%grid%vol(i)
        CALL ErrorStop(global,ERR_INVALID_VALUE ,586,'PPICLF:Broken PE')
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
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,1,vfP(lx,ly,lz,i))
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
        CALL ErrorStop(global,ERR_DEALLOCATE,640,'PPICLF:xGrid')
      END IF ! global%error
    ENDIF

    DEALLOCATE(rhoF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,647,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uxF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,653,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uyF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,659,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uzF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,665,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(csF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,671,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(tpF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,677,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(ppF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,683,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(vfP,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,689,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpxF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,695,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpyF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,701,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpzF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,707,'PPICLF:xGrid')
    END IF ! global%error
        
    DEALLOCATE(SDRX,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,713,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDRY,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,719,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDRZ,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,725,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(rhsR,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,731,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(pGcX,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,737,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(pGcY,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,743,'PPICLF:xGrid')
    END IF !global%error    

    DEALLOCATE(pGcZ,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,749,'PPICLF:xGrid')
    END IF !global%error    

    DEALLOCATE(JFX,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,755,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFXCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,761,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFY,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,767,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFYCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,773,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFZ,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,779,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFZCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,785,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFE,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,791,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFECell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,797,'PPICLF:xGrid')
    END IF ! global%error    

    DEALLOCATE(PhiP,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,803,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(domgdx,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,809,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(domgdy,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,815,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(domgdz,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,821,'PPICLF:xGrid')
    END IF ! global%error

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

