










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
! Purpose: compute interaction source for drag forces on Lagrangian particles.
!
! Description: none.
!
! Input: region  = current region.
!
! Output: region%levels(iLev)%plag%inrtSources
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_CalcDragUnsteady_AMImplicit.F90,v 1.3 2016/02/08 22:26:47 rahul Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_CalcDragUnsteady_AMImplicit( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag
  USE ModError
  USE ModParameters
  USE INRT_ModParameters


  USE RFLU_ModDifferentiationCells, ONLY: RFLU_ComputeGradCellsGGScalar, &
                                          RFLU_ComputeGradCellsGGVector

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: dragUnsteady, iCont, nCont, nPcls
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass

  REAL(RFREAL) :: CdTotal,diamL,factor,gamma,machL,massL,mixtVolR,pi,psiL, &
                  relVelMagL,reyL,tauLR
  REAL(RFREAL),          DIMENSION(3)   :: relVel, accelL
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pDv, pTv
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: dudtMixt,dudtPlag

  TYPE(t_plag)  , POINTER :: pPlag
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! New variables for unsteady force ----------
  INTEGER :: icg,iT,nUnsteadyData
  REAL(RFREAL), DIMENSION(3) :: forceIU,forcePG,forceTotal,forceVU
  REAL(RFREAL) :: A,B,CamEff,dt,fH,kernelVU,mf,mu,nu,refArea,rhoMixt, &
                  speedSound,time,vFrac,vFracCorr,volL,volMixt
  ! Subbu - cyldet case
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  REAL(RFREAL) :: ConvFlux,rho,tBeg,tEnd
  INTEGER :: Coord
  ! Subbu - End cyldet case 

! New variables for augmenting rhs ----------
  REAL(RFREAL) :: coeffIU,contFac,energydotg,energydotp
  REAL(RFREAL) :: drudtMixt,drvdtMixt,drwdtMixt,drudtPlag,drvdtPlag,drwdtPlag 
  REAL(RFREAL), DIMENSION(3) :: ug,up
! -------------------------------------------
!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_CalcDragUnsteady_AMImplicit.F90,v $'

  global => region%global
  pRegion => region

  CALL RegisterFunction( global,'INRT_CalcDragUnsteady_AMImplicit',"../rocinteract/INRT_CalcDragUnsteady_AMImplicit.F90" )

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_CalcDragUnsteady_AMImplicit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_CalcDragUnsteady_AMImplicit.F90,v $
! Revision 1.3  2016/02/08 22:26:47  rahul
! Subtracted work done contribution from F_pg to fluid phase energy equation.
!
! Revision 1.2  2015/12/19 00:26:46  rahul
! Suppressed the computation of pressure gradient force. This is a
! consequence of governing equations' formulation of multiphase AUSM+up
! scheme.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************

