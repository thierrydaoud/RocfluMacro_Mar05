










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
! $Id: INRT_CalcDrag.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_CalcDrag( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: icg, dragModel, iCont, nCont, nPcls
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass

  REAL(RFREAL) :: C1,C2,C3,C4,CdFinal,CdSTD,CdM1,CdM2,CdMcr,DeltaCd,diamL,f1M, &
                  f1M1,f1M2,f2M,f2M1,f2M2,f3M,f3M1,f3M2,f4M,factor,gamma,lre, &
                  Mach,Mach1,Mach2,massL,mixtVolR,pi,psiL,relVelMagL,reyL, &
                  tauLR,vFrac,vFracCorr
  REAL(RFREAL),          DIMENSION(3)   :: relVel, accelL
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pDv, pTv

  TYPE(t_plag)  , POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_CalcDrag.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_CalcDrag',"../rocinteract/INRT_CalcDrag.F90" )

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_CalcDrag

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_CalcDrag.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2009/07/09 20:44:19  mparmar
! Added Parmar drag law
!
! Revision 1.3  2008/12/06 08:43:49  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:01  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:11  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:14  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2007/03/08 14:59:36  fnajjar
! Fixed bug in Sommerfeld drag law for psiL
!
! Revision 1.2  2007/03/07 22:16:47  fnajjar
! Added Sommerfeld drag law
!
! Revision 1.1  2004/12/01 21:56:14  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2004/01/31 03:59:22  haselbac
! Initial integration for Rocflu and Rocpart
!
! Revision 1.4  2003/04/03 21:10:17  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.3  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************

