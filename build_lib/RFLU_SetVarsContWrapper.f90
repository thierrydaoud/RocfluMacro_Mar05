










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
! ******************************************************************************
!
! Purpose: Wrapper for setting continuous-phase variables.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  icgBeg       Beginning cell index
!  icgEnd       Ending cell index
!
! Output: None.
!
! Notes: 
!   1. Split from RFLU_SetVarsWrapper because in RungeKuttaMP need to update
!      continuous-phase variables separately and discrete-phase variables can 
!      only be set once their communication is completed.
!
! ******************************************************************************
!
! $Id: RFLU_SetVarsContWrapper.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetVarsContWrapper(pRegion,icgBeg,icgEnd)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid

  USE RFLU_ModConvertCv, ONLY: RFLU_ScalarConvertCvCons2Prim, & 
                               RFLU_ScalarConvertCvPrim2Cons

  USE ModInterfaces, ONLY: RFLU_SetVars

  USE ModInterfacesSpecies, ONLY: SPEC_UpdateDependentVars
                                 
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: icgBeg,icgEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetVarsContWrapper.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetVarsContWrapper',"../libflu/RFLU_SetVarsContWrapper.F90")

! ******************************************************************************
! Set variables
! ******************************************************************************

! ==============================================================================
! Mixture. NOTE need species conserved state vector in primitive form
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)
  END IF ! global%specUsed

  CALL RFLU_SetVars(pRegion,icgBeg,icgEnd)

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)
  END IF ! global%specUsed

! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_UpdateDependentVars(pRegion,icgBeg,icgEnd)
  END IF ! global%specUsed

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetVarsContWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetVarsContWrapper.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:36  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:50  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.1  2005/11/11 13:26:43  fnajjar
! Initial revision
!
! ******************************************************************************

