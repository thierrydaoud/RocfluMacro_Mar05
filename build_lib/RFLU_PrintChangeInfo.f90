










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
! Purpose: Display minimum and maximum values of change vector for given 
!   local domain.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: N/A.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_PrintChangeInfo.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2000, 2001, 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_PrintChangeInfo(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModTools, ONLY: MakeNonZero

  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region 
  
  USE ModInterfaces, ONLY: RFLU_PrintLocInfo   

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Parameters
! ==============================================================================   

  TYPE(t_region), POINTER :: pRegion 

! ==============================================================================
! Local variables
! ==============================================================================
 
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,ic,iv
  INTEGER :: dummy(1) 
  INTEGER :: loc(CV_MIXT_DENS:CV_MIXT_ENER,MIN_VAL:MAX_VAL)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pRhs
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: rhsn
  TYPE(t_grid), POINTER :: pGrid  
  TYPE(t_global), POINTER :: global
  
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PrintChangeInfo.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PrintChangeInfo',"../rocflu/RFLU_PrintChangeInfo.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing relative change '// & 
                             'information...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal 
  END IF ! global%verbLevel

! ==============================================================================
! Set pointers and allocate memory
! ==============================================================================

  pCv   => pRegion%mixt%cv
  pRhs  => pRegion%mixt%rhs
  pGrid => pRegion%grid

  ALLOCATE(rhsn(CV_MIXT_DENS:CV_MIXT_ENER,1:pGrid%nCells),STAT=errorFlag)
  global%error = errorFlag  
  IF ( global%error /= 0 ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,141,'rhsn')
  END IF ! global

! ==============================================================================
! Compute relative change
! ==============================================================================

  DO ic = 1,pGrid%nCells
    DO iv = CV_MIXT_DENS,CV_MIXT_ENER
      rhsn(iv,ic) = pRhs(iv,ic)/MakeNonZero(pCv(iv,ic))
    END DO ! iv
  END DO ! ic

! ==============================================================================
! Find locations of extrema and print information on extrema: NOTE Asinine 
! coding needed because of poor FORTRAN interface for MINLOC and MAXLOC
! functions... 
! ==============================================================================

  dummy                     = MINLOC(rhsn(CV_MIXT_DENS,1:pGrid%nCells))
  loc(CV_MIXT_DENS,MIN_VAL) = dummy(1)

  dummy                     = MINLOC(rhsn(CV_MIXT_XMOM,1:pGrid%nCells))
  loc(CV_MIXT_XMOM,MIN_VAL) = dummy(1)
  
  dummy                     = MINLOC(rhsn(CV_MIXT_YMOM,1:pGrid%nCells))
  loc(CV_MIXT_YMOM,MIN_VAL) = dummy(1)
  
  dummy                     = MINLOC(rhsn(CV_MIXT_ZMOM,1:pGrid%nCells))
  loc(CV_MIXT_ZMOM,MIN_VAL) = dummy(1)
  
  dummy                     = MINLOC(rhsn(CV_MIXT_ENER,1:pGrid%nCells))
  loc(CV_MIXT_ENER,MIN_VAL) = dummy(1)


  dummy                     = MAXLOC(rhsn(CV_MIXT_DENS,1:pGrid%nCells))
  loc(CV_MIXT_DENS,MAX_VAL) = dummy(1)

  dummy                     = MAXLOC(rhsn(CV_MIXT_XMOM,1:pGrid%nCells))
  loc(CV_MIXT_XMOM,MAX_VAL) = dummy(1)
  
  dummy                     = MAXLOC(rhsn(CV_MIXT_YMOM,1:pGrid%nCells))
  loc(CV_MIXT_YMOM,MAX_VAL) = dummy(1)
  
  dummy                     = MAXLOC(rhsn(CV_MIXT_ZMOM,1:pGrid%nCells))
  loc(CV_MIXT_ZMOM,MAX_VAL) = dummy(1)
  
  dummy                     = MAXLOC(rhsn(CV_MIXT_ENER,1:pGrid%nCells))
  loc(CV_MIXT_ENER,MAX_VAL) = dummy(1)
  

  WRITE(STDOUT,'(A,3X,A,2(1X,E15.8),2(1X,I9))') SOLVER_NAME,'Mass:      ', & 
                MINVAL(rhsn(CV_MIXT_DENS,1:pGrid%nCells)), & 
                MAXVAL(rhsn(CV_MIXT_DENS,1:pGrid%nCells)), & 
                loc(CV_MIXT_DENS,MIN_VAL),loc(CV_MIXT_DENS,MAX_VAL)
  WRITE(STDOUT,'(A,3X,A,2(1X,E15.8),2(1X,I9))') SOLVER_NAME,'X-momentum:', & 
                MINVAL(rhsn(CV_MIXT_XVEL,1:pGrid%nCells)), & 
                MAXVAL(rhsn(CV_MIXT_XVEL,1:pGrid%nCells)), & 
                loc(CV_MIXT_XVEL,MIN_VAL),loc(CV_MIXT_XVEL,MAX_VAL)
  WRITE(STDOUT,'(A,3X,A,2(1X,E15.8),2(1X,I9))') SOLVER_NAME,'Y-momentum:', & 
                MINVAL(rhsn(CV_MIXT_YVEL,1:pGrid%nCells)), & 
                MAXVAL(rhsn(CV_MIXT_YVEL,1:pGrid%nCells)), & 
                loc(CV_MIXT_YVEL,MIN_VAL),loc(CV_MIXT_YVEL,MAX_VAL)
  WRITE(STDOUT,'(A,3X,A,2(1X,E15.8),2(1X,I9))') SOLVER_NAME,'Z-momentum:', & 
                MINVAL(rhsn(CV_MIXT_ZVEL,1:pGrid%nCells)), & 
                MAXVAL(rhsn(CV_MIXT_ZVEL,1:pGrid%nCells)), & 
                loc(CV_MIXT_ZVEL,MIN_VAL),loc(CV_MIXT_ZVEL,MAX_VAL)
  WRITE(STDOUT,'(A,3X,A,2(1X,E15.8),2(1X,I9))') SOLVER_NAME,'Energy:    ', & 
                MINVAL(rhsn(CV_MIXT_ENER,1:pGrid%nCells)), & 
                MAXVAL(rhsn(CV_MIXT_ENER,1:pGrid%nCells)), & 
                loc(CV_MIXT_ENER,MIN_VAL),loc(CV_MIXT_ENER,MAX_VAL)
    
! ==============================================================================
! Print out locations of cells at which extrema occur
! ==============================================================================  
  
  IF ( global%verbLevel /= VERBOSE_NONE ) THEN 
    CALL RFLU_PrintLocInfo(pRegion,loc,CV_MIXT_ENER-CV_MIXT_DENS+1, & 
                           LOCINFO_MODE_SILENT,OUTPUT_MODE_MASTER_ONLY)
  END IF ! global%verbLevel
  
! ==============================================================================
! Deallocate memory
! ==============================================================================  
  
  DEALLOCATE(rhsn,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,229,'rhsn')
  END IF ! global  
  
  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing relative change '// & 
                             'information done.'
  END IF ! global%verbLevel    
  
  CALL DeregisterFunction(global)  
  
! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE RFLU_PrintChangeInfo


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_PrintChangeInfo.F90,v $
!   Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
!   merged rocflu micro and macro
!
!   Revision 1.1.1.1  2014/07/15 14:31:38  brollin
!   New Stable version
!
!   Revision 1.3  2008/12/06 08:43:48  mtcampbe
!   Updated license.
!
!   Revision 1.2  2008/11/19 22:17:00  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.1  2007/04/09 18:49:57  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.1  2007/04/09 18:01:01  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.10  2004/01/22 16:04:33  haselbac
!   Changed declaration to eliminate warning on ALC
!
!   Revision 1.9  2003/06/04 22:43:01  haselbac
!   Adapted call to RFLU_PrintLocInfo
!
!   Revision 1.8  2003/03/15 18:51:54  haselbac
!   Adapted call to function
!
!   Revision 1.7  2003/01/28 14:46:30  haselbac
!   Cosmetics only
!
!   Revision 1.6  2002/10/08 15:49:29  haselbac
!   {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
!   Revision 1.5  2002/09/09 15:51:56  haselbac
!   global now under region
!
!   Revision 1.4  2002/07/25 14:29:34  haselbac
!   Cosmetic changes to output
!
!   Revision 1.3  2002/06/17 13:34:12  haselbac
!   Prefixed SOLVER_NAME to all screen output
!
!   Revision 1.2  2002/06/10 21:31:59  haselbac
!   Now printing relative changes
!
!   Revision 1.1  2002/06/05 18:56:48  haselbac
!   Initial revision
!
! ******************************************************************************

