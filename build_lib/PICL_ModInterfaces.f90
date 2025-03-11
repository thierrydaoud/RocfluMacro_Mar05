










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
! Purpose: set explicit interfaces to subroutines and functions.
!
! Description: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_ModInterfaces.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001-2004 by the University of Illinois
!
!******************************************************************************
  
MODULE PICL_ModInterfaces

  IMPLICIT NONE
  
  INTERFACE

! =============================================================================
! 1 specific library
! =============================================================================

  SUBROUTINE PICL_TEMP_Runge( pRegion)
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE PICL_TEMP_Runge

  SUBROUTINE PICL_TEMP_WriteVTU( pRegion)
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE PICL_TEMP_WriteVTU
       
  SUBROUTINE PICL_TEMP_InitSolver( pRegion)!,global )
    USE ModDataStruct, ONLY : t_region!,t_global
     !TYPE(t_region), INTENT(INOUT) :: pRegion
     TYPE(t_region), POINTER :: pRegion   
     !TYPE(t_global), INTENT(INOUT) :: global
     !TYPE(t_region), DIMENSION(:), POINTER :: regions
  END SUBROUTINE PICL_TEMP_InitSolver

  SUBROUTINE PICL_ReadPiclSection(global)
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE PICL_ReadPiclSection
  

! -----------------------------------------------------------------------------
! 1-specific routines
! -----------------------------------------------------------------------------
  
  
  END INTERFACE

END MODULE PICL_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModInterfaces.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/04/16 23:21:41  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.2  2007/04/15 02:38:04  haselbac
! Added interface for PLAG_RFLU_InitSolSerial_2D.F90
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.46  2007/03/15 22:00:18  haselbac
! Added IFs for PLAG_RFLU_InitSolSerial_{1,3}D.F90
!
! Revision 1.45  2005/05/19 16:01:40  fnajjar
! Removed Interfaces for obsolete routine calls
!
! Revision 1.44  2005/04/27 15:00:08  fnajjar
! Removed interface calls to individual FindCells routines
!
! Revision 1.43  2005/04/25 18:39:10  luca1
! Imposed PDF from file option for random particle ejection
!
! Revision 1.42  2005/03/11 02:23:06  haselbac
! Adapted interfaces for new PLAG_RFLU_FindCellsTrajXYZ routines
!
! Revision 1.41  2005/01/01 21:33:12  haselbac
! Added interface for Apte method, modif interface for PLAG_ReflectParticleData
!
! Revision 1.40  2004/11/05 21:51:01  fnajjar
! Added interfaces to particle-cell search routines
!
! Revision 1.39  2004/11/05 21:09:14  haselbac
! Fixed bad check-in
!
! Revision 1.38  2004/11/05 20:33:57  haselbac
! Adapted interface for PLAG_ReflectParticleData
!
! Revision 1.37  2004/10/11 22:09:22  haselbac
! Renamed procedures
!
! Revision 1.36  2004/10/08 22:09:26  haselbac
! Added entry for PLAG_RFLU_FindParticleCellsBrut
!
! Revision 1.35  2004/08/20 23:27:13  fnajjar
! Added Infrastructure for Plag prep tool
!
! Revision 1.34  2004/06/16 22:58:48  fnajjar
! Renamed injcModel to injcDiamDist for CRE kernel
!
! Revision 1.33  2004/04/09 22:57:45  fnajjar
! Added Interfaces for RFLO-specific routines
!
! Revision 1.32  2004/04/08 01:33:00  haselbac
! Added interface for PLAG_ReflectParticleData
!
! Revision 1.31  2004/03/26 21:26:26  fnajjar
! Added new routines for 1-specific routines
!
! Revision 1.30  2004/03/15 21:05:54  haselbac
! Deleted/added interface
!
! Revision 1.29  2004/03/08 22:17:22  fnajjar
! Included interface calls for 1-specific injection routines
!
! Revision 1.28  2004/03/05 23:17:05  haselbac
! Added interface for PLAG_RFLU_FindParticleCellsTraj
!
! Revision 1.27  2004/02/26 21:02:12  haselbac
! Added common and 1-specific interfaces
!
! Revision 1.26  2004/02/25 21:56:45  fnajjar
! Included generic RKUpdate for PLAG
!
! Revision 1.25  2004/02/06 21:19:36  fnajjar
! Included proper INTENT to Interfaces
!
! Revision 1.24  2004/01/15 21:10:52  fnajjar
! Separated Interfaces for corner-edge cell routines
!
! Revision 1.23  2003/11/21 22:44:35  fnajjar
! Changed global to region for updated Random Number Generator
!
! Revision 1.22  2003/11/12 21:34:36  fnajjar
! Added Corner-Edge cells Interface calls
!
! Revision 1.21  2003/11/03 21:22:20  fnajjar
! Added PLAG_copyFaceVectors
!
! Revision 1.20  2003/09/13 20:14:21  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.19  2003/05/28 15:16:11  fnajjar
! Removed obsolete PLAG_mixt calls as embedded in Rocinteract
!
! Revision 1.18  2003/04/17 00:12:08  fnajjar
! Added INTENT(IN) for ijkNR, ijkNRI, ijkNRJ, ijkNRK
!
! Revision 1.17  2003/04/14 19:02:57  fnajjar
! Bug fix to include POINTER attribute for PLAG_PatchUpdate
!
! Revision 1.16  2003/04/14 14:32:20  fnajjar
! Added PLAG_initInputValues for proper initialization
!
! Revision 1.15  2003/03/28 19:50:39  fnajjar
! Moved interfaces for wrapper routines
!
! Revision 1.14  2003/02/21 17:06:55  fnajjar
! Added Interfaces to Data Send and Recv
!
! Revision 1.13  2003/01/24 22:07:57  f-najjar
! Include Interface call to PLAG_ClearRequests
!
! Revision 1.12  2003/01/23 17:12:14  f-najjar
! Included Interface for PLAG_BufferSizeRecv and PLAG_BufferSizeSend
!
! Revision 1.11  2003/01/23 17:05:06  f-najjar
! Moved Interface call to PLAG_patchBufferSendRecv into ModInterfacesLagrangian
!
! Revision 1.10  2003/01/17 19:32:20  f-najjar
! Included correct INTERFACE for PLAG_PatchRemoveDataOutflow and PLAG_PatchLoadDataBuffers
!
! Revision 1.9  2003/01/16 22:46:18  f-najjar
! Include iReg to calling sequence for PLAG_GetCellIndices
!
! Revision 1.8  2003/01/16 22:30:10  f-najjar
! Included INTERFACE call for PLAG_PatchBufferSendRecv
!
! Revision 1.7  2003/01/16 22:28:40  f-najjar
! Included INTERFACE call for PLAG_AppendDataFromBuffers
!
! Revision 1.6  2003/01/13 19:00:33  f-najjar
! Removed PLAG_allocateDataBuffers
!
! Revision 1.5  2003/01/10 19:37:40  f-najjar
! Included Interface call to PLAG_BoundaryConditionsSet
!
! Revision 1.4  2003/01/10 19:28:04  f-najjar
! Included Interface call for PLAG_PatchExchangeConf
!
! Revision 1.3  2003/01/10 19:09:10  f-najjar
! Included iReg in calling sequence
!
! Revision 1.2  2002/10/25 14:11:12  f-najjar
! Included interface calls to PLAG subroutines
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************

