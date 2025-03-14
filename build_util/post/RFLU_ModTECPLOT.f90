










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
! Purpose: Collect routines to write TECPLOT file.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes:
!   1. Note the various array format required by TECPLOT - e.g., coordinate-
!      array must first be converted.
!   2. Routines were separated so that can open TECPLOT file and write several
!      zones to it. This is useful to visualize parallel results.
!
! ******************************************************************************
!
! $Id: RFLU_ModTECPLOT.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModTECPLOT

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region

!  USE SPEC_ModParameters

  USE RFLU_ModTECPLOTUtils

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PARAMETER, PRIVATE :: &
    RCSIdentString = '$RCSfile: RFLU_ModTECPLOT.F90,v $ $Revision: 1.1.1.1 $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_TEC_BuildDataFieldSurf, &
            RFLU_TEC_BuildDataFieldVol, &
            RFLU_TEC_BuildDataPatch, &
            RFLU_TEC_BuildDataPatchStats, &
            RFLU_TEC_CloseFileField, &
            RFLU_TEC_CloseFilePatch, &
            RFLU_TEC_CloseFilePatchStats, &
            RFLU_TEC_CloseFilePnt, &
            RFLU_TEC_DestroyDataFieldSurf, &
            RFLU_TEC_DestroyDataFieldVol, &
            RFLU_TEC_DestroyDataPatch, &
            RFLU_TEC_Init, &
            RFLU_TEC_OpenFileField, &
            RFLU_TEC_OpenFilePatch, &
            RFLU_TEC_OpenFilePatchStats, &
            RFLU_TEC_OpenFilePnt, & 
            RFLU_TEC_WriteFileFieldSurf, &
            RFLU_TEC_WriteFileFieldVol, &
            RFLU_TEC_WriteFilePatch, &
            RFLU_TEC_WriteFilePnt

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_TEC_BuildHeaderField, &
             RFLU_TEC_BuildHeaderPatch, & 
             RFLU_TEC_BuildHeaderPatchStats 

! ******************************************************************************
! Subroutines and functions
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Collect surface data for writing to TECPLOT field file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   pPatch              Pointer to patch data
!
! Output: None.
!
! Notes:
!   1. Isolated this code so that can easily add data for writing to file.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_BuildDataFieldSurf(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch  
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,icg,ifl,iPatch,iVar,iVarFaceTot,iVarTot,iVarVertTot, & 
             ivg,ivl,nVarsFaceTot,nVarsFaceTotSave,nVarsVertTot, &
             nVarsVertTotSave
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pVar
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECFIL100

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_BuildDataFieldSurf', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Building surface data for TECPLOT field file...'
  END IF ! global%verbLevel


! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Set TECPLOT file context to field file
! ******************************************************************************

  errorFlag = TECFIL100(fileType2CntrTEC(FILE_TYPE_FIELD))
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,226)
  END IF ! global%error

! ******************************************************************************
! Count number of variables
! ******************************************************************************

! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  nVarsTEC = 3

! ==============================================================================
! Solution
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN

! ------------------------------------------------------------------------------
!   Mixture
! ------------------------------------------------------------------------------

    nVarsTEC = nVarsTEC + pRegion%plot%nCvMixt & 
                        + pRegion%plot%nDvMixt &
                        + pRegion%plot%nGvMixt &
                        + pRegion%plot%nPv                   

! ------------------------------------------------------------------------------
!   Physical modules
! ------------------------------------------------------------------------------

    IF ( global%specUsed .EQV. .TRUE. ) THEN
      nVarsTEC = nVarsTEC + pRegion%specInput%nSpecies & 
!                          + pRegion%specInput%nSpeciesEE*EEV_SPEC_NVAR
                          + pRegion%specInput%nSpeciesEE*4
    END IF ! global%specUsed
  END IF ! global%postPlotType

! ******************************************************************************
! Set position of variables
! ******************************************************************************

  ALLOCATE(posTEC(nVarsTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,274,'posTEC')
  END IF ! global%error

! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  iVarTot = 0

  DO iVar = 1,3
    iVarTot = iVarTot + 1
    posTEC(iVarTot) = VAR_POS_VERT
  END DO ! iVar

! ==============================================================================
! Mixture
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
  
! ------------------------------------------------------------------------------
!   Vertex data
! ------------------------------------------------------------------------------  

    IF ( global%postInterpType /= INTERP_TYPE_NONE ) THEN 
      DO iVar = 1,pRegion%plot%nCvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_VERT
      END DO ! iVar

      DO iVar = 1,pRegion%plot%nDvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_VERT
      END DO ! iVar

      DO iVar = 1,pRegion%plot%nGvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_VERT
      END DO ! iVar

      DO iVar = 1,pRegion%plot%nPv
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_VERT
      END DO ! iVar

! ------------------------------------------------------------------------------
!   Face data
! ------------------------------------------------------------------------------  
  
    ELSE 
      DO iVar = 1,pRegion%plot%nCvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_FACE
      END DO ! iVar

      DO iVar = 1,pRegion%plot%nDvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_FACE
      END DO ! iVar    

      DO iVar = 1,pRegion%plot%nGvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_FACE
      END DO ! iVar
      
      DO iVar = 1,pRegion%plot%nPv
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_FACE
      END DO ! iVar    
    END IF ! global%postInterpType
  END IF ! global%postPlotType

! ==============================================================================
! Physical modules
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN

! ------------------------------------------------------------------------------
!   Vertex data
! ------------------------------------------------------------------------------  

    IF ( global%postInterpType /= INTERP_TYPE_NONE ) THEN   
      IF ( global%specUsed .EQV. .TRUE. ) THEN
        DO iVar = 1,pRegion%specInput%nSpecies
          iVarTot = iVarTot + 1
          posTEC(iVarTot) = VAR_POS_VERT
        END DO ! iVar
        
!        DO iVar = 1,pRegion%specInput%nSpeciesEE*EEV_SPEC_NVAR
        DO iVar = 1,pRegion%specInput%nSpeciesEE*4   
          iVarTot = iVarTot + 1
          posTEC(iVarTot) = VAR_POS_VERT
        END DO ! iVar
      END IF ! global%specUsed

! ------------------------------------------------------------------------------
!   Face data
! ------------------------------------------------------------------------------  
  
    ELSE    
      IF ( global%specUsed .EQV. .TRUE. ) THEN
        DO iVar = 1,pRegion%specInput%nSpecies
          iVarTot = iVarTot + 1
          posTEC(iVarTot) = VAR_POS_FACE
        END DO ! iVar
        
!        DO iVar = 1,pRegion%specInput%nSpeciesEE*EEV_SPEC_NVAR
        DO iVar = 1,pRegion%specInput%nSpeciesEE*4   
          iVarTot = iVarTot + 1
          posTEC(iVarTot) = VAR_POS_FACE
        END DO ! iVar
      END IF ! global%specUsed   
    END IF ! global%postInterpType
  END IF ! global%postPlotType

! ******************************************************************************
! Determine how many variables of each position
! ******************************************************************************

  nVarsFaceTEC = 0
  nVarsVertTEC = 0

  DO iVar = 1,nVarsTEC
    IF ( posTEC(iVar) == VAR_POS_FACE ) THEN
      nVarsFaceTEC = nVarsFaceTEC + 1
    ELSE IF ( posTEC(iVar) == VAR_POS_VERT ) THEN
      nVarsVertTEC = nVarsVertTEC + 1
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,407)
    END IF ! posTEC
  END DO ! iVar

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  ALLOCATE(pPatch%varVertTEC(pPatch%nBVertTot,nVarsVertTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,418,'pPatch%varVertTEC')
  END IF ! global%error

  IF ( nVarsFaceTEC > 0 ) THEN
    ALLOCATE(pPatch%varFaceTEC(pPatch%nBFacesTot,nVarsFaceTEC),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,425,'pPatch%varFaceTEC')
    END IF ! global%error
  ELSE
    NULLIFY(pPatch%varFaceTEC)
  END IF ! nVarsFace

! ******************************************************************************
! Assemble data
! ******************************************************************************

  nVarsFaceTot = 0
  nVarsVertTot = 0
  
! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  DO ivl = 1,pPatch%nBVertTot
    ivg = pPatch%bv(ivl)

    pPatch%varVertTEC(ivl,1) = pGrid%xyz(XCOORD,ivg)
    pPatch%varVertTEC(ivl,2) = pGrid%xyz(YCOORD,ivg)
    pPatch%varVertTEC(ivl,3) = pGrid%xyz(ZCOORD,ivg)
  END DO ! ivl

  nVarsVertTot = 3
    
! ==============================================================================
! Variables
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN

! ------------------------------------------------------------------------------
!   Mixture
! ------------------------------------------------------------------------------

! - Vertex data ----------------------------------------------------------------

    IF ( global%postInterpType /= INTERP_TYPE_NONE ) THEN 
      DO ivl = 1,pPatch%nBVertTot
        ivg = pPatch%bv(ivl)

        iVarVertTot = nVarsVertTot

        DO iVar = 1,pRegion%plot%nCvMixt
          pPatch%varVertTEC(ivl,iVarVertTot+iVar) = & 
            pRegion%mixt%cvVert(iVar,ivg)
        END DO ! iVar
      END DO ! ivl        

      IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
        DO ivl = 1,pPatch%nBVertTot
          ivg = pPatch%bv(ivl)

          iVarVertTot = nVarsVertTot + pRegion%plot%nCvMixt

          DO iVar = 1,pRegion%plot%nDvMixt
            pPatch%varVertTEC(ivl,iVarVertTot+iVar) = & 
              pRegion%mixt%dvVert(DV_MIXT_TEMP,ivg)
          END DO ! iVar
        END DO ! ivl        
      ELSE
        DO ivl = 1,pPatch%nBVertTot
          ivg = pPatch%bv(ivl)

          iVarVertTot = nVarsVertTot + pRegion%plot%nCvMixt

          DO iVar = 1,pRegion%plot%nDvMixt
            pPatch%varVertTEC(ivl,iVarVertTot+iVar) = & 
              pRegion%mixt%dvVert(iVar,ivg)
          END DO ! iVar
        END DO ! ivl        
      END IF ! global%solverType

      DO ivl = 1,pPatch%nBVertTot
        ivg = pPatch%bv(ivl)

        iVarVertTot = nVarsVertTot + pRegion%plot%nCvMixt & 
                                   + pRegion%plot%nDvMixt
       
        DO iVar = 1,pRegion%plot%nGvMixt
          pPatch%varVertTEC(ivl,iVarVertTot+iVar) = & 
            pRegion%mixt%gvVert(iVar,ivg)
        END DO ! iVar        
       
        iVarVertTot = nVarsVertTot + pRegion%plot%nCvMixt &
                                   + pRegion%plot%nDvMixt & 
                                   + pRegion%plot%nGvMixt
             
        DO iVar = 1,pRegion%plot%nPv
          pPatch%varVertTEC(ivl,iVarVertTot+iVar) = & 
            pRegion%plot%pvVert(iVar,ivg)
        END DO ! iVar             
      END DO ! ivl        

! - Face data ------------------------------------------------------------------

    ELSE       
      DO ifl = 1,pPatch%nBFacesTot
        icg = pPatch%bf2c(ifl)

        iVarFaceTot = nVarsFaceTot

        DO iVar = 1,pRegion%plot%nCvMixt
          pPatch%varFaceTEC(ifl,iVarFaceTot+iVar) = pRegion%mixt%cv(iVar,icg)
        END DO ! iVar
      END DO ! ifl                                  

      IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
        DO ifl = 1,pPatch%nBFacesTot
          icg = pPatch%bf2c(ifl)

          iVarFaceTot = nVarsFaceTot + pRegion%plot%nCvMixt

          DO iVar = 1,pRegion%plot%nDvMixt
            pPatch%varFaceTEC(ifl,iVarFaceTot+iVar) = &
                                             pRegion%mixt%dv(DV_MIXT_TEMP,icg)
          END DO ! iVar
        END DO ! ifl                                  
      ELSE
        DO ifl = 1,pPatch%nBFacesTot
          icg = pPatch%bf2c(ifl)

          iVarFaceTot = nVarsFaceTot + pRegion%plot%nCvMixt

          DO iVar = 1,pRegion%plot%nDvMixt
            pPatch%varFaceTEC(ifl,iVarFaceTot+iVar) = pRegion%mixt%dv(iVar,icg)
          END DO ! iVar
        END DO ! ifl                                  
      END IF ! global%solverType

      DO ifl = 1,pPatch%nBFacesTot
        icg = pPatch%bf2c(ifl)

        iVarFaceTot = nVarsFaceTot + pRegion%plot%nCvMixt & 
                                   + pRegion%plot%nDvMixt

        DO iVar = 1,pRegion%plot%nGvMixt
          pPatch%varFaceTEC(ifl,iVarFaceTot+iVar) = pRegion%mixt%gv(iVar,icg)
        END DO ! iVar

        iVarFaceTot = nVarsFaceTot + pRegion%plot%nCvMixt &
                                   + pRegion%plot%nDvMixt & 
                                   + pRegion%plot%nGvMixt

        DO iVar = 1,pRegion%plot%nPv
          pPatch%varFaceTEC(ifl,iVarFaceTot+iVar) = pRegion%plot%pv(iVar,icg)
        END DO ! iVar
      END DO ! ifl                                  
    END IF ! global%postInterpType

    IF ( global%postInterpType /= INTERP_TYPE_NONE ) THEN 
      nVarsVertTot = nVarsVertTot + pRegion%plot%nCvMixt & 
                                  + pRegion%plot%nDvMixt &
                                  + pRegion%plot%nGvMixt &  
                                  + pRegion%plot%nPv                                                                      
    ELSE 
      nVarsFaceTot = nVarsFaceTot + pRegion%plot%nCvMixt & 
                                  + pRegion%plot%nDvMixt &
                                  + pRegion%plot%nGvMixt &    
                                  + pRegion%plot%nPv      
    END IF ! global%postInterpType
    
! ------------------------------------------------------------------------------
!   Physical modules
! ------------------------------------------------------------------------------

    nVarsFaceTotSave = nVarsFaceTot
    nVarsVertTotSave = nVarsVertTot
      
    IF ( global%postInterpType /= INTERP_TYPE_NONE ) THEN 
      
! --- Vertex data --------------------------------------------------------------      

      nVarsVertTot = nVarsVertTotSave

      IF ( global%specUsed .EQV. .TRUE. ) THEN
        DO ivl = 1,pPatch%nBVertTot
          ivg = pPatch%bv(ivl)

          iVarVertTot = nVarsVertTot

          DO iVar = 1,pRegion%specInput%nSpecies
            pPatch%varVertTEC(ivl,iVarVertTot+iVar) & 
              = pRegion%spec%cvVert(iVar,ivg)
          END DO ! iVar
        END DO ! ivl

        nVarsVertTot = nVarsVertTot + pRegion%specInput%nSpecies

        DO ivl = 1,pPatch%nBVertTot
          ivg = pPatch%bv(ivl)
          
          DO iVar = 1,pRegion%specInput%nSpeciesEE
            iVarVertTot = nVarsVertTot

!            pPatch%varVertTEC(ivl,iVarVertTot+EEV_SPEC_NVAR*(iVar-1)+1) & 
!              = pRegion%spec%eevVert(EEV_SPEC_XVEL,iVar,ivg)
!            pPatch%varVertTEC(ivl,iVarVertTot+EEV_SPEC_NVAR*(iVar-1)+2) & 
!              = pRegion%spec%eevVert(EEV_SPEC_YVEL,iVar,ivg)
!            pPatch%varVertTEC(ivl,iVarVertTot+EEV_SPEC_NVAR*(iVar-1)+3) & 
!              = pRegion%spec%eevVert(EEV_SPEC_ZVEL,iVar,ivg)
!            pPatch%varVertTEC(ivl,iVarVertTot+EEV_SPEC_NVAR*(iVar-1)+4) & 
!              = pRegion%spec%eevVert(EEV_SPEC_TEMP,iVar,ivg)
            pPatch%varVertTEC(ivl,iVarVertTot+4*(iVar-1)+1) & 
              = pRegion%spec%eevVert(1,iVar,ivg)
            pPatch%varVertTEC(ivl,iVarVertTot+4*(iVar-1)+2) & 
              = pRegion%spec%eevVert(2,iVar,ivg)
            pPatch%varVertTEC(ivl,iVarVertTot+4*(iVar-1)+3) & 
              = pRegion%spec%eevVert(3,iVar,ivg)
            pPatch%varVertTEC(ivl,iVarVertTot+4*(iVar-1)+4) & 
              = pRegion%spec%eevVert(4,iVar,ivg)                                                        
          END DO ! iVar
        END DO ! ivl

!        nVarsVertTot = nVarsVertTot + pRegion%specInput%nSpecies*EEV_SPEC_NVAR
        nVarsVertTot = nVarsVertTot + pRegion%specInput%nSpecies*4 
      END IF ! global%specUsed

! - Face data ------------------------------------------------------------------

    ELSE 
      nVarsFaceTot = nVarsFaceTotSave

      IF ( global%specUsed .EQV. .TRUE. ) THEN
        DO ifl = 1,pPatch%nBFacesTot
          icg = pPatch%bf2c(ifl)

          iVarFaceTot = nVarsFaceTot

          DO iVar = 1,pRegion%specInput%nSpecies
            pPatch%varFaceTEC(ifl,iVarFaceTot+iVar) & 
              = pRegion%spec%cv(iVar,icg)
          END DO ! iVar
        END DO ! ifl
              
        nVarsFaceTot = nVarsFaceTot + pRegion%specInput%nSpecies 
        
        DO ifl = 1,pPatch%nBFacesTot
          icg = pPatch%bf2c(ifl)
          
          DO iVar = 1,pRegion%specInput%nSpeciesEE
            iVarFaceTot = nVarsFaceTot

!            pPatch%varFaceTEC(ifl,iVarFaceTot+EEV_SPEC_NVAR*(iVar-1)+1) & 
!              = pRegion%spec%eev(EEV_SPEC_XVEL,iVar,icg)
!            pPatch%varFaceTEC(ifl,iVarFaceTot+EEV_SPEC_NVAR*(iVar-1)+2) & 
!              = pRegion%spec%eev(EEV_SPEC_YVEL,iVar,icg)
!            pPatch%varFaceTEC(ifl,iVarFaceTot+EEV_SPEC_NVAR*(iVar-1)+3) & 
!              = pRegion%spec%eev(EEV_SPEC_ZVEL,iVar,icg)
!            pPatch%varFaceTEC(ifl,iVarFaceTot+EEV_SPEC_NVAR*(iVar-1)+4) & 
!              = pRegion%spec%eev(EEV_SPEC_TEMP,iVar,icg) 
            pPatch%varFaceTEC(ifl,iVarFaceTot+4*(iVar-1)+1) & 
              = pRegion%spec%eev(1,iVar,icg)
            pPatch%varFaceTEC(ifl,iVarFaceTot+4*(iVar-1)+2) & 
              = pRegion%spec%eev(2,iVar,icg)
            pPatch%varFaceTEC(ifl,iVarFaceTot+4*(iVar-1)+3) & 
              = pRegion%spec%eev(3,iVar,icg)
            pPatch%varFaceTEC(ifl,iVarFaceTot+4*(iVar-1)+4) & 
              = pRegion%spec%eev(4,iVar,icg)                                                        
          END DO ! iVar
        END DO ! ifl

!        nVarsFaceTot = nVarsFaceTot + pRegion%specInput%nSpecies*EEV_SPEC_NVAR
        nVarsFaceTot = nVarsFaceTot + pRegion%specInput%nSpecies*4                  
      END IF ! global%specUsed
    END IF ! global%postInterpType
    
    IF ( global%postInterpType /= INTERP_TYPE_NONE ) THEN 
      nVarsVertTot = nVarsVertTot + pRegion%specInput%nSpecies
    ELSE 
      nVarsFaceTot = nVarsFaceTot + pRegion%specInput%nSpecies
    END IF ! global%postInterpType    
  END IF ! global%postPlotType

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building surface data for '// &
                             'TECPLOT field file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_BuildDataFieldSurf







! ******************************************************************************
!
! Purpose: Collect volume data for writing to TECPLOT field file.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Isolated this code so that can easily add data for writing to file.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_BuildDataFieldVol(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,icg,iVar,iVarTot,iVarCellTot,iVarVertTot,ivg, & 
             nVarsCellTot,nVarsVertTot
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECFIL100

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_BuildDataFieldVol', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Building volume data for TECPLOT field file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Set TECPLOT file context to field file
! ******************************************************************************

  errorFlag = TECFIL100(fileType2CntrTEC(FILE_TYPE_FIELD))
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,796)
  END IF ! global%error

! ******************************************************************************
! Count number of variables
! ******************************************************************************

! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  nVarsTEC = 3

! ==============================================================================
! Solution
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN

! ------------------------------------------------------------------------------
!   Mixture
! ------------------------------------------------------------------------------

    nVarsTEC = nVarsTEC + pRegion%plot%nCvMixt & 
                        + pRegion%plot%nDvMixt &
                        + pRegion%plot%nGvMixt &                        
                        + pRegion%plot%nPv                          

! ------------------------------------------------------------------------------
!   Physical modules
! ------------------------------------------------------------------------------

    IF ( global%specUsed .EQV. .TRUE. ) THEN
      nVarsTEC = nVarsTEC + pRegion%specInput%nSpecies & 
!                          + pRegion%specInput%nSpeciesEE*EEV_SPEC_NVAR
                          + pRegion%specInput%nSpeciesEE*4
    END IF ! global%specUsed
  END IF ! global%postPlotType

! ******************************************************************************
! Set position of variables
! ******************************************************************************

  ALLOCATE(posTEC(nVarsTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,844,'posTEC')
  END IF ! global%error

! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  iVarTot = 0

  DO iVar = 1,3
    iVarTot = iVarTot + 1
    posTEC(iVarTot) = VAR_POS_VERT
  END DO ! iVar

! ==============================================================================
! Mixture
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
  
! ------------------------------------------------------------------------------  
!   Vertex data
! ------------------------------------------------------------------------------  
  
    IF ( global%postInterpType /= INTERP_TYPE_NONE ) THEN
      DO iVar = 1,pRegion%plot%nCvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_VERT
      END DO ! iVar

      DO iVar = 1,pRegion%plot%nDvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_VERT
      END DO ! iVar
      
      DO iVar = 1,pRegion%plot%nGvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_VERT
      END DO ! iVar         
      
      DO iVar = 1,pRegion%plot%nPv
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_VERT
      END DO ! iVar   
      
! ------------------------------------------------------------------------------  
!   Cell data
! ------------------------------------------------------------------------------  
             
    ELSE 
      DO iVar = 1,pRegion%plot%nCvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_CELL
      END DO ! iVar

      DO iVar = 1,pRegion%plot%nDvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_CELL
      END DO ! iVar    

      DO iVar = 1,pRegion%plot%nGvMixt
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_CELL
      END DO ! iVar 

      DO iVar = 1,pRegion%plot%nPv
        iVarTot = iVarTot + 1
        posTEC(iVarTot) = VAR_POS_CELL
      END DO ! iVar          
    END IF ! global%postInterpType
  END IF ! global%postPlotType

! ==============================================================================
! Physical modules
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN

! ------------------------------------------------------------------------------
!   Vertex data
! ------------------------------------------------------------------------------  
    
    IF ( global%postInterpType /= INTERP_TYPE_NONE ) THEN   
      IF ( global%specUsed .EQV. .TRUE. ) THEN
        DO iVar = 1,pRegion%specInput%nSpecies
          iVarTot = iVarTot + 1
          posTEC(iVarTot) = VAR_POS_VERT
        END DO ! iVar
        
!        DO iVar = 1,pRegion%specInput%nSpeciesEE*EEV_SPEC_NVAR
        DO iVar = 1,pRegion%specInput%nSpeciesEE*4   
          iVarTot = iVarTot + 1
          posTEC(iVarTot) = VAR_POS_VERT
        END DO ! iVar
      END IF ! global%specUsed

! ------------------------------------------------------------------------------
!   Vertex data
! ------------------------------------------------------------------------------  
    
    ELSE 
      IF ( global%specUsed .EQV. .TRUE. ) THEN
        DO iVar = 1,pRegion%specInput%nSpecies
          iVarTot = iVarTot + 1
          posTEC(iVarTot) = VAR_POS_CELL
        END DO ! iVar
        
!        DO iVar = 1,pRegion%specInput%nSpeciesEE*EEV_SPEC_NVAR
        DO iVar = 1,pRegion%specInput%nSpeciesEE*4        
          iVarTot = iVarTot + 1
          posTEC(iVarTot) = VAR_POS_CELL
        END DO ! iVar        
      END IF ! global%specUsed        
    END IF ! global%postInterpType
  END IF ! global%postPlotType

! ******************************************************************************
! Determine how many variables of each position
! ******************************************************************************

  nVarsCellTEC = 0
  nVarsVertTEC = 0

  DO iVar = 1,nVarsTEC
    IF ( posTEC(iVar) == VAR_POS_CELL ) THEN
      nVarsCellTEC = nVarsCellTEC + 1
    ELSE IF ( posTEC(iVar) == VAR_POS_VERT ) THEN
      nVarsVertTEC = nVarsVertTEC + 1
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,977)
    END IF ! posTEC
  END DO ! iVar

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  ALLOCATE(pRegion%varVertTEC(pGrid%nVertTot,nVarsVertTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,988,'pRegion%varVertTEC')
  END IF ! global%error

  IF ( nVarsCellTEC > 0 ) THEN
    ALLOCATE(pRegion%varCellTEC(pGrid%nCellsTot,nVarsCellTEC),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,995,'pRegion%varCellTEC')
    END IF ! global%error
  ELSE
    NULLIFY(pRegion%varCellTEC)
  END IF ! nVarsCellTEC

! ******************************************************************************
! Assemble data
! ******************************************************************************

  nVarsCellTot = 0
  nVarsVertTot = 0
  
! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  DO ivg = 1,pGrid%nVertTot
    pRegion%varVertTEC(ivg,1) = pGrid%xyz(XCOORD,ivg)
    pRegion%varVertTEC(ivg,2) = pGrid%xyz(YCOORD,ivg)
    pRegion%varVertTEC(ivg,3) = pGrid%xyz(ZCOORD,ivg)
  END DO ! ivg

  nVarsVertTot = 3

! ==============================================================================
! Variables
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN

! ------------------------------------------------------------------------------
!   Mixture
! ------------------------------------------------------------------------------

! - Vertex data ----------------------------------------------------------------

    IF ( global%postInterpType /= INTERP_TYPE_NONE ) THEN 
      DO ivg = 1,pGrid%nVertTot
        iVarVertTot = nVarsVertTot

        DO iVar = 1,pRegion%plot%nCvMixt
          pRegion%varVertTEC(ivg,iVarVertTot+iVar) & 
            = pRegion%mixt%cvVert(iVar,ivg)
        END DO ! iVar
      END DO ! ivg

      IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
        DO ivg = 1,pGrid%nVertTot
          iVarVertTot = nVarsVertTot + pRegion%plot%nCvMixt

          DO iVar = 1,pRegion%plot%nDvMixt
            pRegion%varVertTEC(ivg,iVarVertTot+iVar) & 
              = pRegion%mixt%dvVert(DV_MIXT_TEMP,ivg)
          END DO ! iVar
        END DO ! ivg
      ELSE
        DO ivg = 1,pGrid%nVertTot
          iVarVertTot = nVarsVertTot + pRegion%plot%nCvMixt

          DO iVar = 1,pRegion%plot%nDvMixt
            pRegion%varVertTEC(ivg,iVarVertTot+iVar) & 
              = pRegion%mixt%dvVert(iVar,ivg)
          END DO ! iVar
        END DO ! ivg
      END IF ! global%solverType

      DO ivg = 1,pGrid%nVertTot
        iVarVertTot = nVarsVertTot + pRegion%plot%nCvMixt &
                                   + pRegion%plot%nDvMixt

        DO iVar = 1,pRegion%plot%nGvMixt
          pRegion%varVertTEC(ivg,iVarVertTot+iVar) & 
            = pRegion%mixt%gvVert(iVar,ivg)
        END DO ! iVar 
        
        iVarVertTot = nVarsVertTot + pRegion%plot%nCvMixt &
                                   + pRegion%plot%nDvMixt &
                                   + pRegion%plot%nGvMixt       
        
        DO iVar = 1,pRegion%plot%nPv
          pRegion%varVertTEC(ivg,iVarVertTot+iVar) & 
            = pRegion%plot%pvVert(iVar,ivg)
        END DO ! iVar
      END DO ! ivg

      nVarsVertTot = nVarsVertTot + pRegion%plot%nCvMixt & 
                                  + pRegion%plot%nDvMixt &
                                  + pRegion%plot%nGvMixt &
                                  + pRegion%plot%nPv                                 

! - Cell data ------------------------------------------------------------------

    ELSE 
      DO icg = 1,pGrid%nCellsTot
        iVarCellTot = nVarsCellTot

        DO iVar = 1,pRegion%plot%nCvMixt
          pRegion%varCellTEC(icg,iVarCellTot+iVar) = pRegion%mixt%cv(iVar,icg)
        END DO ! iVar
      END DO ! icg

      IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
        DO icg = 1,pGrid%nCellsTot
          iVarCellTot = nVarsCellTot + pRegion%plot%nCvMixt

          DO iVar = 1,pRegion%plot%nDvMixt
            pRegion%varCellTEC(icg,iVarCellTot+iVar) = &
                                         pRegion%mixt%dv(DV_MIXT_TEMP,icg)
          END DO ! iVar
        END DO ! icg
      ELSE
        DO icg = 1,pGrid%nCellsTot
          iVarCellTot = nVarsCellTot + pRegion%plot%nCvMixt

          DO iVar = 1,pRegion%plot%nDvMixt
            pRegion%varCellTEC(icg,iVarCellTot+iVar) = pRegion%mixt%dv(iVar,icg)
          END DO ! iVar
        END DO ! icg
      END IF ! global%solverType

      DO icg = 1,pGrid%nCellsTot
        iVarCellTot = nVarsCellTot + pRegion%plot%nCvMixt & 
                                   + pRegion%plot%nDvMixt

        DO iVar = 1,pRegion%plot%nGvMixt
          pRegion%varCellTEC(icg,iVarCellTot+iVar) = pRegion%mixt%gv(iVar,icg)
        END DO ! iVar  

        iVarCellTot = nVarsCellTot + pRegion%plot%nCvMixt & 
                                   + pRegion%plot%nDvMixt &
                                   + pRegion%plot%nGvMixt  

        DO iVar = 1,pRegion%plot%nPv
          pRegion%varCellTEC(icg,iVarCellTot+iVar) = pRegion%plot%pv(iVar,icg)
        END DO ! iVar            
      END DO ! icg

      nVarsCellTot = nVarsCellTot + pRegion%plot%nCvMixt & 
                                  + pRegion%plot%nDvMixt &
                                  + pRegion%plot%nGvMixt &   
                                  + pRegion%plot%nPv
    END IF ! global%postInterpType 

! ------------------------------------------------------------------------------
!   Physical modules
! ------------------------------------------------------------------------------

! - Vertex data ----------------------------------------------------------------

    IF ( global%postInterpType /= INTERP_TYPE_NONE ) THEN 
      IF ( global%specUsed .EQV. .TRUE. ) THEN
        DO ivg = 1,pGrid%nVertTot
          DO iVar = 1,pRegion%specInput%nSpecies
            iVarVertTot = nVarsVertTot

            pRegion%varVertTEC(ivg,iVarVertTot+iVar) & 
              = pRegion%spec%cvVert(iVar,ivg)
          END DO ! iVar
        END DO ! ivg

        nVarsVertTot = nVarsVertTot + pRegion%specInput%nSpecies
        
        DO ivg = 1,pGrid%nVertTot
          DO iVar = 1,pRegion%specInput%nSpeciesEE
            iVarVertTot = nVarsVertTot

!            pRegion%varVertTEC(ivg,iVarVertTot+EEV_SPEC_NVAR*(iVar-1)+1) & 
!              = pRegion%spec%eevVert(EEV_SPEC_XVEL,iVar,ivg)
!            pRegion%varVertTEC(ivg,iVarVertTot+EEV_SPEC_NVAR*(iVar-1)+2) & 
!              = pRegion%spec%eevVert(EEV_SPEC_YVEL,iVar,ivg)
!            pRegion%varVertTEC(ivg,iVarVertTot+EEV_SPEC_NVAR*(iVar-1)+3) & 
!              = pRegion%spec%eevVert(EEV_SPEC_ZVEL,iVar,ivg)
!            pRegion%varVertTEC(ivg,iVarVertTot+EEV_SPEC_NVAR*(iVar-1)+4) & 
!              = pRegion%spec%eevVert(EEV_SPEC_TEMP,iVar,ivg)
            pRegion%varVertTEC(ivg,iVarVertTot+4*(iVar-1)+1) & 
              = pRegion%spec%eevVert(1,iVar,ivg)
            pRegion%varVertTEC(ivg,iVarVertTot+4*(iVar-1)+2) & 
              = pRegion%spec%eevVert(2,iVar,ivg)
            pRegion%varVertTEC(ivg,iVarVertTot+4*(iVar-1)+3) & 
              = pRegion%spec%eevVert(3,iVar,ivg)
            pRegion%varVertTEC(ivg,iVarVertTot+4*(iVar-1)+4) & 
              = pRegion%spec%eevVert(4,iVar,ivg)                                                        
          END DO ! iVar
        END DO ! ivg

!        nVarsVertTot = nVarsVertTot + pRegion%specInput%nSpecies*EEV_SPEC_NVAR
        nVarsVertTot = nVarsVertTot + pRegion%specInput%nSpecies*4                
      END IF ! global%specUsed
    ELSE 
      IF ( global%specUsed .EQV. .TRUE. ) THEN
        DO icg = 1,pGrid%nCellsTot
          DO iVar = 1,pRegion%specInput%nSpecies
            iVarCellTot = nVarsCellTot

            pRegion%varCellTEC(icg,iVarCellTot+iVar) = pRegion%spec%cv(iVar,icg)
          END DO ! iVar
        END DO ! icg

        nVarsCellTot = nVarsCellTot + pRegion%specInput%nSpecies
        
        DO icg = 1,pGrid%nCellsTot
          DO iVar = 1,pRegion%specInput%nSpeciesEE
            iVarCellTot = nVarsCellTot

!            pRegion%varCellTEC(icg,iVarCellTot+EEV_SPEC_NVAR*(iVar-1)+1) & 
!              = pRegion%spec%eev(EEV_SPEC_XVEL,iVar,icg)
!            pRegion%varCellTEC(icg,iVarCellTot+EEV_SPEC_NVAR*(iVar-1)+2) & 
!              = pRegion%spec%eev(EEV_SPEC_YVEL,iVar,icg)
!            pRegion%varCellTEC(icg,iVarCellTot+EEV_SPEC_NVAR*(iVar-1)+3) & 
!              = pRegion%spec%eev(EEV_SPEC_ZVEL,iVar,icg)
!            pRegion%varCellTEC(icg,iVarCellTot+EEV_SPEC_NVAR*(iVar-1)+4) & 
!              = pRegion%spec%eev(EEV_SPEC_TEMP,iVar,icg) 
            pRegion%varCellTEC(icg,iVarCellTot+4*(iVar-1)+1) & 
              = pRegion%spec%eev(1,iVar,icg)
            pRegion%varCellTEC(icg,iVarCellTot+4*(iVar-1)+2) & 
              = pRegion%spec%eev(2,iVar,icg)
            pRegion%varCellTEC(icg,iVarCellTot+4*(iVar-1)+3) & 
              = pRegion%spec%eev(3,iVar,icg)
            pRegion%varCellTEC(icg,iVarCellTot+4*(iVar-1)+4) & 
              = pRegion%spec%eev(4,iVar,icg)  
          END DO ! iVar
        END DO ! icg

!        nVarsCellTot = nVarsCellTot + pRegion%specInput%nSpecies*EEV_SPEC_NVAR
        nVarsCellTot = nVarsCellTot + pRegion%specInput%nSpecies*4               
      END IF ! global%specUsed
    END IF ! global%postInterpType

  END IF ! global%postPlotType

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Building volume data for TECPLOT field file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_BuildDataFieldVol







! ******************************************************************************
!
! Purpose: Collect patch data for writing to TECPLOT patch file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   pPatch              Pointer to patch data
!
! Output: None.
!
! Notes:
!   1. Isolated this code so that can easily add data for writing to file.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_BuildDataPatch(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,ifg,ifl,iPatch,iVar,iVarFaceTot,iVarVertTot,iVarTot, &
             ivg,ivl,nVarsFaceTot,nVarsVertTot
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECFIL100

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_BuildDataPatch', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Building data for TECPLOT patch file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Set TECPLOT file context to patch file
! ******************************************************************************

  errorFlag = TECFIL100(fileType2CntrTEC(FILE_TYPE_PATCH))
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,1323)
  END IF ! global%error

! ******************************************************************************
! Count number of variables
! ******************************************************************************

! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  nVarsTEC = 3

! ==============================================================================
! Mixture
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
    nVarsTEC = nVarsTEC + 5
  END IF ! global%postPlotType

! ==============================================================================
! Physical modules
! ==============================================================================

! ******************************************************************************
! Set position of variables
! ******************************************************************************

  ALLOCATE(posTEC(nVarsTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,1355,'posTEC')
  END IF ! global%error

! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  iVarTot = 0

  DO iVar = 1,3
    iVarTot = iVarTot + 1
    posTEC(iVarTot) = VAR_POS_VERT
  END DO ! iVar

! ==============================================================================
! Mixture
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
    DO iVar = 1,5
      iVarTot = iVarTot + 1
      posTEC(iVarTot) = VAR_POS_FACE
    END DO ! iVar
  END IF ! global%postPlotType
  
! ==============================================================================
! Physical modules
! ==============================================================================

! ******************************************************************************
! Determine how many variables of each position
! ******************************************************************************

  nVarsFaceTEC = 0
  nVarsVertTEC = 0

  DO iVar = 1,nVarsTEC
    IF ( posTEC(iVar) == VAR_POS_FACE ) THEN
      nVarsFaceTEC = nVarsFaceTEC + 1
    ELSE IF ( posTEC(iVar) == VAR_POS_VERT ) THEN
      nVarsVertTEC = nVarsVertTEC + 1
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,1397)
    END IF ! posTEC
  END DO ! iVar

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  ALLOCATE(pPatch%varVertTEC(pPatch%nBVert,nVarsVertTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,1408,'pPatch%varVertTEC')
  END IF ! global%error

  IF ( nVarsFaceTEC > 0 ) THEN
    ALLOCATE(pPatch%varFaceTEC(pPatch%nBFaces,nVarsFaceTEC),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,1415,'pPatch%varFaceTEC')
    END IF ! global%error
  ELSE
    NULLIFY(pPatch%varFaceTEC)
  END IF ! nVarsFaceTEC

! ******************************************************************************
! Assemble data
! ******************************************************************************

! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  DO ivl = 1,pPatch%nBVert
    ivg = pPatch%bv(ivl)

    pPatch%varVertTEC(ivl,1) = pGrid%xyz(XCOORD,ivg)
    pPatch%varVertTEC(ivl,2) = pGrid%xyz(YCOORD,ivg)
    pPatch%varVertTEC(ivl,3) = pGrid%xyz(ZCOORD,ivg)
  END DO ! ivl

  nVarsFaceTot = 0
  nVarsVertTot = 3

! ==============================================================================
! Variables
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
    nVarsFaceTot = 0

! ------------------------------------------------------------------------------
!   Mixture
! ------------------------------------------------------------------------------

    DO ifl = 1,pPatch%nBFaces
      iVarFaceTot = nVarsFaceTot

      pPatch%varFaceTEC(ifl,iVarFaceTot+1) = pPatch%cp(ifl)
      pPatch%varFaceTEC(ifl,iVarFaceTot+2) = pPatch%cf(XCOORD,ifl)
      pPatch%varFaceTEC(ifl,iVarFaceTot+3) = pPatch%cf(YCOORD,ifl)
      pPatch%varFaceTEC(ifl,iVarFaceTot+4) = pPatch%cf(ZCOORD,ifl)
      pPatch%varFaceTEC(ifl,iVarFaceTot+5) = pPatch%ch(ifl)
    END DO ! ifl

    nVarsFaceTot = nVarsFaceTot + 5

! ------------------------------------------------------------------------------
!   Physical modules
! ------------------------------------------------------------------------------

  END IF ! global%postPlotType

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Building data for TECPLOT patch file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_BuildDataPatch









! ******************************************************************************
!
! Purpose: Collect patch statistics data for writing to TECPLOT patch file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   pPatch              Pointer to patch data
!
! Output: None.
!
! Notes:
!   1. Isolated this code so that can easily add data for writing to file.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_BuildDataPatchStats(pRegion,pPatch)


  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,ifg,ifl,iPatch,iVar,iVarFaceTot,iVarVertTot,iVarTot, &
             ivg,ivl,nVarsFaceTot,nVarsVertTot
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECFIL100

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_BuildDataPatchStats', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Building data for TECPLOT patch statistics file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Set TECPLOT file context to patch statistics file
! ******************************************************************************

  errorFlag = TECFIL100(fileType2CntrTEC(FILE_TYPE_PATCH_STATS))
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,1572)
  END IF ! global%error

! ******************************************************************************
! Count number of variables
! ******************************************************************************

! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  nVarsTEC = 3

! ==============================================================================
! Mixture
! ==============================================================================

! ==============================================================================
! Physical modules
! ==============================================================================

! ------------------------------------------------------------------------------
! Lagrangian particle statistics
! ------------------------------------------------------------------------------


! ******************************************************************************
! Set position of variables
! ******************************************************************************

  ALLOCATE(posTEC(nVarsTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,1614,'posTEC')
  END IF ! global%error

! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  iVarTot = 0

  DO iVar = 1,3
    iVarTot = iVarTot + 1
    posTEC(iVarTot) = VAR_POS_VERT
  END DO ! iVar

! ==============================================================================
! Mixture
! ==============================================================================
  
! ==============================================================================
! Physical modules
! ==============================================================================

! ------------------------------------------------------------------------------
! Lagrangian particle statistics
! ------------------------------------------------------------------------------

  
! ******************************************************************************
! Determine how many variables of each position
! ******************************************************************************

  nVarsFaceTEC = 0
  nVarsVertTEC = 0

  DO iVar = 1,nVarsTEC
    IF ( posTEC(iVar) == VAR_POS_FACE ) THEN
      nVarsFaceTEC = nVarsFaceTEC + 1
    ELSE IF ( posTEC(iVar) == VAR_POS_VERT ) THEN
      nVarsVertTEC = nVarsVertTEC + 1
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,1662)
    END IF ! posTEC
  END DO ! iVar

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  ALLOCATE(pPatch%varVertTEC(pPatch%nBVert,nVarsVertTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,1673,'pPatch%varVertTEC')
  END IF ! global%error

  IF ( nVarsFaceTEC > 0 ) THEN
    ALLOCATE(pPatch%varFaceTEC(pPatch%nBFaces,nVarsFaceTEC),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,1680,'pPatch%varFaceTEC')
    END IF ! global%error
  ELSE
    NULLIFY(pPatch%varFaceTEC)
  END IF ! nVarsFaceTEC

! ******************************************************************************
! Assemble data
! ******************************************************************************

! ==============================================================================
! Coordinates. NOTE do not add any data before coordinates!
! ==============================================================================

  DO ivl = 1,pPatch%nBVert
    ivg = pPatch%bv(ivl)

    pPatch%varVertTEC(ivl,1) = pGrid%xyz(XCOORD,ivg)
    pPatch%varVertTEC(ivl,2) = pGrid%xyz(YCOORD,ivg)
    pPatch%varVertTEC(ivl,3) = pGrid%xyz(ZCOORD,ivg)
  END DO ! ivl

  nVarsFaceTot = 0
  nVarsVertTot = 3

! ==============================================================================
! Variables
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
    nVarsFaceTot = 0

! ------------------------------------------------------------------------------
!   Mixture
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   Physical modules
! ------------------------------------------------------------------------------

! - Lagrangian particle statistics ---------------------------------------------

  END IF ! global%postPlotType

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Building data for TECPLOT patch statistics file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_BuildDataPatchStats







! ******************************************************************************
!
! Purpose: Build header for field data.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_BuildHeaderField(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iPv,iPv2,iVar
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_BuildHeaderField', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Building header for TECPLOT field file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Build header
! ******************************************************************************

! ==============================================================================
! Write grid only
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_ONLY ) THEN
    WRITE(headerTEC,'(3(1X,A))') 'x','y','z'
  ELSE
  
! ==============================================================================
! Write grid and solution depending on fluid model
! ==============================================================================
  
    SELECT CASE ( pRegion%mixtInput%fluidModel ) 
    
! ------------------------------------------------------------------------------
!     Incompressible fluid model
! ------------------------------------------------------------------------------    
    
      CASE ( FLUID_MODEL_INCOMP ) 
        WRITE(headerTEC,'(3(1X,A))') 'x','y','z'

        IF ( global%postPlotMixtCvFlag .EQV. .TRUE. ) THEN 
          WRITE(headerTEC,'(A,4(1X,A))') TRIM(headerTEC),'u','v','w','p'        
        END IF ! global%postPlotMixtCvFlag
      
! ------------------------------------------------------------------------------
!     Compressible fluid model
! ------------------------------------------------------------------------------      
      
      CASE ( FLUID_MODEL_COMP ) 
        WRITE(headerTEC,'(3(1X,A))') 'x','y','z'
        
        IF ( global%postPlotMixtCvFlag .EQV. .TRUE. ) THEN 
          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            WRITE(headerTEC,'(A,5(1X,A))') TRIM(headerTEC),'r','u','v','w','p'
          ELSE
            WRITE(headerTEC,'(A,5(1X,A))') TRIM(headerTEC),'r','ru','rv', &
                                                           'rw','rE'
          END IF ! global%solverType
        END IF ! global%postPlotMixtCvFlag

        IF ( global%postPlotMixtDvFlag .EQV. .TRUE. ) THEN 
          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            IF ( pRegion%plot%nDvMixt > 0 ) THEN
              WRITE(headerTEC,'(A,1(1X,A))') TRIM(headerTEC),'T'
            END IF ! pRegion%plot%nDvMixt
          ELSE
            IF ( pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL ) THEN
              WRITE(headerTEC,'(A,5(1X,A))') &
                                   TRIM(headerTEC),'p','T','a','eA','eB'
            ELSE
              WRITE(headerTEC,'(A,3(1X,A))') TRIM(headerTEC),'p','T','a'
            END IF ! pRegion%mixtInput%gasModel            
          END IF ! global%solverType
        END IF ! global%postPlotMixtDvFlag

        DO iVar = 1,pRegion%plot%nGvMixt
          WRITE(headerTEC,'(A,1X,A,I2.2)') TRIM(headerTEC),'GV',iVar
        END DO ! iVar  

        DO iPv = 1,pRegion%plot%nPv
          iPv2 = pRegion%plot%pvi2pv(iPv)
          
          WRITE(headerTEC,'(A,1X,A)') TRIM(headerTEC), &
                                      TRIM(pRegion%plot%pvNameShort(iPv2))
        END DO ! iPv   

        DO iVar = 1,pRegion%specInput%nSpecies
          WRITE(headerTEC,'(A,1X,A,I2.2)') TRIM(headerTEC),'rY',iVar
        END DO ! iVar
        
        DO iVar = 1,pRegion%specInput%nSpeciesEE
          WRITE(headerTEC,'(A,1X,A,I2.2)') TRIM(headerTEC),'u',iVar
          WRITE(headerTEC,'(A,1X,A,I2.2)') TRIM(headerTEC),'v',iVar
          WRITE(headerTEC,'(A,1X,A,I2.2)') TRIM(headerTEC),'w',iVar
          WRITE(headerTEC,'(A,1X,A,I2.2)') TRIM(headerTEC),'T',iVar                                      
        END DO ! iVar

! ------------------------------------------------------------------------------
!     Default
! ------------------------------------------------------------------------------

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,1894)
    END SELECT ! pRegion%mixtInput%fluidModel 
  END IF ! global%postPlotType

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Building header for TECPLOT field file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_BuildHeaderField






! ******************************************************************************
!
! Purpose: Build header for patch data.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_BuildHeaderPatch(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iVar
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_BuildHeaderPatch', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Building header for TECPLOT patch file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Build header
! ******************************************************************************

! ==============================================================================
! Mixture 
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_ONLY ) THEN
    WRITE(headerTEC,'(3(1X,A))') 'x','y','z'
  ELSE
    WRITE(headerTEC,'(8(1X,A))') 'x','y','z', &
                                 'Cp','Cfx','Cfy','Cfz','Ch'
  END IF ! global%postPlotType

! ==============================================================================
! Physical modules
! ==============================================================================

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Building header for TECPLOT patch file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_BuildHeaderPatch








! ******************************************************************************
!
! Purpose: Build header for patch statistics data.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_BuildHeaderPatchStats(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iVar
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_BuildHeaderPatchStats', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Building header for TECPLOT patch statistics file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Build header
! ******************************************************************************

! ==============================================================================
! Mixture 
! ==============================================================================

  WRITE(headerTEC,'(3(1X,A))') 'x','y','z'

! ==============================================================================
! Physical modules
! ==============================================================================

! ------------------------------------------------------------------------------
! Lagrangian particle statistics
! ------------------------------------------------------------------------------


! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Building header for TECPLOT patch statistics file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_BuildHeaderPatchStats








! ******************************************************************************
!
! Purpose: Close TECPLOT field file.
!
! Description: None.
!
! Input:
!   global                Pointer to global type
!
! Output: None.
!
! Notes:
!   1. Isolated into separate routine so that can write several zones to same
!      TECPLOT file.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_CloseFileField(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECFIL100

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_TEC_CloseFileField', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing TECPLOT field file...'
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Be patient, this may take a while...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set TECPLOT file context to field file
! ******************************************************************************

  errorFlag = TECFIL100(fileType2CntrTEC(FILE_TYPE_FIELD))
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,2169)
  END IF ! global%error

! ******************************************************************************
! Close file
! ******************************************************************************

  CALL RFLU_TEC_CloseFile(global)

  fileCntrTEC = fileCntrTEC - 1

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing TECPLOT field file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_CloseFileField







! ******************************************************************************
!
! Purpose: Close TECPLOT patch file.
!
! Description: None.
!
! Input:
!   global                Pointer to global type
!
! Output: None.
!
! Notes:
!   1. Isolated into separate routine so that can write several zones to same
!      TECPLOT file.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_CloseFilePatch(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECFIL100

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_TEC_CloseFilePatch', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing TECPLOT patch file...'
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Be patient, this may take a while...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set TECPLOT file context to field file
! ******************************************************************************

  errorFlag = TECFIL100(fileType2CntrTEC(FILE_TYPE_PATCH))
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,2258)
  END IF ! global%error

! ******************************************************************************
! Close file
! ******************************************************************************

  CALL RFLU_TEC_CloseFile(global)

  fileCntrTEC = fileCntrTEC - 1

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing TECPLOT patch file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_CloseFilePatch






! ******************************************************************************
!
! Purpose: Close TECPLOT patch statistics file.
!
! Description: None.
!
! Input:
!   global                Pointer to global type
!
! Output: None.
!
! Notes:
!   1. Isolated into separate routine so that can write several zones to same
!      TECPLOT file.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_CloseFilePatchStats(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECFIL100

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_TEC_CloseFilePatchStats', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Closing TECPLOT patch statistics file...'
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Be patient, this may take a while...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set TECPLOT file context to field file
! ******************************************************************************

  errorFlag = TECFIL100(fileType2CntrTEC(FILE_TYPE_PATCH_STATS))
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,2347)
  END IF ! global%error

! ******************************************************************************
! Close file
! ******************************************************************************

  CALL RFLU_TEC_CloseFile(global)

  fileCntrTEC = fileCntrTEC - 1

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Closing TECPLOT patch statistics file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_CloseFilePatchStats







! ******************************************************************************
!
! Purpose: Close TECPLOT point file.
!
! Description: None.
!
! Input:
!   global                Pointer to global type
!
! Output: None.
!
! Notes:
!   1. Isolated into separate routine so that can write several zones to same
!      TECPLOT file.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_CloseFilePnt(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_TEC_CloseFilePnt', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing TECPLOT point file...'
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Be patient, this may take a while...'
  END IF ! global%verbLevel

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_PLOT,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,2431)
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing TECPLOT point file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_CloseFilePnt








! ******************************************************************************
!
! Purpose: Destroy surface field data for writing to TECPLOT file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   pPatch              Pointer to patch data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_DestroyDataFieldSurf(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch  
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_DestroyDataFieldSurf', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Destroying TECPLOT surface field data...'
  END IF ! global%verbLevel

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(posTEC,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,2512,'posTEC')
  END IF ! global%error

  DEALLOCATE(pPatch%varVertTEC,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,2518,'pPatch%varVertTEC')
  END IF ! global%error

  IF ( ASSOCIATED(pPatch%varFaceTEC) .EQV. .TRUE. ) THEN
    DEALLOCATE(pPatch%varFaceTEC,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,2525,'pPatch%varFaceTEC')
    END IF ! global%error
  END IF ! ASSOCIATED

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Destroying TECPLOT surface field data done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_DestroyDataFieldSurf







! ******************************************************************************
!
! Purpose: Destroy volume field data for writing to TECPLOT file.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_DestroyDataFieldVol(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_DestroyDataFieldVol', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Destroying TECPLOT volume field data...'
  END IF ! global%verbLevel

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(posTEC,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,2605,'posTEC')
  END IF ! global%error

  DEALLOCATE(pRegion%varVertTEC,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,2611,'pRegion%varVertTEC')
  END IF ! global%error

  IF ( ASSOCIATED(pRegion%varCellTEC) .EQV. .TRUE. ) THEN
    DEALLOCATE(pRegion%varCellTEC,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,2618,'pRegion%varCellTEC')
    END IF ! global%error
  END IF ! ASSOCIATED

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Destroying TECPLOT volume field data done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_DestroyDataFieldVol







! ******************************************************************************
!
! Purpose: Destroy patch data for writing to TECPLOT file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   pPatch              Pointer to patch data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_DestroyDataPatch(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,iPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_DestroyDataPatch', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying TECPLOT patch data...'
  END IF ! global%verbLevel

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(posTEC,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,2699,'posTEC')
  END IF ! global%error

  DEALLOCATE(pPatch%varVertTEC,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,2705,'pPatch%varVertTEC')
  END IF ! global%error

  IF ( ASSOCIATED(pPatch%varFaceTEC) .EQV. .TRUE. ) THEN
    DEALLOCATE(pPatch%varFaceTEC,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,2712,'pPatch%varFaceTEC')
    END IF ! global%error
  END IF ! ASSOCIATED

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Destroying TECPLOT patch data done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_DestroyDataPatch






! ******************************************************************************
!
! Purpose: Initialize TECPLOT interface
!
! Description: None.
!
! Input:
!   global                Pointer to global data
!
! Output: None.
!
! Notes:
!   1. Needed because TECPLOT interface has some global data which needs to be
!      initialized once.
!   2. Must be called before first TECPLOT file opened.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_Init(global)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_TEC_Init', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing TECPLOT...'
  END IF ! global%verbLevel

! ******************************************************************************
! Initialize data
! ******************************************************************************

  fileCntrTEC = 0

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing TECPLOT done...'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_Init






! ******************************************************************************
!
! Purpose: Open TECPLOT field file
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Isolated into separate routine so that can write several zones to same
!      TECPLOT file.
!   2. Need pointer to region because need to know how many variables are
!      to be written to file, and hence need to know about number of species
!      (as an example). Passing the global does not allow for this kind of
!      information because the user input is under the region data structure.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_OpenFileField(pRegion)

  USE ModBuildFileNames, ONLY: BuildFileNamePlainSteady, &
                               BuildFileNamePlainUnsteady

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,title
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_OpenFileField', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening TECPLOT field file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Build header
! ******************************************************************************

  CALL RFLU_TEC_BuildHeaderField(pRegion)

! ******************************************************************************
! Open file
! ******************************************************************************

  title = global%caseName

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNamePlainUnsteady(global,FILEDEST_INDIR,'.plt', &
                                    global%currentTime,iFileName)
  ELSE IF ( global%flowType == FLOW_STEADY ) THEN
    CALL BuildFileNamePlainSteady(global,FILEDEST_INDIR,'.plt', &
                                  global%currentIter,iFileName)
  ELSE ! defensive coding
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,2878)
  END IF ! global%flowType

  fileCntrTEC = fileCntrTEC + 1

  IF ( fileCntrTEC > FILE_CNTR_TEC_MAX ) THEN
    CALL ErrorStop(global,ERR_TECPLOT_FILECNTR,2884)
  END IF ! fileCntrTEC

  CALL RFLU_TEC_OpenFile(global,title,iFileName)

  fileType2CntrTEC(FILE_TYPE_FIELD) = fileCntrTEC

! ==============================================================================
! Write info
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'File name:',TRIM(iFileName)
    WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'File counter:',fileCntrTEC
  END IF ! global%verbLevel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening TECPLOT field file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_OpenFileField






! ******************************************************************************
!
! Purpose: Open TECPLOT patch file
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Isolated into separate routine so that can write several zones to same
!      TECPLOT file.
!   2. Need pointer to region because need to know how many variables are
!      to be written to file, and hence need to know about number of species
!      (as an example). Passing the global does not allow for this kind of
!      information because the user input is under the region data structure.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_OpenFilePatch(pRegion)

  USE ModBuildFileNames, ONLY: BuildFileNamePlainSteady, &
                               BuildFileNamePlainUnsteady

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,title
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_OpenFilePatch', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening TECPLOT patch file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Build header
! ******************************************************************************

  CALL RFLU_TEC_BuildHeaderPatch(pRegion)

! ******************************************************************************
! Open file
! ******************************************************************************

  title = global%caseName

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNamePlainUnsteady(global,FILEDEST_INDIR,'.pat.plt', &
                                    global%currentTime,iFileName)
  ELSE IF ( global%flowType == FLOW_STEADY ) THEN
    CALL BuildFileNamePlainSteady(global,FILEDEST_INDIR,'.pat.plt', &
                                  global%currentIter,iFileName)
  ELSE ! defensive coding
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,2995)
  END IF ! global%flowType

  fileCntrTEC = fileCntrTEC + 1

  IF ( fileCntrTEC > FILE_CNTR_TEC_MAX ) THEN
    CALL ErrorStop(global,ERR_TECPLOT_FILECNTR,3001)
  END IF ! fileCntrTEC

  CALL RFLU_TEC_OpenFile(global,title,iFileName)

  fileType2CntrTEC(FILE_TYPE_PATCH) = fileCntrTEC

! ==============================================================================
! Write info
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'File name:',TRIM(iFileName)
    WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'File counter:',fileCntrTEC
  END IF ! global%verbLevel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening TECPLOT patch file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_OpenFilePatch








! ******************************************************************************
!
! Purpose: Open TECPLOT patch statistics file
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Isolated into separate routine so that can write several zones to same
!      TECPLOT file.
!   2. Need pointer to region because need to know how many variables are
!      to be written to file, and hence need to know about number of species
!      (as an example). Passing the global does not allow for this kind of
!      information because the user input is under the region data structure.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_OpenFilePatchStats(pRegion)

  USE ModBuildFileNames, ONLY: BuildFileNamePlainSteady, &
                               BuildFileNamePlainUnsteady

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,title
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_OpenFilePatchStats', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Opening TECPLOT patch statistics file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Build header
! ******************************************************************************

  CALL RFLU_TEC_BuildHeaderPatchStats(pRegion)

! ******************************************************************************
! Open file
! ******************************************************************************

  title = global%caseName

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNamePlainUnsteady(global,FILEDEST_INDIR,'.stats.plt', &
                                    global%currentTime,iFileName)
  ELSE IF ( global%flowType == FLOW_STEADY ) THEN
    CALL BuildFileNamePlainSteady(global,FILEDEST_INDIR,'.stats.plt', &
                                  global%currentIter,iFileName)
  ELSE ! defensive coding
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,3115)
  END IF ! global%flowType

  fileCntrTEC = fileCntrTEC + 1

  IF ( fileCntrTEC > FILE_CNTR_TEC_MAX ) THEN
    CALL ErrorStop(global,ERR_TECPLOT_FILECNTR,3121)
  END IF ! fileCntrTEC

  CALL RFLU_TEC_OpenFile(global,title,iFileName)

  fileType2CntrTEC(FILE_TYPE_PATCH_STATS) = fileCntrTEC

! ==============================================================================
! Write info
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'File name:',TRIM(iFileName)
    WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'File counter:',fileCntrTEC
  END IF ! global%verbLevel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Opening TECPLOT patch statistics file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_OpenFilePatchStats






! ******************************************************************************
!
! Purpose: Open TECPLOT point file
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Isolated into separate routine so that can write several zones to same
!      TECPLOT file.
!   2. Need pointer to region because need to know how many variables are
!      to be written to file, and hence need to know about number of species
!      (as an example). Passing the global does not allow for this kind of
!      information because the user input is under the region data structure.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_OpenFilePnt(pRegion)

  USE ModBuildFileNames, ONLY: BuildFileNamePlainUnsteady

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName
  INTEGER :: errorFlag,iFile
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_OpenFilePnt', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening TECPLOT point file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Open file
! ******************************************************************************

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNamePlainUnsteady(global,FILEDEST_INDIR,'.plag.dat', &
                                    global%currentTime,iFileName)
  ELSE ! defensive coding
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,3221)
  END IF ! global%flowType

  iFile = IF_PLOT
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,3229,TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! Write header
! ******************************************************************************

  WRITE(iFile,'(1X,A)') 'TITLE="'//TRIM(global%casename)//'"' 
  WRITE(iFile,'(1X,A)') 'VARIABLES="x" "y" "z" "up" "vp" "wp" "diam"'

! ******************************************************************************
! Write info
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'File name:',TRIM(iFileName)
  END IF ! global%verbLevel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening TECPLOT point file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_OpenFilePnt







! ******************************************************************************
!
! Purpose: Write surface data to TECPLOT field file.
!
! Description: None.
!
! Input:
!   pRegion        Pointer to region data
!   pPatch         Pointer to patch data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_WriteFileFieldSurf(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_WriteFileFieldSurf', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Writing surface data to TECPLOT file...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal
    WRITE(STDOUT,'(A,3X,A,1X,I5)') SOLVER_NAME,'Global patch:', &
                                   pPatch%iPatchGlobal                               
  END IF ! global%verbLevel

! ******************************************************************************
! Triangles
! ******************************************************************************

  IF ( pPatch%nBTris > 0 ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Actual triangles...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneSurf(pRegion,pPatch,FACE_TYPE_TRI,FACE_KIND_AB)
  ENDIF ! pPatch%nBTris

  IF ( (global%postZoneVirtFlag .EQV. .TRUE.) .AND. & 
       (pPatch%nBTrisTot > pPatch%nBTris) ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Virtual triangles...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneSurf(pRegion,pPatch,FACE_TYPE_TRI,FACE_KIND_VB)
  ENDIF ! global%postZoneVirtFlag

! ******************************************************************************
!   Quadrilaterals
! ******************************************************************************

  IF ( pPatch%nBQuads > 0 ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Actual quadrilaterals...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneSurf(pRegion,pPatch,FACE_TYPE_QUAD,FACE_KIND_AB)
  END IF ! pPatch%nBQuads

  IF ( (global%postZoneVirtFlag .EQV. .TRUE.) .AND. & 
       (pPatch%nBQuadsTot > pPatch%nBQuads) ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Virtual quadrilaterals...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneSurf(pRegion,pPatch,FACE_TYPE_QUAD,FACE_KIND_VB)
  END IF ! global%postZoneVirtFlag

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Writing surface data to TECPLOT field file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_WriteFileFieldSurf







! ******************************************************************************
!
! Purpose: Write volume data to TECPLOT field file
!
! Description: None.
!
! Input:
!   pRegion        Pointer to region data
!
! Output: None.
!
! Notes:
!   1. Write partition boundary faces here because they access field volume
!      data.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_WriteFileFieldVol(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_WriteFileFieldVol', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Writing volume data to TECPLOT field file...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointer to grid
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Tetrahedra
! ******************************************************************************

  IF ( pGrid%nTets > 0 ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Actual tetrahedra...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneVol(pRegion,CELL_TYPE_TET,CELL_KIND_ACTUAL)
  ENDIF ! pGrid%nTets

  IF ( (global%postZoneVirtFlag .EQV. .TRUE.) .AND. & 
       (pGrid%nTetsTot > pGrid%nTets) ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Virtual tetrahedra...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneVol(pRegion,CELL_TYPE_TET,CELL_KIND_VIRTUAL)
  ENDIF ! global%postZoneVirtFlag

! ******************************************************************************
! Hexahedra
! ******************************************************************************

  IF ( pGrid%nHexs > 0 ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Actual hexahedra...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneVol(pRegion,CELL_TYPE_HEX,CELL_KIND_ACTUAL)
  ENDIF ! pGrid%nHexs

  IF ( (global%postZoneVirtFlag .EQV. .TRUE.) .AND. & 
       (pGrid%nHexsTot > pGrid%nHexs) ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Virtual hexahedra...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneVol(pRegion,CELL_TYPE_HEX,CELL_KIND_VIRTUAL)
  ENDIF ! global%postZoneVirtFlag

! ******************************************************************************
! Prisms
! ******************************************************************************

  IF ( pGrid%nPris > 0 ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Actual prisms...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneVol(pRegion,CELL_TYPE_PRI,CELL_KIND_ACTUAL)
  ENDIF ! pGrid%nPris

  IF ( (global%postZoneVirtFlag .EQV. .TRUE.) .AND. & 
       (pGrid%nPrisTot > pGrid%nPris) ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Virtual prisms...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneVol(pRegion,CELL_TYPE_PRI,CELL_KIND_VIRTUAL)
  ENDIF ! global%postZoneVirtFlag

! ******************************************************************************
! Pyramids
! ******************************************************************************

  IF ( pGrid%nPyrs > 0 ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Actual pyramids...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneVol(pRegion,CELL_TYPE_PYR,CELL_KIND_ACTUAL)
  ENDIF ! pGrid%nPyrs

  IF ( (global%postZoneVirtFlag .EQV. .TRUE.) .AND. & 
       (pGrid%nPyrsTot > pGrid%nPyrs) ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Virtual pyramids...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneVol(pRegion,CELL_TYPE_PYR,CELL_KIND_VIRTUAL)
  ENDIF ! global%postZoneVirtFlag

! ******************************************************************************
! Writing special cells and faces
! ******************************************************************************

  IF ( pGrid%nCellsSpecial > 0 ) THEN
    IF ( global%verbLevel > VERBOSE_NONE  ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Special cells...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneCellsSpecial(pRegion)
  END IF ! pGrid%nCellsSpecial
  
  IF ( pGrid%nFacesSpecial > 0 ) THEN
    IF ( global%verbLevel > VERBOSE_NONE  ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Special faces...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneFacesSpecial(pRegion)
  END IF ! pGrid%nFacesSpecial  

! ******************************************************************************
! Partition boundaries. NOTE can only write these faces if have no cell data
! because cannot easily access cell data for these faces. 
! ******************************************************************************

  IF ( (pGrid%nFaces /= pGrid%nFacesTot) .AND. (nVarsCellTEC == 0) ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Interpartition faces...'
    END IF ! global%verbLevel

    CALL RFLU_TEC_WriteZoneInterf(pRegion)
  END IF ! pGrid%nFaces

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Writing volume data to TECPLOT field file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_WriteFileFieldVol







! ******************************************************************************
!
! Purpose: Write data to TECPLOT patch file.
!
! Description: None.
!
! Input:
!   pRegion        Pointer to region data
!   pPatch         Pointer to patch data
!
! Output: None.
!
! Notes:
!   1. Only write data for actual faces because patch data only stored for
!      actual faces.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_WriteFilePatch(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_WriteFilePatch', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing data to TECPLOT patch file...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal
    WRITE(STDOUT,'(A,3X,A,1X,I5)') SOLVER_NAME,'Global patch:', &
                                   pPatch%iPatchGlobal                               
  END IF ! global%verbLevel

! ******************************************************************************
! Loop over patches
! ******************************************************************************

  IF ( pPatch%nBFaces > 0 ) THEN
    CALL RFLU_TEC_WriteZoneSurfMixed(pRegion,pPatch)
  END IF ! pPatch%nBFaces

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Writing data to TECPLOT patch file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_WriteFilePatch







! ******************************************************************************
!
! Purpose: Write point data to TECPLOT file
!
! Description: None.
!
! Input:
!   pRegion        Pointer to region data
!
! Output: None.
!
! Notes:
!   1. At present hard-wired for Rocpart. Will be extended eventually to allow
!      arbitrary point data to be written.
!
! ******************************************************************************

SUBROUTINE RFLU_TEC_WriteFilePnt(pRegion)


  USE ModTools, ONLY: FloatEqual


  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,iPcl
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_WriteFilePnt', &
                        "../../utilities/post/RFLU_ModTECPLOT.F90")


! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing to TECPLOT point file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_WriteFilePnt







END MODULE RFLU_ModTECPLOT

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModTECPLOT.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.8  2008/12/06 08:43:58  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:12  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2008/01/07 14:04:42  haselbac
! Added IFs on postZoneVirtFlag
!
! Revision 1.5  2007/11/30 13:52:53  mparmar
! Removed temporary DEBUG statement
!
! Revision 1.4  2007/11/28 23:05:53  mparmar
! Added logic to add temperature in tecplot file for SOLV_IMPLICIT_HM
!
! Revision 1.3  2007/09/04 13:21:34  haselbac
! Added capability to include subset of variables in plt files
!
! Revision 1.2  2007/04/12 17:58:07  haselbac
! Added capability to write only a subset of particles to file
!
! Revision 1.1  2007/04/09 18:58:09  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.36  2007/03/19 21:44:24  haselbac
! Adapted to changes related to plotting variables
!
! Revision 1.35  2006/05/02 18:23:33  fnajjar
! Reverting to original code for correctness
!
! Revision 1.34  2006/05/02 17:48:09  fnajjar
! Allowed surface statistics to be written on active patches
!
! Revision 1.33  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.32  2006/03/26 20:22:34  haselbac
! Removed error trap for GL model
!
! Revision 1.31  2005/12/13 23:30:53  haselbac
! Cosmetics
!
! Revision 1.30  2005/12/13 23:11:46  fnajjar
! Added if statement to write Tecplot header for non-null nPcls
!
! Revision 1.29  2005/11/27 02:01:53  haselbac
! Added support for EEv
!
! Revision 1.28  2005/11/17 22:33:21  haselbac
! Bug fixes to allow only geometry to be postprocd with 1
!
! Revision 1.27  2005/11/10 02:52:23  haselbac
! Added writing of gv for variable properties cases
!
! Revision 1.26  2005/09/23 19:01:31  haselbac
! Added capability to write patch stats files
!
! Revision 1.25  2005/08/10 00:38:48  haselbac
! Modified writing of PV labels into header
!
! Revision 1.24  2005/08/09 01:11:47  haselbac
! Rewrote field surf routines to operate on one patch at a time
!
! Revision 1.23  2005/05/18 22:24:57  fnajjar
! ACH: Adapted point files to multiple regions
!
! Revision 1.22  2005/05/01 14:23:10  haselbac
! Added processing of plotting vars
!
! Revision 1.21  2005/01/06 04:42:51  haselbac
! Now write partition boundaries only if have no cell data
!
! Revision 1.20  2004/12/27 23:34:41  haselbac
! Added writing of field cell and face data
!
! Revision 1.19  2004/12/21 15:09:54  fnajjar
! Added PLAG surface statistics to Tecplot file
!
! Revision 1.18  2004/11/14 19:59:31  haselbac
! Added code for incompressible fluid model
!
! Revision 1.17  2004/09/27 01:42:24  haselbac
! Added call to write zone with special faces
!
! Revision 1.16  2004/07/28 15:29:21  jferry
! created global variable for spec use
!
! Revision 1.15  2004/07/20 03:11:58  haselbac
! Bug fix when writing surface data with species
!
! Revision 1.14  2004/07/07 01:00:27  haselbac
! Added NULLIFY statements to fix problems on Blue Pacific
!
! Revision 1.13  2004/07/02 03:05:56  haselbac
! Added message for closing patch file
!
! Revision 1.12  2004/06/16 20:01:43  haselbac
! Complete rewrite of module, adaptation to TEC10 and surface data
!
! Revision 1.11  2004/03/05 23:23:23  haselbac
! Added RFLU_WriteFileTECPLOTPoint, removed superfluous dupList
!
! Revision 1.10  2003/11/25 21:03:59  haselbac
! Extended module to deal with arbitrary number of variables
!
! Revision 1.9  2003/08/07 15:36:39  haselbac
! Changed var names
!
! Revision 1.8  2003/05/05 18:42:51  haselbac
! Replaced plotTypeTEC by global%plotType
!
! Revision 1.7  2003/05/02 21:44:50  haselbac
! Fixed bug: pGrid not set in surface data routine
!
! Revision 1.6  2003/05/02 16:41:12  haselbac
! Added check for zero boundary vertices for surface data file
!
! Revision 1.5  2003/04/28 22:46:04  haselbac
! Added function to write out surface data only
!
! Revision 1.4  2003/04/01 16:41:34  haselbac
! Removed getting of special cells
!
! Revision 1.3  2003/03/25 19:18:44  haselbac
! Fixed bug, now case with single face (only AV face) works
!
! Revision 1.2  2003/03/20 20:07:19  haselbac
! Modified RegFun call to avoid probs with
! long "../../utilities/post/RFLU_ModTECPLOT.F90" names
!
! Revision 1.1  2003/03/15 19:16:54  haselbac
! Initial revision
!
! ******************************************************************************

