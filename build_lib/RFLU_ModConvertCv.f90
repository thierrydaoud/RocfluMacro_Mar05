










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
! Purpose: Suite of routines to convert state vector.
!
! Description: None.
!
! Notes: 
!   1. The routines which convert scalar conserved vectors assume that that 
!      form contains variables per unit volume (i.e., density times variable
!      per unit mass). 
!
! ******************************************************************************
!
! $Id: RFLU_ModConvertCv.F90,v 1.3 2016/02/08 20:02:15 fred Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModConvertCv

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
    

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_ConvertCvCons2Prim, & 
            RFLU_ConvertCvPrim2Cons, & 
            RFLU_ScalarConvertCvCons2Prim, &
            RFLU_ScalarConvertCvPrim2Cons
   
  SAVE    
     
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS


  
! ******************************************************************************
!
! Purpose: Convert conserved state vector to primitive variables.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to data of current region
!   cvStateFuture       Future state of conserved variables
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ConvertCvCons2Prim(pRegion,cvStateFuture)

    USE ModInterfaces, ONLY: MixtPerf_R_M, &
                             MixtPerf_T_DPR, &
                             MixtPerf_G_CpR

   USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
   USE ModSpecies,           ONLY: t_spec_input

   USE RFLU_ModJWL !Fred - allowing for JWL EOS use in conversion
! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: cvStateFuture
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,indMol
    REAL(RFREAL) :: gc,ir,mw,p,r,mwAir,cpAir,gcAir,gAir,T,a,b
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global


     LOGICAL :: scalarConvFlag
     INTEGER :: iCvSpecProducts,iCvSpecAir
     REAL(RFREAL) :: Yproducts
     REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
     TYPE(t_spec_input), POINTER :: pSpecInput


! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvertCvCons2Prim',"../modflu/RFLU_ModConvertCv.F90")


! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid
    pCv   => pRegion%mixt%cv
    pDv   => pRegion%mixt%dv
    pGv   => pRegion%mixt%gv

    indMol = pRegion%mixtInput%indMol


    pSpecInput => pRegion%specInput
    pCvSpec => pRegion%spec%cv


! ******************************************************************************
!   Actual conversion
! ******************************************************************************

    SELECT CASE ( pRegion%mixt%cvState )

! ==============================================================================
!     Convert from conservative to primitive form 
! ==============================================================================

      CASE ( CV_MIXT_STATE_CONS ) 
        SELECT CASE ( cvStateFuture )

! ------------------------------------------------------------------------------      
!         Convert to duvwp form        
! ------------------------------------------------------------------------------

          CASE ( CV_MIXT_STATE_DUVWP )
            pRegion%mixt%cvState = CV_MIXT_STATE_DUVWP          
          
            DO icg = 1,pGrid%nCellsTot
              ir = 1.0_RFREAL/pCv(CV_MIXT_DENS,icg)

              pCv(CV_MIXT_XVEL,icg) = ir*pCv(CV_MIXT_XMOM,icg)
              pCv(CV_MIXT_YVEL,icg) = ir*pCv(CV_MIXT_YMOM,icg)
              pCv(CV_MIXT_ZVEL,icg) = ir*pCv(CV_MIXT_ZMOM,icg)
              pCv(CV_MIXT_PRES,icg) = pDv(DV_MIXT_PRES,icg)
            END DO ! icg

! ------------------------------------------------------------------------------          
!         Convert to duvwt form             
! ------------------------------------------------------------------------------

          CASE (CV_MIXT_STATE_DUVWT)
            pRegion%mixt%cvState = CV_MIXT_STATE_DUVWT

            SELECT CASE ( pRegion%mixtInput%fluidModel ) 

! ----------- Compressible fluid model -----------------------------------------

              CASE ( FLUID_MODEL_COMP )              
                SELECT CASE ( pRegion%mixtInput%gasModel )
                
! --------------- TC perfect gas or mixture thereof, pseudo-gas                 
                
                  CASE ( GAS_MODEL_TCPERF, &
                         GAS_MODEL_MIXT_TCPERF, & 
                         GAS_MODEL_MIXT_PSEUDO )
         
                    DO icg = 1,pGrid%nCellsTot                
                      r  = pCv(CV_MIXT_DENS,icg)
                      p  = pDv(DV_MIXT_PRES,icg)                
                      ir = 1.0_RFREAL/r

                      pCv(CV_MIXT_XVEL,icg) = ir*pCv(CV_MIXT_XMOM,icg)
                      pCv(CV_MIXT_YVEL,icg) = ir*pCv(CV_MIXT_YMOM,icg)
                      pCv(CV_MIXT_ZVEL,icg) = ir*pCv(CV_MIXT_ZMOM,icg)

                      mw = pGv(GV_MIXT_MOL,indMol*icg)
                      gc = MixtPerf_R_M(mw)


                      IF ( global%piclUsed .EQV. .TRUE. ) THEN
                        r  = r/(1.0_RFREAL - pRegion%mixt%piclVF(icg))
                      END IF

                      pCv(CV_MIXT_TEMP,icg) = MixtPerf_T_DPR(r,p,gc)
                    END DO ! icg 
                    
! --------------- Gas-liquid mixture                    
                    
                  CASE ( GAS_MODEL_MIXT_GASLIQ ) 
                    DO icg = 1,pGrid%nCellsTot
                      ir = 1.0_RFREAL/pCv(CV_MIXT_DENS,icg)

                      pCv(CV_MIXT_XVEL,icg) = ir*pCv(CV_MIXT_XMOM,icg)
                      pCv(CV_MIXT_YVEL,icg) = ir*pCv(CV_MIXT_YMOM,icg)
                      pCv(CV_MIXT_ZVEL,icg) = ir*pCv(CV_MIXT_ZMOM,icg)

                      pCv(CV_MIXT_TEMP,icg) = pDv(DV_MIXT_TEMP,icg) 
                    END DO ! icg

!---------------- Mixture of perfect gas and detonation products, Fred 9/10/15

                  CASE (GAS_MODEL_MIXT_JWL)
       IF (global%specUsed .EQV. .TRUE.) THEN
         IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
             scalarConvFlag = .FALSE.
         ELSE
             scalarConvFlag = .TRUE.
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
         END IF

                    DO icg = 1,pGrid%nCellsTot
                    r = pCv(CV_MIXT_DENS,icg)
                    ir = 1.0_RFREAL/r
                    p = pDv(DV_MIXT_PRES,icg)

                    pCv(CV_MIXT_XVEL,icg) = ir*pCv(CV_MIXT_XMOM,icg)
                    pCv(CV_MIXT_YVEL,icg) = ir*pCv(CV_MIXT_YMOM,icg)
                    pCv(CV_MIXT_ZVEL,icg) = ir*pCv(CV_MIXT_ZMOM,icg)

              iCvSpecAir = SPEC_GetSpeciesIndex(global,pSpecInput,'AIR')

                    mwAir = pSpecInput%specType(iCvSpecAir)%pMaterial%molW
                    cpAir = pSpecInput%specType(iCvSpecAir)%pMaterial%spht
                    gcAir = MixtPerf_R_M(mwAir)
                    gAir = MixtPerf_G_CpR(cpAir,gcAir)


                      IF ( global%piclUsed .EQV. .TRUE. ) THEN
                         r  = r/(1.0_RFREAL - pRegion%mixt%piclVF(icg))
                        END IF
                  iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
                  Yproducts =  pCvSpec(iCvSpecProducts,icg) !spec%cv is primitive here

      CALL RFLU_JWL_ComputeEnergyMixt(pRegion,icg,gAir,gcAir,p,r,Yproducts,a,b,T)

                     pCv(CV_MIXT_TEMP,icg) = T
                    END DO

        IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertcvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
        END IF
      END IF
                  CASE DEFAULT
                    CALL ErrorStop(global,ERR_REACHED_DEFAULT,345)
                END SELECT ! pRegion%mixtInput%gasModel              

! ----------- Default ----------------------------------------------------------

              CASE DEFAULT
                CALL ErrorStop(global,ERR_REACHED_DEFAULT,351)     
            END SELECT ! pRegion%mixtInput%fluidModel 
             
! ------------------------------------------------------------------------------          
!         Default            
! ------------------------------------------------------------------------------
  
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,359)   
        END SELECT ! cvStateFuture

! ==============================================================================
!     Error - invalid input
! ==============================================================================

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,367)
    END SELECT ! pRegion%mixt%cvStateFuture

! ******************************************************************************
!   End
! ******************************************************************************


    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ConvertCvCons2Prim




! ******************************************************************************
!
! Purpose: Convert primitive state vector to consverved variables.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to data of current region
!   cvStateFuture       Future state of conserved variables
!
! Output: None.
!
! Notes: 
!   1. Strictly speaking, cvStateFuture is not needed (there is only one
!      state for conserved variables), but kept for consistency with 
!      RFLU_ConvertCvCons2Prim.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ConvertCvPrim2Cons(pRegion,cvStateFuture)

    USE ModInterfaces, ONLY: MixtPerf_Cv_CpR, &
                             MixtGasLiq_Eo_CvmTVm2, &
                             MixtPerf_Eo_DGPUVW, &
                             MixtPerf_G_CpR, &
                             MixtPerf_R_M

    USE ModSpecies, ONLY: t_spec_input
    USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
    USE SPEC_RFLU_ModPBAUtils, ONLY: SPEC_RFLU_PBA_GetSolution

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: cvStateFuture
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,indCp,indMol
    REAL(RFREAL) :: cp,cpProducts,Cvm,cvg,cvl,cvv,d,dummyVal,e,g,gc, &
                    gcProducts,gProducts,mw,mwProducts,nTol,p,po,r,ro,Rg,Rv, &
                    rYg,rYl,rYv,t,u,ud,v,Vm2,vFrac,w,YExplosive,YProducts

    REAL(RFREAL) :: r_pure !Varaible to hold pure (non-vf weighted) density for conversions/EoS calls

    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
    LOGICAL :: scalarConvFlag
    INTEGER :: iCvSpecExplosive,iCvSpecProducts
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
    TYPE(t_spec_input), POINTER :: pSpecInput 
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvertCvPrim2Cons',"../modflu/RFLU_ModConvertCv.F90")


! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid
    pCv   => pRegion%mixt%cv
    pDv   => pRegion%mixt%dv
    pGv   => pRegion%mixt%gv

    pCvSpec    => pRegion%spec%cv
    pSpecInput => pRegion%specInput

    indCp  = pRegion%mixtInput%indCp
    indMol = pRegion%mixtInput%indMol

    !indVFracE = 0

    nTol = 1.0E-14_RFREAL

! ******************************************************************************
!   Actual conversion
! ******************************************************************************
  
    IF ( pRegion%mixt%cvState == CV_MIXT_STATE_DUVWP .OR. & 
         pRegion%mixt%cvState == CV_MIXT_STATE_DUVWT ) THEN 

! ==============================================================================
!     Convert from primitive to conservative form 
! ==============================================================================

      SELECT CASE ( cvStateFuture )
        CASE ( CV_MIXT_STATE_CONS )
          pRegion%mixt%cvState = CV_MIXT_STATE_CONS          

          SELECT CASE ( pRegion%mixtInput%gasModel )

! ------------------------------------------------------------------------------
!           Thermally and calorically perfect gas or pseudo-gas
! ------------------------------------------------------------------------------

            CASE ( GAS_MODEL_TCPERF, & 
                   GAS_MODEL_MIXT_TCPERF, & 
                   GAS_MODEL_MIXT_PSEUDO )
              IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
                DO icg = 1,pGrid%nCellsTot
                  r = pCv(CV_MIXT_DENS,icg)
                  u = pCv(CV_MIXT_XVEL,icg)
                  v = pCv(CV_MIXT_YVEL,icg)
                  w = pCv(CV_MIXT_ZVEL,icg)
                  p = pCv(CV_MIXT_PRES,icg)

                  pCv(CV_MIXT_XMOM,icg) = r*u
                  pCv(CV_MIXT_YMOM,icg) = r*v
                  pCv(CV_MIXT_ZMOM,icg) = r*w

                  g  = global%refGamma

                      pCv(CV_MIXT_ENER,icg) = r*MixtPerf_Eo_DGPUVW(r,g,p,u,v,w)

                      IF ( global%piclUsed .EQV. .TRUE. ) THEN
                       r_pure = r/(1.0_RFREAL-pRegion%mixt%piclVF(icg))
                       pCv(CV_MIXT_ENER,icg) = r*MixtPerf_Eo_DGPUVW(r_pure,g,p,u,v,w)
                      END IF


                END DO ! icg 
              ELSE
! ------------- No program burn ------------------------------------------------
                IF ( (global%pbaFlag .EQV. .FALSE.) .OR. &
                     (global%pbaBurnFlag .EQV. .FALSE.) ) THEN
                  DO icg = 1,pGrid%nCellsTot
                    r = pCv(CV_MIXT_DENS,icg)
                    u = pCv(CV_MIXT_XVEL,icg)
                    v = pCv(CV_MIXT_YVEL,icg)
                    w = pCv(CV_MIXT_ZVEL,icg)
                    p = pDv(DV_MIXT_PRES,icg)

                    pCv(CV_MIXT_XMOM,icg) = r*u
                    pCv(CV_MIXT_YMOM,icg) = r*v
                    pCv(CV_MIXT_ZMOM,icg) = r*w

                    cp = pGv(GV_MIXT_CP,indCp*icg)
                    mw = pGv(GV_MIXT_MOL,indMol*icg)                  
                    gc = MixtPerf_R_M(mw)
                    g  = MixtPerf_G_CpR(cp,gc)

                      pCv(CV_MIXT_ENER,icg) = r*MixtPerf_Eo_DGPUVW(r,g,p,u,v,w)

                      IF ( global%piclUsed .EQV. .TRUE. ) THEN
                       r_pure = r/(1.0_RFREAL-pRegion%mixt%piclVF(icg))
                       pCv(CV_MIXT_ENER,icg) =r*MixtPerf_Eo_DGPUVW(r_pure,g,p,u,v,w)
                      END IF

                  END DO ! icg 
! ------------- Program burn ---------------------------------------------------
                ELSE
                  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_TCPERF ) &
                    THEN
                    CALL ErrorStop(global,ERR_GASMODEL_INVALID,594, & 
                                   'Program burn needs mixture gas model.')              
                  END IF ! gasModel  

                  IF ( global%specUsed .EQV. .FALSE. ) THEN
                    CALL ErrorStop(global,ERR_REACHED_DEFAULT,599) ! Defensive coding
                  END IF ! global%specUsed

                  iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                          'EXPLOSIVE')
                  ud = pSpecInput%specType(iCvSpecExplosive)%pMaterial%detonVel
              
                  iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                         'PRODUCTS')
                  mwProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%molw
                  cpProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%spht
                  gcProducts = MixtPerf_R_M(mwProducts)
                  gProducts  = MixtPerf_G_CpR(cpProducts,gcProducts)

                  DO icg = 1,pGrid%nCellsTot
                    r = pCv(CV_MIXT_DENS,icg)
                    u = pCv(CV_MIXT_XVEL,icg)
                    v = pCv(CV_MIXT_YVEL,icg)
                    w = pCv(CV_MIXT_ZVEL,icg)
                    p = pDv(DV_MIXT_PRES,icg)

                    pCv(CV_MIXT_XMOM,icg) = r*u
                    pCv(CV_MIXT_YMOM,icg) = r*v
                    pCv(CV_MIXT_ZMOM,icg) = r*w

                    cp = pGv(GV_MIXT_CP,indCp*icg)
                    mw = pGv(GV_MIXT_MOL,indMol*icg)                  
                    gc = MixtPerf_R_M(mw)
                    g  = MixtPerf_G_CpR(cp,gc)

                    !vFrac = (1.0_RFREAL - pRegion%plag%vFracE(1,indVFracE*icg))


                      IF ( global%piclUsed .EQV. .TRUE. ) THEN
                       r_pure = r/(1.0_RFREAL-pRegion%mixt%piclVF(icg))
                      END IF

                    YExplosive = pCvSpec(iCvSpecExplosive,icg)/r
                    YProducts  = pCvSpec(iCvSpecProducts,icg)/r

! DEBUG: Manoj-PBA1D, Notes: To match Jianghui's implementation
!                         1. No special code for Program Burn here
!                Changing above back to what Manoj implemented
!                         2. Use same material properties for both species
! END DEBUG                    
                    IF ( ABS(YExplosive) < nTol ) THEN
                     pCv(CV_MIXT_ENER,icg) = r*MixtPerf_Eo_DGPUVW(r,g,p,u,v,w)

                      IF ( global%piclUsed .EQV. .TRUE. ) THEN
                       pCv(CV_MIXT_ENER,icg) =r*MixtPerf_Eo_DGPUVW(r_pure,g,p,u,v,w)
                      END IF
                    ELSEIF ( ABS(1.0_RFREAL-YExplosive) < nTol ) THEN
                      ro = pRegion%mixtInput%prepRealVal6
                      po = pRegion%mixtInput%prepRealVal7

                      CALL SPEC_RFLU_PBA_GetSolution(gProducts,gcProducts,ro,po,ud, &
                                                     YProducts,e,p,d,dummyVal,dummyVal)

                      pCv(CV_MIXT_ENER,icg) = d*e+0.5_RFREAL*d*(u*u+v*v+w*w)
                    ELSE
                      p                     = p/YProducts
                      pCv(CV_MIXT_ENER,icg) = r*MixtPerf_Eo_DGPUVW(r,g,p,u,v,w)
                      !Fred Note - This implemntation is most likely correct now
                      !if running program burn + plag...but should be double
                      !checked - .05/13/2022

                      !Brad Note - I concur with Fred but for program burn +
                      !picl - .04/19/2023  

                      IF ( global%piclUsed .EQV. .TRUE. ) THEN
                       pCv(CV_MIXT_ENER,icg) =r*MixtPerf_Eo_DGPUVW(r_pure,gProducts,p,u,v,w)
                      END IF
                    END IF ! YExplosive
                  END DO ! icg 
                END IF ! pbaFlag
              END IF ! global%solverType

! ------------------------------------------------------------------------------
!           Mixture of gas, liquid, and vapor
! ------------------------------------------------------------------------------

            CASE ( GAS_MODEL_MIXT_GASLIQ )  
              Rg  = MixtPerf_R_M(pRegion%specInput%specType(1)% &
                                 pMaterial%molw)
              cvg = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)% &
                                    pMaterial%spht,Rg)
              Rv  = MixtPerf_R_M(pRegion%specInput%specType(2)% &
                                 pMaterial%molw)
              cvv = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)% &
                                    pMaterial%spht,Rv)
              cvl = global%refCvLiq

              DO icg = 1,pGrid%nCellsTot
                r = pCv(CV_MIXT_DENS,icg)
                u = pCv(CV_MIXT_XVEL,icg)
                v = pCv(CV_MIXT_YVEL,icg)
                w = pCv(CV_MIXT_ZVEL,icg)
                t = pDv(DV_MIXT_TEMP,icg)

                pCv(CV_MIXT_XMOM,icg) = r*u
                pCv(CV_MIXT_YMOM,icg) = r*v
                pCv(CV_MIXT_ZMOM,icg) = r*w

                rYg = pCvSpec(1,icg)
                rYv = pCvSpec(2,icg)
                rYl = r - rYg - rYv 

                Cvm = (rYl*cvl + rYg*cvg + rYv*cvv)/r
                Vm2 = u*u + v*v + w*w

                pCv(CV_MIXT_ENER,icg) = r*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
              END DO ! icg

! ------------------------------------------------------------------------------
!           Mixture of thermally and calorically perfect gas and detonation
!           products
! ------------------------------------------------------------------------------

            CASE ( GAS_MODEL_MIXT_JWL )
              IF ( global%specUsed .EQV. .TRUE. ) THEN
                IF ( pRegion%spec%cvState == CV_MIXT_STATE_PRIM ) THEN 
                  scalarConvFlag = .FALSE.
                ELSE
                  scalarConvFlag = .TRUE.

                  CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv, & 
                                                     pRegion%spec%cvState)
                END IF ! pRegion%spec%cvState  

                iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                       'PRODUCTS')

                DO icg = 1,pGrid%nCellsTot
                  r = pCv(CV_MIXT_DENS,icg)
                  u = pCv(CV_MIXT_XVEL,icg)
                  v = pCv(CV_MIXT_YVEL,icg)
                  w = pCv(CV_MIXT_ZVEL,icg)
                  p = pDv(DV_MIXT_PRES,icg)

                  pCv(CV_MIXT_XMOM,icg) = r*u
                  pCv(CV_MIXT_YMOM,icg) = r*v
                  pCv(CV_MIXT_ZMOM,icg) = r*w
                  
                  ! Subbu - Correction
                  ! Since pCvSpec has been transformed to Prim state
                  ! YProducts should not be divided by density
                  ! YProducts = pCvSpec(iCvSpecProducts,icg)/r
                  YProducts = pCvSpec(iCvSpecProducts,icg)
                  ! Subbu - End Correction

                  e = YProducts*pDv(DV_MIXT_EJWL,icg) &
                    + (1.0_RFREAL-YProducts)*pDv(DV_MIXT_EPERF,icg) 
                  
                  pCv(CV_MIXT_ENER,icg) = r*(e + 0.5_RFREAL*(u*u+v*v+w*w))

                END DO ! icg 
                
                IF ( scalarConvFlag .EQV. .TRUE. ) THEN 
                  CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv, & 
                                                     pRegion%spec%cvState)
                END IF ! scalarConvFlag                
              END IF ! global%specUsed  

! ------------------------------------------------------------------------------
!           Other or invalid gas models             
! ------------------------------------------------------------------------------

            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,807)
          END SELECT ! pRegion%mixtInput
          
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,811)
      END SELECT ! cvStateFuture

! ==============================================================================
!   Error - invalid input
! ==============================================================================

    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,819)
    END IF ! pRegion%mixt%cvState                   
                         
! ******************************************************************************
!   End
! ******************************************************************************


    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ConvertCvPrim2Cons




! ******************************************************************************
!
! Purpose: Convert conserved state vector to primitive variables.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to data of current region
!   cvScal              State vector of conserved variables
!   cvScalStateCurrent  Current state of conserved variables
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ScalarConvertCvCons2Prim(pRegion,cvScal,cvScalStateCurrent)

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(INOUT) :: cvScalStateCurrent
    REAL(RFREAL), DIMENSION(:,:) :: cvScal
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,iScal,nScal
    REAL(RFREAL) :: ir
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ScalarConvertCvCons2Prim',"../modflu/RFLU_ModConvertCv.F90")

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid
    pCv   => pRegion%mixt%cv

    nScal = SIZE(cvScal,1)

! ******************************************************************************
!   Actual conversion
! ******************************************************************************

    SELECT CASE ( cvScalStateCurrent )

! ==============================================================================
!     Conserved, so convert to primitive variables     
! ==============================================================================

      CASE ( CV_MIXT_STATE_CONS )
        cvScalStateCurrent = CV_MIXT_STATE_PRIM         

        DO icg = 1,pGrid%nCellsTot
          ir = 1.0_RFREAL/pCv(CV_MIXT_DENS,icg)

          DO iScal = 1,nScal
            cvScal(iScal,icg) = ir*cvScal(iScal,icg)
          END DO ! iScal
        END DO ! icg

! ==============================================================================
!     Error - invalid input
! ==============================================================================

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,921)
    END SELECT ! cvScalStateCurrent

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ScalarConvertCvCons2Prim




! ******************************************************************************
!
! Purpose: Convert primitive state vector to conserved variables.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to data of current region
!   cvScal              State vector of conserved variables
!   cvScalStateCurrent  Current state of conserved variables
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ScalarConvertCvPrim2Cons(pRegion,cvScal,cvScalStateCurrent)

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(INOUT) :: cvScalStateCurrent
    REAL(RFREAL), DIMENSION(:,:) :: cvScal
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,iScal,nScal
    REAL(RFREAL) :: r
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ScalarConvertCvPrim2Cons',"../modflu/RFLU_ModConvertCv.F90")

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid
    pCv   => pRegion%mixt%cv

    nScal = SIZE(cvScal,1)

! ******************************************************************************
!   Actual conversion
! ******************************************************************************

    SELECT CASE ( cvScalStateCurrent )

! ==============================================================================
!     Primitive, so convert to conserved variables     
! ==============================================================================

      CASE ( CV_MIXT_STATE_PRIM )
        cvScalStateCurrent = CV_MIXT_STATE_CONS         

        DO icg = 1,pGrid%nCellsTot
          r = pCv(CV_MIXT_DENS,icg)

          DO iScal = 1,nScal
            cvScal(iScal,icg) = r*cvScal(iScal,icg)
          END DO ! iScal
        END DO ! icg

! ==============================================================================
!     Error - invalid input
! ==============================================================================

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,1019)
    END SELECT ! cvScalStateCurrent

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ScalarConvertCvPrim2Cons





! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModConvertCv


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModConvertCv.F90,v $
! Revision 1.3  2016/02/08 20:02:15  fred
! Fixing non-1 flag compiling error
!
! Revision 1.2  2016/02/04 19:57:48  fred
! Adding iterative JWL EOS capabilities for the cylindrical detonation case
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2010/03/14 23:45:20  mparmar
! Added support for SOLV_IMPLICIT_HM
!
! Revision 1.3  2008/12/06 08:43:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:40  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.8  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.7  2006/03/26 20:22:00  haselbac
! Added support for GL model, cosmetics
!
! Revision 1.6  2005/11/14 16:59:08  haselbac
! Added support for pseudo-gas model
!
! Revision 1.5  2005/11/10 02:26:54  haselbac
! Cosmetics only
!
! Revision 1.4  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.3  2005/07/07 22:45:01  haselbac
! Added profiling calls
!              
! Revision 1.2  2004/01/29 22:57:35  haselbac  
! Added routines for scalars                   
!
! Revision 1.1  2002/09/09 15:17:58  haselbac  
! Initial revision                             
!
! ******************************************************************************

