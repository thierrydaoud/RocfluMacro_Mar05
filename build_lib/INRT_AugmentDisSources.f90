










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
! Purpose: converts the transfers of primary quantities along Edges for each
!          particle into augmentations of RHS terms for all quantities, for
!          an interaction involving Lagrangian particles
!
! Description: none.
!
! Input: iInrt = index of interaction
!
! Output: augments region%levels(iLev)%...%rhs structures
!
! Notes:
!
!   The RHS structures use opposite sign as the input source structure
!
!   For efficiency, this routine requires Nodes to be stored in this order:
!     Mixture, Lagrangian particle, Eulerian particle, Internal
!
!   The energy corresponding to a gas mass is actually taken to be the
!   enthalpy because energy added in reactions are measured as enthalpies
!
!******************************************************************************
!
! $Id: INRT_AugmentDisSources.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_AugmentDisSources( region,iInrt )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract,   ONLY : t_inrt_input,t_inrt_interact
  USE ModMixture,    ONLY : t_mixt
  USE ModPartLag,    ONLY : t_plag
  USE ModSpecies,    ONLY : t_spec
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region
  INTEGER,        INTENT(IN)            :: iInrt

! ... loop variables
  INTEGER :: iPcls,iPlag,iPeul,iEdge

! ... local variables
  INTEGER, PARAMETER :: MAX_NODES = 101

  CHARACTER(CHRLEN)  :: RCSIdentString

  LOGICAL :: computeAux

  INTEGER :: nPcls,nPlag,nPeul,nIntl,nInputEdges,nNodes,nEdges
  INTEGER :: indMixt,indPlag0,indPeul0,indIntl
  INTEGER :: indPlagVapor,indPlagn,indPeuln
  INTEGER :: indCp,gasModel,ic,iNod,off1beg,ic1beg
  INTEGER :: tEdge,iNode(2),token(2)
  INTEGER, POINTER :: pCvPlagMass(:), aiv(:,:)

  REAL(RFREAL) :: factorLimitForce,factorImpulseX,factorImpulseY,factorImpulseZ
  REAL(RFREAL) :: temp1,temp2,spht,massdot,hcapdot,enerdot
  REAL(RFREAL) :: kinedot1,kinedot2,thrmdot1,thrmdot2,enerdot1,enerdot2
  REAL(RFREAL) :: intlMass,intlEner,intlTemp,intlHcap,contFac
  REAL(RFREAL) :: sphtPlag(MAX_NODES),sphtPeul(MAX_NODES)
  REAL(RFREAL), DIMENSION(3) :: velo1,velo2,intlVelo,intlMome
  REAL(RFREAL), DIMENSION(3) :: momedot,momedot1,momedot2
  REAL(RFREAL)               :: src(MAX_NODES,5)
  REAL(RFREAL), DIMENSION(:),   POINTER :: p1begMixtTemp,pPlagTemp
  REAL(RFREAL), DIMENSION(:,:), POINTER :: p1begMixtVelo,pPlagVelo
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pMixtRhs,pPlagRhs,pPeulRhs
  REAL(RFREAL), DIMENSION(:,:), POINTER :: gv,arv,primary

  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_mixt),          POINTER :: pMixt
  TYPE(t_plag),          POINTER :: pPlag
  TYPE(t_spec),          POINTER :: pPeul
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_AugmentDisSources.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_AugmentDisSources',"../rocinteract/INRT_AugmentDisSources.F90" )

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_AugmentDisSources

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_AugmentDisSources.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
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
! Revision 1.3  2006/02/15 20:18:30  wasistho
! put peul within ifdef
!
! Revision 1.2  2005/10/31 21:09:37  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/12/01 21:56:10  fnajjar
! Initial revision after changing case
!
! Revision 1.15  2004/07/28 15:42:12  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.14  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.13  2004/03/25 21:14:53  jferry
! fixed pointer offset bug
!
! Revision 1.12  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.11  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.10  2004/01/31 03:59:22  haselbac
! Initial integration for Rocflu and Rocpart
!
! Revision 1.9  2003/09/19 20:35:26  jferry
! Implemented oxidizer species for burning interaction
!
! Revision 1.8  2003/05/08 17:17:14  jferry
! changed energy associated with mass to enthalpy
!
! Revision 1.7  2003/05/07 15:13:10  jferry
! Rearranged for efficiency
!
! Revision 1.6  2003/04/09 15:02:39  jferry
! removed erroneous volume normalization for continuum rhs
!
! Revision 1.5  2003/04/03 21:10:17  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.4  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.3  2003/03/24 23:30:52  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.2  2003/03/11 16:05:54  jferry
! Created data type for material properties
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************

