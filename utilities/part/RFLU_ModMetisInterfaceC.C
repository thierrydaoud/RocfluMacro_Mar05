//*********************************************************************
//* Illinois Open Source License                                      *
//*                                                                   *
//* University of Illinois/NCSA                                       * 
//* Open Source License                                               *
//*                                                                   *
//* Copyright@2008, University of Illinois.  All rights reserved.     *
//*                                                                   *
//*  Developed by:                                                    *
//*                                                                   *
//*     Center for Simulation of Advanced Rockets                     *
//*                                                                   *
//*     University of Illinois                                        *
//*                                                                   *
//*     www.csar.uiuc.edu                                             *
//*                                                                   *
//* Permission is hereby granted, free of charge, to any person       *
//* obtaining a copy of this software and associated documentation    *
//* files (the "Software"), to deal with the Software without         *
//* restriction, including without limitation the rights to use,      *
//* copy, modify, merge, publish, distribute, sublicense, and/or      *
//* sell copies of the Software, and to permit persons to whom the    *
//* Software is furnished to do so, subject to the following          *
//* conditions:                                                       *
//*                                                                   *
//*                                                                   *
//* @ Redistributions of source code must retain the above copyright  * 
//*   notice, this list of conditions and the following disclaimers.  *
//*                                                                   * 
//* @ Redistributions in binary form must reproduce the above         *
//*   copyright notice, this list of conditions and the following     *
//*   disclaimers in the documentation and/or other materials         *
//*   provided with the distribution.                                 *
//*                                                                   *
//* @ Neither the names of the Center for Simulation of Advanced      *
//*   Rockets, the University of Illinois, nor the names of its       *
//*   contributors may be used to endorse or promote products derived * 
//*   from this Software without specific prior written permission.   *
//*                                                                   *
//* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
//* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
//* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
//* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
//* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
//* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
//* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
//* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
//*********************************************************************
//* Please acknowledge The University of Illinois Center for          *
//* Simulation of Advanced Rockets in works and publications          *
//* resulting from this software or its derivatives.                  *
//*********************************************************************
//*******************************************************************************
//
// Purpose: 
//
// Description: None.
//
// see RFLU_ModMetisInterface.F90 for an explanation of why this file exists. 
// More options can be added easily as needed. 
//
// need to include extern "C" becuase rocflu uses the c++ compiler
//
// Notes: 
//   1. The routine which builds the access lists MUST be called after the face 
//      lists have been built. 
//
//*******************************************************************************
//
// $Id: RFLU_ModBFaceGradAccessList.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
//
// Copyright: (c) 2003-2004 by the University of Illinois
//
//*******************************************************************************

#include <metis.h>


extern "C" void rflu_part_metis_noptions(int *output){
     *output = METIS_NOPTIONS;
     return;
}

extern "C" void rflu_part_metis_option_numbering(int *output){
     *output = METIS_OPTION_NUMBERING;
     return;
}

extern "C" void rflu_part_metis_option_contig(int *output){
     *output = METIS_OPTION_CONTIG;
     return;
}

extern "C" void rflu_part_metis_option_dbglvl(int *output){
     *output = METIS_OPTION_DBGLVL;
     return;
}

