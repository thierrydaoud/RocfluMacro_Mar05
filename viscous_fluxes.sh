#!/bin/bash

# this script opens the viscous flux related files in separate vim tabs

file1=libfloflu/RungeKuttaMP.F90
file2=libfloflu/ConvectiveFluxes.F90
file3=modflu/RFLU_ModViscousFlux.F90
file4=libflu/RFLU_DerivedInputValues.F90
file5=modflu/RFLU_ModAUSMPlusUpFlux.F90
file6=libfloflu/SourceTerms.F90
file7=rocpicl/PICL_TEMP_Runge.F90

vim -p $file1 $file2 $file3 $file4 $file5 $file6 $file7
