#!/bin/bash

# this script opens the gradient related subroutines in separate vim tabs

file1=libfloflu/CellGradientsMP.F90
file2=modfloflu/ModMixture.F90
file3=modflu/RFLU_ModAllocateMemory.F90
file4=modflu/RFLU_ModDeallocateMemory.F90
file5=modflu/RFLU_ModReadWriteFlow.F90

vim -p $file1 $file2 $file3 $file4 $file5
