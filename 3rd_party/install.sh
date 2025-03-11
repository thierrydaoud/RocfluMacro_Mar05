#!/bin/bash

SRC=$(pwd)/local

cd GKlib
set -e
## make config prefix= CONFIG_FLAGS="-D CMAKE_C_FLAGS=-D_POSIX_C_SOURCE=199309L -DCMAKE_INSTALL_PREFIX=$SRC"
make config prefix=$SRC
make
make install
cd ..

cd METIS
make config prefix=$SRC gklib_path=$SRC
make install
cd ..
