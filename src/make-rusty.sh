#!/bin/bash

#module add modules
module add intel-oneapi-compilers
module add intel-oneapi-mpi
module add cuda

make clean
make -j4


#./nbodyplus.exe -f nbody.dat >stdout 2>stderr
#./nbodyplus.exe -f binary.dat >stdout 2>stderr
