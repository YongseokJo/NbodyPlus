#!/bin/bash

#module add modules
module add modules/2.2-20230808
#module add modules/2.1-20230222
module add intel-oneapi-compilers
module add intel-oneapi-mpi
#module add cuda
#module add gcc
module add openmpi



#make clean
#make -j8 -f ./Makefile_no_cuda
#mv ./abyss_mpi.exe abyss_mpi_no_cuda.exe


#module add modules/2.2-20230808
module add modules/2.3-20240529
module add intel-oneapi-compilers
module add intel-oneapi-mpi
module add cuda
module add openmpi #/cuda

make clean
make -j8
cp ./abyss.exe abyss-stable.exe


#./nbodyplus.exe -f nbody.dat >stdout 2>stderr
#./nbodyplus.exe -f binary.dat >stdout 2>stderr
