#!/bin/bash

module add modules
module add intel-oneapi-compilers
module add intel-oneapi-mpi

#./nbodyplus.exe
./nbodyplus.exe -f binary.dat -t 0.01 -dt 0.001 >stdout 2>stderr &
