#!/bin/bash

module add modules/2.1.1-20230405
module add intel-oneapi-compilers
module add intel-oneapi-mpi

make clean
make
./nbodyplus -f binary.dat >stdout 2>stderr
