export PATH=$PATH:/appl/intel/oneapi/compiler/2021.4.0/linux/bin/intel64
which icc
#module load icc/latest
make clean
make -j4
cp nbodyplus_binary.exe /data1/wispedia/nbody/purenbody
date
