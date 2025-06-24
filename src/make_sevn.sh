export LD_LIBRARY_PATH=/usr/local/cuda-11.8/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
module add icc/latest
export PATH=/home/vinicius/install/openmpi-4.0.5/bin:$PATH
export LD_LIBRARY_PATH=/home/vinicius/install/openmpi-4.0.5/lib:$LD_LIBRARY_PATH

#module add mpi/latest
make clean
make sevn -j8
#cp nbodyplus_binary.exe /data1/wispedia/nbody/purenbody
#cp abyss_mpi.exe ../../purenbody/
date
