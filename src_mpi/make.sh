#module add mpi/latest
# make clean
# make -j8

export LD_LIBRARY_PATH=/data/vinicius/sevn_n14/build/lib64/sevn:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/data/vinicius/sevn/build/lib64/sevn:$LD_LIBRARY_PATH

make clean USE_SEVN=1
make -j8 USE_SEVN=1

date
