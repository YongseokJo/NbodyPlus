export PATH=$PATH:/home/vinicius/install/mpich-3.3.2/mpich3_15/bin
#export PATH=$PATH:/home/vinicius/install/mpich-3.3.2/mpich3/bin
export USE_CUDA=1
#export USE_SEVN=1
make clean
make -j8
date
