module load intel21
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export USE_CUDA=1
#export USE_SEVN=1
make clean
make -j8
date
