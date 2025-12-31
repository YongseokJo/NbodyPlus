#!/bin/bash

#SBATCH --job-name=abyss


#SBATCH -A b1094
###SBATCH --partition ciera-himem
#SBATCH --partition ciera-gpu

###SBATCH --partition ciera-std
###SBATCH --partition short
###SBATCH --partition ciera-std

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --gpus=1
###SBATCH --mem=1T

#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=g.kerex@gmail.com     # Where to send mail	
###SBATCH --time=120:00:00               # Time limit hrs:min:sec
#SBATCH --time=48:00:00               # Time limit hrs:min:sec


echo "JobID: $SLURM_JOB_ID"
echo "Node list: $SLURM_NODELIST"
echo "CPUs/task: $SLURM_CPUS_PER_TASK"
echo "Mem/CPU:   $SLURM_MEM_PER_CPU"
echo "Mem/node:  $SLURM_MEM_PER_NODE"

#free -h 
#grep -E 'MemTotal|MemFree|MemAvailable' /proc/meminfo


#lscpu

pwd; hostname; date
source ~/.bash_profile

module purge
# Use GCC + OpenMPI + CUDA (compatible toolchain)
module add mpi/openmpi-4.1.1-gcc.10.2.0
module add cuda/11.4.4-gcc-10.4.0

# Set CUDA paths explicitly
export CUDA_HOME=/hpc/software/spack_v20d1/spack/opt/spack/linux-rhel7-x86_64/gcc-10.4.0/cuda-11.4.4-e6jupats6qhpjoykra7r6bov6gfaol3w
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:$LD_LIBRARY_PATH"

# Suppress harmless InfiniBand transport warnings
export OMPI_MCA_mtl=^psm

cd $(pwd)

mpirun -np $SLURM_NTASKS ./abyss.exe -c config.txt >stdout 2>stderr



