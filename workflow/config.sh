#!/usr/bin/env bash
# ABYSS workflow configuration (sourced by scripts in workflow/bin)
#
# Copy to workflow/config.local.sh to override locally.

# Which environment are we running on?
SYSTEM_NAME=${SYSTEM_NAME:-"ciera"}

# Scheduler: slurm | pbs | local
SCHEDULER=${SCHEDULER:-"slurm"}

# Job resources (scheduler-aware)
ACCOUNT=${ACCOUNT:-"b1094"}
PARTITION=${PARTITION:-"ciera-gpu"}   # Slurm partition
QUEUE=${QUEUE:-""}                    # PBS queue (optional)

NODES=${NODES:-1}
NTASKS=${NTASKS:-16}
CPUS_PER_TASK=${CPUS_PER_TASK:-1}
GPUS=${GPUS:-1}
WALLTIME=${WALLTIME:-"12:00:00"}       # Slurm HH:MM:SS, PBS HH:MM:SS

# Test selection
TEST_DIR=${TEST_DIR:-"test/test1"}
TEST_CONFIG=${TEST_CONFIG:-"config.txt"}

# Build options
# USE_CUDA: 1 to build CUDA variant; 0 CPU-only. "auto" tries to detect GPUs.
USE_CUDA=${USE_CUDA:-1}

# Preferred MPI compiler wrapper
CXX=${CXX:-"mpicxx"}

# Optional extra compile flags passed through to Makefile:
#  - Use to override MaxNumParticle / MaxNumNeighbor for small test runs.
EXTRA_CXXFLAGS=${EXTRA_CXXFLAGS:-""}
EXTRA_NVCCFLAGS=${EXTRA_NVCCFLAGS:-""}

# Preferred launcher inside allocations
# slurm: srun recommended; pbs: mpirun typically
MPI_LAUNCHER=${MPI_LAUNCHER:-"srun"}   # srun | mpirun | none
SRUN_MPI=${SRUN_MPI:-"pmix_v4"}

# Modules to load (optional). Keep empty if your cluster doesn't use modules.
# Example:
# MODULES=("cuda/cuda-12.1.0-openmpi-4.1.4")
MODULES=(${MODULES[@]:-})

# OpenMPI shared-memory backing file relocation
OMPI_BACKING_DIR=${OMPI_BACKING_DIR:-"/tmp"}

# Output locations
WORK_DIR=${WORK_DIR:-"workflow/runs"}  # run artifacts go here
