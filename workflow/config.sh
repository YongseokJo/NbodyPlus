#!/usr/bin/env bash

# ABYSS workflow config (sourceable)
# Override locally by creating workflow/config.local.sh (gitignored).

# Scheduler: slurm | pbs | local
SCHEDULER="slurm"

# Slurm defaults (Ciera GPU)
ACCOUNT="b1094"
PARTITION="ciera-gpu"
WALLTIME="12:00:00"
NODES="1"
NTASKS="16"
CPUS_PER_TASK="1"
GPUS="1"
SRUN_MPI="pmix_v4"

# Build toggles
USE_CUDA="1"
USE_SEVN="0"

# What to run
TEST_DIR="test/test1"
RUN_CONFIG="config.toml"

# Analysis
# Python interpreter used by workflow/bin/analyze.sh when running tools/*.py
# Override this in workflow/config.local.sh, e.g.:
PYTHON=/gpfs/home/vjl4366/pyenv/venv/bin/python

# Which tool scripts to run (relative to repo root). Empty means "summary only".
# Examples:
#   ANALYZE_TOOLS=(tools/analyze_profiling.py tools/analyze_energy.py)
ANALYZE_TOOLS=(
  tools/analyze_profiling.py
  tools/analyze_energy.py
)

# Executable path produced by build (relative to repo root)
EXECUTABLE="src/abyss.exe"

# Environment discovery (edit for your cluster)
MPI_CANDIDATES=(
  "/hpc/software/spack_v20d1/spack/opt/spack/linux-rhel7-x86_64/gcc-12.3.0/openmpi-4.1.6rc2-4jmm4uy2wrgvpfqpoc6geocztczhw6yx"
  "/hpc/software/2025/spack/opt/spack/linux-rhel8-x86_64/gcc-12.4.0/openmpi-4.1.7-wn5igvodjfykt4bpgrf5n64iilky7udv"
)

CUDA_CANDIDATES=(
  "/hpc/software/spack_v20d1/spack/opt/spack/linux-rhel7-x86_64/gcc-12.3.0/cuda-12.4.1-bmnxx2e3tuup6fgzp3e7o4i4wumixq5j"
  "/usr/local/cuda"
  "/usr/local/cuda-12.0"
)

HDF5_CANDIDATES=(
  "/hpc/software/2025/spack/opt/spack/linux-rhel8-x86_64/gcc-12.4.0/hdf5-1.14.5-yke7pax3ot3kyspoqxw6lc4qolvpk43t"
  "/usr/local"
  "/usr"
)
