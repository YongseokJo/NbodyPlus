#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

JOBS="${BUILD_JOBS:-8}"

if command -v module >/dev/null 2>&1; then
  module add modules/2.3-20240529
  module add intel-oneapi-compilers
  module add intel-oneapi-mpi
  module add cuda
  module add openmpi
else
  echo "module command not found; ensure CUDA and MPI toolchains are on PATH." >&2
fi

if [[ -z "${SLURM_JOB_ID:-}" && -z "${CUDA_VISIBLE_DEVICES:-}" ]]; then
  echo "No GPU allocation detected." >&2
  echo "If your cluster requires a GPU node, start one with:" >&2
  echo "  srun --pty -p <gpu-partition> -G 1 -c 8 --mem=0 --time=02:00:00 bash" >&2
fi

export USE_CUDA=1
make clean
make -j"${JOBS}"
