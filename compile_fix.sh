#!/bin/bash
# Compilation script for ABYSS after snake_case refactoring fixes
# Run from the ABYSS root directory: ./compile_fix.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== ABYSS Compilation Script ==="
echo "Date: $(date)"

# Set up environment
export PATH="/hpc/software/spack_v20d1/spack/opt/spack/linux-rhel7-x86_64/gcc-12.3.0/openmpi-4.1.6rc2-4jmm4uy2wrgvpfqpoc6geocztczhw6yx/bin:$PATH"
export LD_LIBRARY_PATH="/hpc/software/spack_v20d1/spack/opt/spack/linux-rhel7-x86_64/gcc-12.3.0/openmpi-4.1.6rc2-4jmm4uy2wrgvpfqpoc6geocztczhw6yx/lib:${LD_LIBRARY_PATH:-}"
export HDF5_DIR="/hpc/software/2025/spack/opt/spack/linux-rhel8-x86_64/gcc-12.4.0/hdf5-1.14.5-yke7pax3ot3kyspoqxw6lc4qolvpk43t"

echo "mpicxx: $(which mpicxx)"
echo "HDF5_DIR: $HDF5_DIR"

cd src

echo ""
echo "=== Cleaning previous build ==="
make clean

echo ""
echo "=== Compiling (without CUDA) ==="
make -j8 2>&1 | tee ../build_output.log

echo ""
echo "=== Compilation complete ==="
echo "Build output saved to: $SCRIPT_DIR/build_output.log"
