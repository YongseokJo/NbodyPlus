#!/bin/bash
#
# ABYSS Build Script
# ------------------
# Compiles ABYSS with HDF5 output and TOML configuration support.
#
# Usage:
#   ./build.sh [options]
#
# Options:
#   --cuda       Enable CUDA GPU acceleration (requires CUDA toolkit)
#   --sevn       Enable SEVN stellar evolution library
#   --clean      Clean before building
#   --slurm      Submit build as a SLURM job (GPU nodes)
#   --test       Run test after successful build
#   --help       Show this help message
#
# SLURM Configuration:
#   Account:   b1094
#   Partition: ciera-gpu
#   Default:   1 GPU, 16 tasks, 12 hour walltime
#

set -e

# Get the directory where this script is located (works when run locally)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Parse arguments
USE_CUDA=""
USE_SEVN=""
DO_CLEAN=""
USE_SLURM=""
RUN_TEST=""

for arg in "$@"; do
    case $arg in
        --cuda)
            USE_CUDA="USE_CUDA=1"
            ;;
        --sevn)
            USE_SEVN="USE_SEVN=1"
            ;;
        --clean)
            DO_CLEAN="1"
            ;;
        --slurm)
            USE_SLURM="1"
            USE_CUDA="USE_CUDA=1"  # SLURM implies GPU/CUDA build
            ;;
        --test)
            RUN_TEST="1"
            ;;
        --help|-h)
            head -22 "$0" | tail -20
            exit 0
            ;;
    esac
done

# Ensure temp directory exists (needed for compilation)
export TMPDIR=${TMPDIR:-$HOME/tmp}
mkdir -p "$TMPDIR"

#=============================================================================
# SLURM Job Submission
#=============================================================================
if [ -n "$USE_SLURM" ]; then
    echo "========================================"
    echo "   ABYSS SLURM Build Submission"
    echo "========================================"
    echo ""
    echo "ABYSS Dir: $SCRIPT_DIR"
    echo "Account:   b1094"
    echo "Partition: ciera-gpu"
    echo ""

    # Determine test flag
    if [ -n "$RUN_TEST" ]; then
        TEST_FLAG="1"
    else
        TEST_FLAG="0"
    fi

    # Create SLURM job script with hardcoded ABYSS path
    SLURM_SCRIPT="$SCRIPT_DIR/slurm_build_gpu.sh"

    cat > "$SLURM_SCRIPT" << SLURM_EOF
#!/bin/bash
#SBATCH --account=b1094
#SBATCH --partition=ciera-gpu
#SBATCH --job-name=abyss_build
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=${SCRIPT_DIR}/logs/build_gpu_%j.log
#SBATCH --error=${SCRIPT_DIR}/logs/build_gpu_%j.log

set -e
set -o pipefail

# Hardcoded ABYSS directory (set at script generation time)
ABYSS_DIR="${SCRIPT_DIR}"

echo "========================================"
echo "   ABYSS GPU Build - SLURM Job"
echo "========================================"
echo ""
echo "Job ID:     \$SLURM_JOB_ID"
echo "Node:       \$SLURM_NODELIST"
echo "ABYSS Dir:  \$ABYSS_DIR"
echo "Started:    \$(date)"
echo ""

cd "\$ABYSS_DIR"

# Set temp directory
export TMPDIR=\${TMPDIR:-\$HOME/tmp}
mkdir -p "\$TMPDIR"

# Ensure a logs directory exists on the compute node
mkdir -p "\$ABYSS_DIR/logs"

#-----------------------------------------------------------------------------
# Load modules / Set environment
#-----------------------------------------------------------------------------
echo "[1/5] Setting up environment..."

# Set MPI path FIRST (before any checks) - use modern compiler
MPI_CANDIDATES=(
    "/hpc/software/spack_v20d1/spack/opt/spack/linux-rhel7-x86_64/gcc-12.3.0/openmpi-4.1.6rc2-4jmm4uy2wrgvpfqpoc6geocztczhw6yx"
    "/hpc/software/2025/spack/opt/spack/linux-rhel8-x86_64/gcc-12.4.0/openmpi-4.1.7-wn5igvodjfykt4bpgrf5n64iilky7udv"
)
for mpi_path in "\${MPI_CANDIDATES[@]}"; do
    if [ -f "\$mpi_path/bin/mpicxx" ]; then
        export PATH="\$mpi_path/bin:\$PATH"
        export LD_LIBRARY_PATH="\$mpi_path/lib:\$LD_LIBRARY_PATH"
        echo "  Using MPI: \$mpi_path"
        break
    fi
done

# Set CUDA path
CUDA_CANDIDATES=(
    "/hpc/software/spack_v20d1/spack/opt/spack/linux-rhel7-x86_64/gcc-12.3.0/cuda-12.4.1-bmnxx2e3tuup6fgzp3e7o4i4wumixq5j"
    "/usr/local/cuda"
    "/usr/local/cuda-12.0"
)
for cuda_path in "\${CUDA_CANDIDATES[@]}"; do
    if [ -f "\$cuda_path/bin/nvcc" ]; then
        export CUDA_HOME="\$cuda_path"
        export PATH="\$CUDA_HOME/bin:\$PATH"
        export LD_LIBRARY_PATH="\$CUDA_HOME/lib64:\$LD_LIBRARY_PATH"
        break
    fi
done

# Ensure nvcc uses a modern host compiler (needs C++11 support)
if [ -z "\${CUDAHOSTCXX:-}" ]; then
    if command -v g++ &>/dev/null; then
        export CUDAHOSTCXX="\$(command -v g++)"
    fi
fi

if [ -n "\${CUDAHOSTCXX:-}" ]; then
    echo "  CUDAHOSTCXX: \$CUDAHOSTCXX"
    # Basic sanity check: does the host compiler accept -std=c++11?
    if ! echo 'int main(){return 0;}' | "\$CUDAHOSTCXX" -std=c++11 -x c++ -c -o /tmp/abyss_cuda_hostcxx_test.o - 2>/dev/null; then
        echo "ERROR: CUDA host compiler (\$CUDAHOSTCXX) does not support -std=c++11."
        echo "Load a newer GCC/G++ module or set CUDAHOSTCXX to a C++11-capable compiler before building."
        exit 1
    fi
    rm -f /tmp/abyss_cuda_hostcxx_test.o
else
    echo "  CUDAHOSTCXX: not set"
fi

# Set HDF5 path
HDF5_CANDIDATES=(
    "/hpc/software/2025/spack/opt/spack/linux-rhel8-x86_64/gcc-12.4.0/hdf5-1.14.5-yke7pax3ot3kyspoqxw6lc4qolvpk43t"
    "/usr/local"
)
for hdf5_path in "\${HDF5_CANDIDATES[@]}"; do
    if [ -f "\$hdf5_path/include/H5Cpp.h" ]; then
        export HDF5_DIR="\$hdf5_path"
        export LD_LIBRARY_PATH="\$HDF5_DIR/lib:\$LD_LIBRARY_PATH"
        break
    fi
done

echo "  CUDA_HOME: \${CUDA_HOME:-not set}"
echo "  HDF5_DIR:  \${HDF5_DIR:-not set}"
echo "  TMPDIR:    \$TMPDIR"

#-----------------------------------------------------------------------------
# Verify compilers
#-----------------------------------------------------------------------------
echo ""
echo "[2/5] Checking compilers..."

if ! command -v mpicxx &> /dev/null; then
    echo "ERROR: mpicxx not found!"
    exit 1
fi
echo "  mpicxx: \$(which mpicxx)"

if ! command -v nvcc &> /dev/null; then
    echo "ERROR: nvcc not found!"
    exit 1
fi
echo "  nvcc:   \$(which nvcc)"
nvcc --version | head -4

# Show GPU info
echo ""
echo "[3/5] GPU Information..."
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader 2>/dev/null || echo "  nvidia-smi not available"

#-----------------------------------------------------------------------------
# Build
#-----------------------------------------------------------------------------
echo ""
echo "[4/5] Building ABYSS with CUDA..."
cd "\$ABYSS_DIR/src"

# Clean first
echo "  Cleaning previous build..."
make clean

# Build with CUDA
echo "  Compiling with USE_CUDA=1..."
BUILD_LOG="\$ABYSS_DIR/logs/build_output_\$(date +%Y%m%d_%H%M%S).log"
make USE_CUDA=1 2>&1 | tee "\$BUILD_LOG"
test \${PIPESTATUS[0]} -eq 0
echo "  Build log: \$BUILD_LOG"

# Check result
if [ -f abyss.exe ]; then
    echo ""
    echo "  Build successful!"
    ls -lh abyss.exe
else
    echo ""
    echo "  Build FAILED! Check build_output.log for details."
    exit 1
fi

#-----------------------------------------------------------------------------
# Test
#-----------------------------------------------------------------------------
RUN_TEST_FLAG="${TEST_FLAG}"
if [ "\$RUN_TEST_FLAG" = "1" ]; then
    echo ""
    echo "[5/5] Running test..."

    TEST_DIR="\$ABYSS_DIR/test/test1"
    if [ -d "\$TEST_DIR" ]; then
        cd "\$TEST_DIR"

        # Copy executable
        cp "\$ABYSS_DIR/src/abyss.exe" ./abyss_gpu.exe

        # Create test config with new TOML format
        cat > config_gpu_test.toml << 'TESTEOF'
# ABYSS Test Configuration (new TOML format)
Filename = "nbody.dat"
StopTime = 1.0e6
OutputDirectory = "output_gpu"

[numerics]
eta = 0.01
FixNumNeighbor = 100
InitialRadius = 0.2
RSearch = 2.5e-4
TSearch = 1e-6

[output]
dtOutput = 1.0e5
Compression = true
CompressionLevel = 6

[restart]
Enabled = false
TESTEOF

        # Remove old output
        rm -rf output_gpu

        # Run with 4 MPI processes
        echo "  Running: mpirun -np 4 ./abyss_gpu.exe -c config_gpu_test.toml"
        mpirun -np 4 ./abyss_gpu.exe -c config_gpu_test.toml 2>&1 | tee "\$ABYSS_DIR/logs/gpu_test_\$(date +%Y%m%d_%H%M%S).log"

        echo ""
        if [ -d output_gpu ]; then
            echo "  Test completed! Output directory created."
            ls -la output_gpu/
        else
            echo "  Test may have failed - no output directory."
        fi
    else
        echo "  Test directory not found: \$TEST_DIR"
    fi
else
    echo ""
    echo "[5/5] Skipping test (use --test to enable)"
fi

echo ""
echo "========================================"
echo "   Build Complete"
echo "========================================"
echo "Finished: \$(date)"
SLURM_EOF

    chmod +x "$SLURM_SCRIPT"

    echo "SLURM script created: $SLURM_SCRIPT"
    echo ""
    echo "Submitting job..."

    # Submit the job
    JOB_OUTPUT=$(sbatch "$SLURM_SCRIPT" 2>&1)
    echo "$JOB_OUTPUT"

    # Extract job ID
    JOB_ID=$(echo "$JOB_OUTPUT" | grep -oP 'Submitted batch job \K\d+' || echo "")

    if [ -n "$JOB_ID" ]; then
        echo ""
        echo "========================================"
        echo "   Job Submitted Successfully"
        echo "========================================"
        echo ""
        echo "Job ID:     $JOB_ID"
        echo "Log file:   ${SCRIPT_DIR}/build_gpu_${JOB_ID}.log"
        echo ""
        echo "Monitor with:"
        echo "  squeue -j $JOB_ID"
        echo "  tail -f ${SCRIPT_DIR}/build_gpu_${JOB_ID}.log"
        echo ""
        echo "Cancel with:"
        echo "  scancel $JOB_ID"
    fi

        # If submission succeeded, remove the temporary SLURM script to avoid clutter.
        if [ -n "$JOB_ID" ]; then
            rm -f "$SLURM_SCRIPT" && echo "Removed temporary SLURM script: $SLURM_SCRIPT"
        else
            echo "SLURM submission failed; keeping $SLURM_SCRIPT for debugging."
        fi

    exit 0
fi

#=============================================================================
# Local Build (non-SLURM)
#=============================================================================
echo "========================================"
echo "        ABYSS Build Script"
echo "========================================"
echo ""

# Try to load modules if available
if command -v module &> /dev/null; then
    echo "[1/4] Loading modules..."
    module load hdf5/1.14.1-2-openmpi-gcc-12.3.0 2>/dev/null || {
        module load hdf5 2>/dev/null || true
    }
    if [ -n "$USE_CUDA" ]; then
        module load cuda 2>/dev/null || true
    fi
fi

# Set paths for HDF5 and MPI
if [ -z "$HDF5_DIR" ]; then
    HDF5_CANDIDATES=(
        "/hpc/software/2025/spack/opt/spack/linux-rhel8-x86_64/gcc-12.4.0/hdf5-1.14.5-yke7pax3ot3kyspoqxw6lc4qolvpk43t"
        "/usr/local"
        "/usr"
    )
    for hdf5_path in "${HDF5_CANDIDATES[@]}"; do
        if [ -f "$hdf5_path/include/H5Cpp.h" ]; then
            export HDF5_DIR="$hdf5_path"
            break
        fi
    done
fi

# Set MPI path if mpicxx is not in PATH
if ! command -v mpicxx &> /dev/null; then
    MPI_CANDIDATES=(
        "/hpc/software/2025/spack/opt/spack/linux-rhel8-x86_64/gcc-12.4.0/openmpi-4.1.7-wn5igvodjfykt4bpgrf5n64iilky7udv"
    )
    for mpi_path in "${MPI_CANDIDATES[@]}"; do
        if [ -f "$mpi_path/bin/mpicxx" ]; then
            export PATH="$mpi_path/bin:$PATH"
            export LD_LIBRARY_PATH="$mpi_path/lib:$LD_LIBRARY_PATH"
            break
        fi
    done
fi

# Add HDF5 library path
if [ -n "$HDF5_DIR" ]; then
    export LD_LIBRARY_PATH="$HDF5_DIR/lib:$LD_LIBRARY_PATH"
fi

echo "[2/4] Checking environment..."
echo "  HDF5_DIR: ${HDF5_DIR:-not set}"
echo "  TMPDIR:   $TMPDIR"

# Ensure a local logs directory exists for non-SLURM builds
mkdir -p "$SCRIPT_DIR/logs"

# Verify compilers are available
if ! command -v mpicxx &> /dev/null; then
    echo ""
    echo "ERROR: mpicxx not found!"
    echo "Please load an MPI module or set PATH to include MPI binaries."
    exit 1
fi
echo "  mpicxx:   $(which mpicxx)"

if [ -n "$USE_CUDA" ]; then
    if ! command -v nvcc &> /dev/null; then
        echo ""
        echo "ERROR: nvcc not found!"
        echo "Please load CUDA module for GPU builds."
        exit 1
    fi
    echo "  nvcc:     $(which nvcc)"

    # Ensure nvcc uses a C++11-capable host compiler.
    if [ -z "${CUDAHOSTCXX:-}" ] && command -v g++ &>/dev/null; then
        export CUDAHOSTCXX="$(command -v g++)"
    fi
    if [ -n "${CUDAHOSTCXX:-}" ]; then
        echo "  CUDAHOSTCXX: $CUDAHOSTCXX"
        if ! echo 'int main(){return 0;}' | "$CUDAHOSTCXX" -std=c++11 -x c++ -c -o /tmp/abyss_cuda_hostcxx_test.o - 2>/dev/null; then
            echo ""
            echo "ERROR: CUDA host compiler ($CUDAHOSTCXX) does not support -std=c++11."
            echo "Load a newer GCC/G++ module or set CUDAHOSTCXX to a C++11-capable compiler before building."
            exit 1
        fi
        rm -f /tmp/abyss_cuda_hostcxx_test.o
    else
        echo "  CUDAHOSTCXX: not set"
    fi
fi

# Build
echo ""
echo "[3/4] Building ABYSS..."
cd "$SCRIPT_DIR/src"

if [ -n "$DO_CLEAN" ]; then
    echo "  Cleaning..."
    make clean
fi

echo "  Compiling..."
BUILD_LOG_LOCAL="$SCRIPT_DIR/logs/build_output_local_$(date +%Y%m%d_%H%M%S).log"
make $USE_CUDA $USE_SEVN 2>&1 | tee "$BUILD_LOG_LOCAL"
test ${PIPESTATUS[0]} -eq 0
echo "  Build log: $BUILD_LOG_LOCAL"

# Check result
echo ""
echo "[4/4] Checking result..."
if [ -f abyss.exe ]; then
    echo ""
    echo "========================================"
    echo "    Build successful!"
    echo "========================================"
    echo ""
    echo "Executable: $SCRIPT_DIR/src/abyss.exe"
    ls -lh abyss.exe

    # Run test if requested
    if [ -n "$RUN_TEST" ]; then
        echo ""
        echo "========================================"
        echo "    Running Test"
        echo "========================================"
        TEST_DIR="$SCRIPT_DIR/test/test1"
        if [ -d "$TEST_DIR" ]; then
            cd "$TEST_DIR"
            cp "$SCRIPT_DIR/src/abyss.exe" ./abyss_test.exe

            # Create test config with new TOML format
            cat > config_quick.toml << 'TESTEOF'
# ABYSS Test Configuration (new TOML format)
Filename = "nbody.dat"
StopTime = 1.0e6
OutputDirectory = "output_test"

[numerics]
eta = 0.01
FixNumNeighbor = 100
InitialRadius = 0.2
RSearch = 2.5e-4
TSearch = 1e-6

[output]
dtOutput = 1.0e5
Compression = true
CompressionLevel = 6

[restart]
Enabled = false
TESTEOF
            rm -rf output_test
            echo ""
            echo "Running: mpirun -np 2 ./abyss_test.exe -c config_quick.toml"
            mpirun -np 2 ./abyss_test.exe -c config_quick.toml
            echo ""
            if [ -d output_test ]; then
                echo "Test completed! Output:"
                ls -la output_test/
            fi
        else
            echo "Test directory not found: $TEST_DIR"
        fi
    else
        echo ""
        echo "To run (example with 4 MPI processes):"
        echo "  cd test/test1"
        echo "  mpirun -np 4 ../../src/abyss.exe -c config.toml"
    fi
else
    echo ""
    echo "========================================"
    echo "    Build FAILED"
    echo "========================================"
    echo ""
    echo "Check the error messages above for details."
    exit 1
fi
