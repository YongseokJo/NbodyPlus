# ABYSS Technology Stack

## Languages

| Language | Version | Usage |
|----------|---------|-------|
| C++ | C++11 | Core simulation code |
| CUDA | 12.x | GPU acceleration for force calculations |
| Python | 3.x | Analysis tools and post-processing |

## Build System

- **Makefile**: Primary build system (`src/Makefile`)
- **Compiler**: `mpicxx` (MPI C++ wrapper)
- **Flags**: `-O2 -march=native -ftree-vectorize -std=c++11 -Wall -Wextra -g`

## Core Dependencies

### Required

| Dependency | Purpose | Discovery |
|------------|---------|-----------|
| MPI (OpenMPI) | Distributed parallelism (root/worker model) | `mpicxx` wrapper |
| HDF5 | Hierarchical data output | `HDF5_DIR` or `HDF5_ROOT` env vars |

### Optional (Compile-time Toggles)

| Dependency | Toggle | Purpose |
|------------|--------|---------|
| CUDA | `USE_CUDA=1` | GPU-accelerated force calculations |
| SEVN | `USE_SEVN=1` | Stellar evolution library |
| NVTX | `NSIGHT` define | NVIDIA profiling markers |

## GPU Stack (when USE_CUDA=1)

- **CUDA Runtime**: `cuda_runtime.h`
- **cuBLAS**: `cublas_v2.h` (for matrix operations)
- **NVTX**: `nvToolsExt.h` (profiling)
- **Target GPU**: NVIDIA A100 (tuned parameters in `src/def.h`)

## Configuration Formats

| Format | Usage |
|--------|-------|
| TOML | New configuration format (`config.toml`) |
| Plain text | Legacy configuration format (`config.txt`) |

## Python Analysis Stack

Located in `tools/`:

- **numpy**: Numerical operations
- **h5py**: HDF5 file reading
- **matplotlib**: Visualization

## Environment Configuration

Primary configuration: `workflow/config.sh`

```bash
SCHEDULER="slurm"
USE_CUDA="1"
USE_SEVN="0"
```

## Key Compile-Time Defines

| Define | Purpose |
|--------|---------|
| `CUDA` | Enable CUDA code paths |
| `SEVN` | Enable stellar evolution |
| `FEWBODY` | Enable few-body dynamics (SDAR) |
| `AR_TTL` | SDAR time-transformed leapfrog |
| `AR_SLOWDOWN_TREE` | SDAR slowdown tree integration |
| `PERFORMANCETRACE` | Enable performance profiling |
| `MULTIMAP` | Use multimap for regular time scheduling |

## Executable

- **Output**: `src/abyss.exe`
- **Usage**: `abyss.exe -c config.toml` or `abyss.exe --config config.txt`
