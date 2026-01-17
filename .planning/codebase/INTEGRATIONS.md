# ABYSS External Integrations

## MPI (Message Passing Interface)

### Purpose
Distributed parallelism with root/worker model.

### Usage
```cpp
#include <mpi.h>

// Initialization in main.cpp
initializeMPI(argc, argv);

// Shared memory windows for particle data
MPI_Win win;
MPI_Win_allocate_shared(size, sizeof(Particle), MPI_INFO_NULL, shared_comm, &particles, &win);
```

### Key APIs Used
- `MPI_Init`, `MPI_Finalize`
- `MPI_Comm_rank`, `MPI_Comm_size`
- `MPI_Win_allocate_shared`, `MPI_Win_free`
- `MPI_Type_create_struct` (custom datatypes)
- `MPI_Comm_split_type` (shared memory communicator)

## HDF5 (Hierarchical Data Format)

### Purpose
Simulation output storage.

### Usage
Output written via `src/read_write.cpp`:
```cpp
#include <hdf5.h>

// Write timestep data to HDF5
writeParticle(global_time, output_num);
```

### Output Structure
```
output.h5
├── Step_0/
│   ├── Time_Myr (attribute)
│   ├── Mass_Msun, X_pc, Y_pc, Z_pc (datasets)
│   ├── Vx_km_s, Vy_km_s, Vz_km_s (datasets)
│   └── E_binary, E_merger, E_PN (attributes)
```

### Configuration
```toml
[output]
Compression = true
CompressionLevel = 6
```

## CUDA (GPU Acceleration)

### Purpose
GPU-accelerated N-body force calculations.

### Toggle
```bash
USE_CUDA=1  # in workflow/config.sh or make USE_CUDA=1
```

### Key Files
| File | Purpose |
|------|---------|
| `src/cuda/cuda_acceleration.cu` | Multi-GPU force dispatch |
| `src/cuda/cuda_kernels.cu` | CUDA kernels |
| `src/cuda/cuda_routines.cu` | GPU initialization |

### GPU Parameters (tuned for A100)
```cpp
// src/def.h
#define NBODY_MAX       100000000
#define BATCH_SIZE      64
#define GRID_DIM_Y      32
#define NNB_PER_BLOCK   128
```

### APIs Used
- CUDA Runtime API (`cuda_runtime.h`)
- cuBLAS (`cublas_v2.h`)
- NVTX for profiling (`nvToolsExt.h`)

## SDAR (Slow-Down Algorithmic Regularization)

### Purpose
Regularized integration for close encounters (binaries, multiples).

### Location
Vendored in `SDAR/` directory.

### Integration
```cpp
// src/FewBody/group.h
#include "AR/symplectic_integrator.h"

struct Group {
    AR::TimeTransformedSymplecticIntegrator<...> sym_int;
    AR::TimeTransformedSymplecticManager<...> manager;
};
```

### Key Features
- Time-transformed symplectic integration
- Binary tree structure for hierarchical systems
- Slowdown factors for efficient integration
- Interrupt handling for mergers/collisions

### Compile Flags
```makefile
-DAR_TTL -DAR_SLOWDOWN_TREE -DAR_SLOWDOWN_TIMESCALE -DFEWBODY
```

## SEVN (Stellar EVolution N-body)

### Purpose
Stellar evolution library for tracking star types, radii, and masses.

### Toggle
```bash
USE_SEVN=1
SEVN_DIR=/path/to/sevn
```

### Integration
```cpp
#ifdef SEVN
#include "sevn.h"
#include "star.h"
#include "binstar.h"

// In Particle struct
StarSEVN* stellar_evolution;
Binstar* binary_evolution;
#endif
```

### Features
- Star phase tracking (main sequence, giants, remnants)
- Mass loss and radius evolution
- Binary star evolution
- Supernova and compact remnant formation

## Enzo Integration

### Purpose
Coupling with Enzo hydrodynamics code (Enzo-Abyss).

### Status
- Separate repository: https://github.com/YongseokJo/enzo-nbody
- ABYSS provides N-body dynamics within Enzo simulations
- `enzo_time_step` variable for time coordination

### Related Code
```cpp
extern double enzo_time_step;  // Time step from Enzo
```

## Python Analysis Tools

### Dependencies
```python
import numpy as np
import h5py
import matplotlib.pyplot as plt
```

### Tools
| Script | Purpose |
|--------|---------|
| `tools/analyze_energy.py` | Energy conservation analysis |
| `tools/analyze_profiling.py` | Performance profiling |
| `tools/read_hdf5_output.py` | HDF5 data extraction |
| `tools/summarize_run.py` | Run summary generation |

## Configuration Formats

### TOML (Modern)
```toml
Filename = "nbody.dat"
StopTime = 1.0e7

[numerics]
eta = 0.01

[output]
dtOutput = 1.0e6
```

### Plain Text (Legacy)
```
Filename = nbody.dat
eta = 0.01
StopTime = 1.0e7
```

## Cluster Integration

### SLURM
Primary scheduler supported via `workflow/config.sh`:
```bash
SCHEDULER="slurm"
ACCOUNT="b1094"
PARTITION="ciera-gpu"
```

### Module System
```bash
module add git
module load cuda openmpi hdf5
```
