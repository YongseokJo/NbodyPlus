# ABYSS Architecture

## System Overview

ABYSS is a direct N-body simulation code for stellar dynamics with:
- 4th-order Hermite integrator with block timesteps
- GPU-accelerated force calculations
- MPI-based parallelism (root/worker model)
- Few-body subsystem integration via SDAR
- Optional stellar evolution via SEVN

## Architectural Pattern

**Root/Worker MPI Model**:
- **Root Process (rank 0)**: Orchestrates simulation, manages global state, handles I/O
- **Worker Processes (rank 1..N)**: Execute force calculations and few-body integrations

## Core Components

### 1. Main Simulation Loop (`src/root_routines.cpp`)

```
┌─────────────────────────────────────────────────────────┐
│                    Root Process                          │
├─────────────────────────────────────────────────────────┤
│  ┌─────────────────────────────────────────────────┐   │
│  │ 1. Update next regular time                      │   │
│  │ 2. Irregular routines (short timestep particles)│   │
│  │ 3. Regular routines (long timestep particles)   │   │
│  │ 4. Stellar evolution (SEVN)                     │   │
│  │ 5. Output (HDF5)                                │   │
│  └─────────────────────────────────────────────────┘   │
│                         │                               │
│                         ▼                               │
│  ┌─────────────────────────────────────────────────┐   │
│  │ Dispatch tasks to workers via MPI shared memory  │   │
│  └─────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│                   Worker Processes                       │
├─────────────────────────────────────────────────────────┤
│  - Compute irregular/regular accelerations              │
│  - GPU force calculations (if CUDA enabled)             │
│  - Few-body integrations (SDAR)                         │
└─────────────────────────────────────────────────────────┘
```

### 2. Particle System (`src/particle.h`)

The `Particle` struct contains:
- **Core properties**: `pid`, `position[3]`, `velocity[3]`, `mass`
- **Acceleration arrays**: `acc_total`, `acc_regular`, `acc_irregular` (4th-order Hermite)
- **Neighbor info**: `num_neighbors`, `neighbors_offset`
- **Time stepping**: Separate irregular and regular timesteps with block time management
- **Binary state**: `binary_state` for tracking binary interactions
- **SEVN integration**: Optional stellar evolution pointers

### 3. Integration Scheme

**Block Timestep Hermite**:
- Particles have individual timesteps quantized to powers of 2
- **Regular timesteps**: Far-field interactions (GPU-accelerated)
- **Irregular timesteps**: Near-field interactions (neighbors)

### 4. Few-Body Dynamics (`src/FewBody/`)

SDAR (Slow-Down Algorithmic Regularization) for close encounters:
- `Group` struct manages binary/multiple systems
- `AR::TimeTransformedSymplecticIntegrator` for regularized integration
- Post-Newtonian corrections (PN1.0, PN2.0, PN2.5 for GW energy loss)
- Merger handling (GW-driven, TDE, stellar collisions)

### 5. GPU Acceleration (`src/cuda/`)

- `cuda_acceleration.cu`: Multi-GPU force calculation with kernel launches
- `cuda_kernels.cu`: CUDA kernels for N-body force computation
- Optimized for NVIDIA A100 (batch size, grid dimensions in `def.h`)

## Data Flow

```
┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│  Input Data  │────▶│ Particle     │────▶│  Simulation  │
│  (nbody.dat) │     │ Initialization│     │    Loop      │
└──────────────┘     └──────────────┘     └──────────────┘
                                                  │
                     ┌────────────────────────────┘
                     ▼
        ┌────────────────────────┐
        │   Time Step Selection  │
        └────────────────────────┘
                     │
        ┌────────────┴────────────┐
        ▼                         ▼
┌──────────────────┐    ┌──────────────────┐
│ Irregular Update │    │  Regular Update  │
│ (neighbors)      │    │  (all particles) │
└──────────────────┘    └──────────────────┘
        │                         │
        ▼                         ▼
┌──────────────────┐    ┌──────────────────┐
│ Few-body Check   │    │ GPU Acceleration │
│ (SDAR groups)    │    │ (force calc)     │
└──────────────────┘    └──────────────────┘
                     │
                     ▼
        ┌────────────────────────┐
        │   HDF5 Output          │
        │   (Step_N groups)      │
        └────────────────────────┘
```

## Key Abstractions

| Abstraction | File | Purpose |
|-------------|------|---------|
| `Particle` | `src/particle.h` | Core particle data structure |
| `Group` | `src/FewBody/group.h` | Few-body subsystem (binaries, multiples) |
| `GlobalVariable` | `src/global_state.h` | Shared simulation state |
| `QueueScheduler` | `src/queue_scheduler.h` | Task distribution to workers |
| `Worker` | `src/worker.h` | Worker process state |

## Entry Points

| Entry Point | File | Description |
|-------------|------|-------------|
| `main()` | `src/main.cpp` | Program entry, MPI init, root/worker dispatch |
| `RootRoutines()` | `src/root_routines.cpp` | Main simulation loop (root process) |
| `WorkerRoutines()` | `src/worker_routines.cpp` | Worker process loop |

## MPI Communication

- **Shared Memory Windows**: `MPI_Win` for particle array, global state, neighbor lists
- **Custom MPI Types**: `queue_type_mpi`, `iparticle_type_mpi`, `jparticle_type_mpi`
- **Tags**: `TASK_TAG`, `PTCL_TAG`, `TIME_TAG`, `QUEUE_TAG`, `TERMINATE_TAG`
