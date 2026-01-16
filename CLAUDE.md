# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ABYSS is an N-body code for self-consistent gravitational dynamics of stars and black holes. It uses a 4th-order Hermite integrator with block time steps, splitting forces into regular (far) and irregular (near-neighbor) components. The code supports MPI parallelization with a root-worker architecture and optional CUDA GPU acceleration.

## Build Commands

```bash
# Use the build script (recommended)
./build.sh              # Standard CPU build
./build.sh --cuda       # With CUDA GPU acceleration
./build.sh --clean      # Clean before building

# Or manual build
cd src && make          # Standard build
cd src && make USE_CUDA=1    # With CUDA
cd src && make clean    # Clean
```

**Dependencies:**
- MPI compiler (mpicxx or mpiicpc)
- HDF5 with C++ bindings (set `HDF5_DIR` environment variable)
- toml11 library (header-only, in `toml11/include`)
- CUDA toolkit (optional, for GPU acceleration)
- SEVN library (optional, for stellar evolution)

**NYU Greene HPC:**
The `build.sh` script auto-detects paths on Greene. Otherwise:
```bash
export HDF5_DIR=/hpc/software/2025/spack/opt/spack/linux-rhel8-x86_64/gcc-12.4.0/hdf5-1.14.5-yke7pax3ot3kyspoqxw6lc4qolvpk43t
export PATH=/hpc/software/2025/spack/opt/spack/linux-rhel8-x86_64/gcc-12.4.0/openmpi-4.1.7-wn5igvodjfykt4bpgrf5n64iilky7udv/bin:$PATH
```

## Running

```bash
mpirun -np <N> ./abyss.exe -c config.toml
```

Config file (TOML format) parameters:
- `Filename` - initial conditions file (space-separated: x y z vx vy vz mass)
- `eta` - timestepping parameter
- `FixNumNeighbor` - number of neighbors (sqrt(N) recommended)
- `InitialRadius` - initial neighbor search radius (pc)
- `StopTime` - simulation duration (yr)
- `dtOutput` - output interval (yr)
- `OutputDirectory` - output directory name

## Architecture

### Core Components

**MPI Structure (main.cpp, RootRoutines.cpp, WorkerRoutines.cpp):**
- Rank 0 (Root): Manages simulation loop, time synchronization, particle output
- Ranks 1+ (Workers): Process irregular force calculations for assigned particles

**Time Integration:**
- Block timesteps with hierarchical levels (`TimeBlockIrr`, `TimeBlockReg`)
- Regular timestep: all particles updated at synchronized intervals
- Irregular timestep: neighbor-based updates at individual particle timescales
- `IrregularRoutines.cpp` / `RegularRoutines.cpp` implement the respective force calculations

**Particle Data (particle.h, Particle struct):**
- Position/Velocity in 3D
- Hermite acceleration derivatives up to 4th order (`a_tot`, `a_reg`, `a_irr`)
- Neighbor list management via offsets into shared `Neighbors` array
- Block-based time tracking (`CurrentBlockIrr`, `CurrentBlockReg`, `TimeBlockIrr`, `TimeBlockReg`)

### Few-Body Dynamics (src/FewBody/)

Uses SDAR (Slow-Down Algorithmic Regularization) from `SDAR/src/` for close encounters:
- `Group` struct wraps SDAR's `TimeTransformedSymplecticIntegrator`
- Particles in close binaries/multiples get center-of-mass particle (`isCMptcl`)
- Members tracked via `Members[]` array and `NumberOfMember`

### Key Data Structures

- `GlobalVariable` - shared state across MPI ranks (via MPI windows)
- `QueueScheduler` - task distribution for worker assignment
- `SkipList` - for efficient particle time ordering

### Compile-Time Options (defined in Makefile)

- `CUDA` - enable GPU acceleration
- `SEVN` - enable stellar evolution
- `FEWBODY` - enable few-body regularization (default on)
- `PERFORMANCETRACE` - enable timing instrumentation
- `MULTIMAP` - alternative regular timestep tracking

### Physical Units (def.h)

Code uses internal units where G=1:
- Position: 4 pc
- Time: 1e10 yr
- Velocity: 4e-10 pc/yr
- Mass: 0.0001424198 Msun

## I/O

- **Input:** Plain text, one particle per line (x y z vx vy vz mass in kpc, km/s, 1e-9 Msun)
- **Output:** HDF5 with timestep groups (`/Step_N`), containing arrays for PID, Mass, Position, Velocity
- Config: TOML format (using toml11 library)
